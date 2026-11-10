#!/usr/bin/env python3
"""Run isolated, TTbar PUPPI throughput jobs after cmsenv."""
import argparse
import csv
import json
import os
from pathlib import Path
import re
import subprocess
import time

# Plot labels are cases; the GPU producer has only three algorithm names.
CASES = {
    'AK4': ('ak4', 'SeededCone', 16, 'scPhase2AK4PFDemo'),
    'SC4_16': ('sc4Alpaka', 'SeededCone', 16, 'scPhase2SC4PFAlpaka'),
    'SC4_unconstrained': ('sc4Alpaka', 'SeededCone', 1024, 'scPhase2SC4PFAlpaka'),
    'SCNMS4': ('sc4Alpaka', 'SCNMSW', 0, 'scPhase2SC4PFAlpaka'),
    'SCNMSW4': ('sc4Alpaka', 'SCNMSW', 0, 'scPhase2SC4PFAlpaka'),
    'LinkTree': ('sc4Alpaka', 'LinkTree', 0, 'scPhase2SC4PFAlpaka'),
    'SCNMSW': ('sc4Alpaka', 'SCNMSW', 0, 'scPhase2SC4PFAlpaka'),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', default=f'/mnt/ramdisk/{os.environ["USER"]}/raw')
    parser.add_argument('--output-dir', default='throughput-results')
    parser.add_argument('--orbits', type=int, default=1000, help='Number of orbits; -1 reads the complete input file')
    parser.add_argument('--synchronize', choices=('True','False'), default='True',
                        help='Wait at Alpaka module boundaries for accurate FastTimerService module times')
    parser.add_argument('--max-rounds', type=int, default=1)
    parser.add_argument('--alpha-seed', type=float, default=None, help='SCNMSW exponent; default 0 for SCNMS4, otherwise 2')
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--streams', type=int, default=1)
    parser.add_argument('--cpu', type=int, help='Pin cmsRun to this CPU (single-thread runs only)')
    parser.add_argument('--gpu', default='0', help='CUDA_VISIBLE_DEVICES value')
    parser.add_argument('--expect-gpu', default='L4', help='Refuse a different GPU; use L4 for local validation')
    parser.add_argument('--backend', choices=('cuda_async', 'serial_sync'), default='cuda_async')
    parser.add_argument('--cases', nargs='+', choices=list(CASES), default=list(CASES)[:6])
    args = parser.parse_args()
    if min(args.repeats, args.threads, args.streams) <= 0 or args.orbits == 0 or args.orbits < -1:
        parser.error('repeats, threads and streams must be positive; orbits must be positive or -1 for all')
    if args.cpu is not None and args.threads != 1:
        parser.error('--cpu is intended for single-thread measurements')
    raw = Path(args.input_dir).resolve()/'run000039/run000039_ls0001_index000001_stream00.raw'
    if not raw.is_file():
        parser.error(f'input missing: {raw}')
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, CUDA_VISIBLE_DEVICES=args.gpu)
    inventory = ''
    if args.backend == 'cuda_async':
        inventory = subprocess.check_output(['nvidia-smi', '--query-gpu=index,name,uuid,driver_version', '--format=csv,noheader'], text=True)
        # For an index or UUID, verify the selected card rather than any card on the host.
        selected = subprocess.check_output(['nvidia-smi', '-i', args.gpu, '--query-gpu=name', '--format=csv,noheader'], text=True).strip()
        if args.expect_gpu.lower() not in selected.lower():
            parser.error(f'selected GPU is {selected}, expected {args.expect_gpu}')
    (output/'hardware.txt').write_text(inventory + subprocess.check_output(['lscpu'], text=True))
    config = Path(__file__).resolve().with_name('runScoutingPhase2PF_cfg.py')
    (output/'study.json').write_text(json.dumps(vars(args), indent=2)+'\n')
    rows = []
    for case in args.cases:
        run, algo, njets, jet_label = CASES[case]
        alpha = args.alpha_seed if args.alpha_seed is not None else (0.0 if case == 'SCNMS4' else 2.0)
        for repeat in range(1, args.repeats + 1):
            jobdir = output/f'{case}-{repeat}'
            jobdir.mkdir()  # refuse to overwrite an existing result
            command = ['cmsRun', str(config), f'run={run}', f'jetAlgo={algo}', f'njets={njets}',
                       f'backend={args.backend}', f'synchronizeGPU={args.synchronize}', 'runNumber=39', f'buBaseDir={raw.parent.parent}',
                       f'fuBaseDir={jobdir / "fu"}', 'buNumStreams=1', 'timeslices=1', 'tmuxPeriod=1',
                       f'numThreads={args.threads}', f'numFwkStreams={args.streams}', f'maxEvents={args.orbits}',
                       'jetR=0.4', 'RSeed=0.4', 'RCen=0.4', 'RClu=0.4', 'RLink=0.36',
                       f'alphaSeed={alpha}', f'maxRounds={args.max_rounds}', 'minSeedPt=0', 'dumpClusters=False' ]
            if args.cpu is not None:
                command = ['taskset', '-c', str(args.cpu)] + command
            (jobdir/'command.json').write_text(json.dumps(command, indent=2)+'\n')
            print(f'{case} repeat {repeat}: {jobdir}', flush=True)
            start = time.monotonic()
            with (jobdir/'cmsRun.log').open('w') as log:
                subprocess.run(command, cwd=jobdir, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
            wall = time.monotonic() - start
            logtext = (jobdir/'cmsRun.log').read_text()
            if args.backend == 'cuda_async' and run == 'sc4Alpaka' and 'AlpakaServiceCudaAsync succesfully initialised' not in logtext:
                raise RuntimeError(f'{case}: CUDA backend was not initialized')
            summary = json.loads(next(jobdir.glob('resources.*.json')).read_text())
            modules = {module['label']: module for module in summary['modules']}
            jet = modules[jet_label]
            events = summary['total']['events']
            counts = re.findall(r'Processed (\d+) orbits, of which (\d+) truncated, and (\d+) events\.', logtext)
            if not counts:
                raise RuntimeError(f'{case}: no unpacker BX count in log')
            orbits, truncated, bx = map(int, counts[-1])
            if events <= 0 or orbits != events or (args.orbits != -1 and events != args.orbits) or truncated or bx != events * 3564 or jet['events'] != events:
                raise RuntimeError(f'{case}: incomplete/filtered input: events={events}, unpacker={counts[-1]}, jet={jet["events"]}')
            # FastTimerService JSON reports cumulative times in milliseconds.
            unpack_label = 'scPhase2PFRawToDigiStruct' if run == 'ak4' else 'scPhase2PFRawToDigiAlpaka'
            jet_ms = jet['time_real'] / jet['events']
            path_ms = sum(modules[label]['time_real'] / modules[label]['events'] for label in (unpack_label, 'goodOrbitsByNBX', jet_label))
            event_loop_s = float(re.search(r'event loop Real/event\s*=\s*([\d.eE+-]+)', logtext).group(1))
            rows.append(dict(case=case, repeat=repeat, algorithm=algo, max_rounds=args.max_rounds, alpha_seed=alpha, synchronize=args.synchronize, backend='cpu' if run == 'ak4' else args.backend, threads=args.threads, streams=args.streams, cpu=args.cpu if args.cpu is not None else 'unrestricted', orbits=events, bx=bx,
                             raw_to_jet_bx_s=3564*1000/path_ms, jet_module_bx_s=3564*1000/jet_ms,
                             event_loop_bx_s=3564/event_loop_s, event_loop_ms_per_orbit=event_loop_s*1000, process_wall_s=wall,
                             jet_ms_per_orbit=jet_ms, raw_to_jet_ms_per_orbit=path_ms))
            with (output/'throughput.csv').open('w', newline='') as dest:
                writer = csv.DictWriter(dest, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            print(f'  Event loop {rows[-1]["event_loop_bx_s"]:.0f} BX/s; RAW-to-jet {rows[-1]["raw_to_jet_bx_s"]:.0f} BX/s; jet module {rows[-1]["jet_module_bx_s"]:.0f} BX/s', flush=True)


if __name__ == '__main__':
    main()
