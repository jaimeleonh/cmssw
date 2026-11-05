#!/usr/bin/env python3
"""Run isolated cmsRun cases and retain FastTimer logs plus sampled GPU/RSS memory."""
import argparse
import ctypes as C
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import time

PACKAGE=Path(__file__).resolve().parents[2]
class Memory(C.Structure):
    _fields_=[('total',C.c_ulonglong),('free',C.c_ulonglong),('used',C.c_ulonglong)]

def parse_events(log):
    events=[]
    current=None
    in_modules=False
    for line in log.splitlines():
        if ' Event Summary ' in line and line.startswith('FastReport'):
            current={}; events.append(current); in_modules=True
        elif 'Run Summary' in line or 'Job Summary' in line:
            current=None; in_modules=False
        elif 'Processes and Paths' in line:
            in_modules=False
        elif current is not None and in_modules:
            m=re.match(r'FastReport\s+([\d.]+) ms\s+([\d.]+) ms\s+[+-]\d+ kB\s+[+-]\d+ kB\s+(.+)$',line)
            if m: current[m[3].strip()]={'cpu_ms':float(m[1]),'real_ms':float(m[2])}
    return events

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--engine',choices=['graph','compiled','torch','torch_original','none'],required=True)
    p.add_argument('--model',choices=['deepsets','part'],default='deepsets')
    p.add_argument('--batch',type=int,default=0)
    p.add_argument('--graph')
    p.add_argument('--events',type=int,default=100)
    p.add_argument('--threads',type=int,default=1)
    p.add_argument('--streams',type=int,default=1)
    p.add_argument('--sync',action='store_true')
    p.add_argument('--validate',action='store_true')
    p.add_argument('--tf32',action='store_true')
    p.add_argument('--observer',action='store_true')
    p.add_argument('--force-ieee',action='store_true',help='Set NVIDIA_TF32_OVERRIDE=0 for the entire cmsRun process')
    p.add_argument('--timeout',type=float,default=1200)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); prefix=a.output.resolve(); prefix.parent.mkdir(parents=True,exist_ok=True)
    env=dict(os.environ,CUDA_VISIBLE_DEVICES='0',BENCH_ENGINE=a.engine,BENCH_MODEL=a.model,
             BENCH_OUTPUT=str(prefix),BENCH_SYNC=str(int(a.sync)),BENCH_VALIDATE=str(int(a.validate)),
             BENCH_TF32=str(int(a.tf32)),BENCH_OBSERVER=str(int(a.observer)))
    if a.graph:env['BENCH_GRAPH']=a.graph
    if a.force_ieee:env['NVIDIA_TF32_OVERRIDE']='0'
    command=['cmsRun',str(PACKAGE/'test/cmssw_benchmark/run_cfg.py'),
        'buBaseDir=/mnt/ramdisk/jleonhol/raw','fuBaseDir=/tmp/particle_transformer_cmssw_fu',
        'buNumStreams=2','runNumber=41',f'maxEvents={a.events}','run=sc4Alpaka','backend=cuda_async',
        f'numThreads={a.threads}',f'numFwkStreams={a.streams}',f'batchSize={a.batch}']
    nvml=C.CDLL('libnvidia-ml.so.1'); assert nvml.nvmlInit_v2()==0
    handle=C.c_void_p(); assert nvml.nvmlDeviceGetHandleByIndex_v2(0,C.byref(handle))==0
    mem=Memory(); nvml.nvmlDeviceGetMemoryInfo(handle,C.byref(mem)); baseline=mem.used
    samples=[]; start=time.monotonic()
    with prefix.with_suffix('.log').open('w') as out:
        proc=subprocess.Popen(command,env=env,stdout=out,stderr=subprocess.STDOUT)
        while proc.poll() is None:
            if time.monotonic()-start>a.timeout:
                proc.terminate()
                try:proc.wait(timeout=10)
                except subprocess.TimeoutExpired:proc.kill();proc.wait()
                break
            nvml.nvmlDeviceGetMemoryInfo(handle,C.byref(mem))
            rss=0
            try:
                text=Path(f'/proc/{proc.pid}/status').read_text()
                m=re.search(r'^VmRSS:\s+(\d+)',text,re.M)
                if m:rss=int(m[1])*1024
            except FileNotFoundError:pass
            samples.append([round(time.monotonic()-start,4),mem.used,rss]); time.sleep(.02)
    elapsed=time.monotonic()-start
    log=prefix.with_suffix('.log').read_text(); event=parse_events(log)
    steady=event[5:]
    labels=sorted(set().union(*(x.keys() for x in steady))) if steady else []
    stats={}
    for label in labels:
        v=[x[label]['real_ms'] for x in steady if label in x]
        stats[label]={'mean_ms':statistics.mean(v),'median_ms':statistics.median(v),
                      'min_ms':min(v),'max_ms':max(v),'events':len(v)}
    data={'settings':vars(a)|{'output':str(prefix)},'command':command,'nvidia_tf32_override':env.get('NVIDIA_TF32_OVERRIDE'),'returncode':proc.returncode,
          'wall_seconds':elapsed,'events':event,'steady_skip':5,'steady':stats,
          'gpu_idle_bytes':baseline,'gpu_peak_sampled_bytes':max((x[1] for x in samples),default=0),
          'rss_peak_sampled_bytes':max((x[2] for x in samples),default=0),
          'memory_sampling_interval_seconds':.02,'memory_samples':samples}
    data['settings']['output']=str(prefix)
    completion=Path(str(prefix)+'.completion.jsonl')
    if a.observer and completion.exists():
        complete=sorted((json.loads(line) for line in completion.read_text().splitlines()),key=lambda x:x['ns'])
        warmup=max(5,5*a.streams)
        selected=complete[warmup-1:]
        if len(selected)>1:
            seconds=(selected[-1]['ns']-selected[0]['ns'])/1e9
            data['completion_throughput']={'warmup_events':warmup,'intervals':len(selected)-1,
                'seconds':seconds,'events_per_second':(len(selected)-1)/seconds,
                'ms_per_event':seconds*1000/(len(selected)-1)}

    timer=Path(str(prefix)+'.timer.json')
    if timer.exists(): data['timer_summary']=json.loads(timer.read_text())
    prefix.with_suffix('.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({k:v for k,v in data.items() if k not in ('events','memory_samples','command')},indent=2),flush=True)
    if proc.returncode:raise SystemExit(proc.returncode)
    count=len(complete) if a.observer and completion.exists() else len(event)
    if count!=a.events:raise RuntimeError(f'Expected {a.events} completed events, found {count}')
if __name__=='__main__':main()
