#!/usr/bin/env bash
set -euo pipefail
result_dir=${1:-test/performance_results/cmssw}
runner=test/cmssw_benchmark/run.py
best_batch() {
    python3 - "$result_dir" "$1" <<'PY'
import sys,json,pathlib
files=pathlib.Path(sys.argv[1]).glob(sys.argv[2]+'_torch_b*_sync*.json')
rows=[json.load(p.open()) for p in files]
rows=[x for x in rows if x.get('returncode')==0 and x.get('steady',{}).get('total')]
print(min(rows,key=lambda x:x['steady']['total']['mean_ms'])['settings']['batch'])
PY
}
run_case() {
    local name=$1
    shift
    if [[ -f "$result_dir/$name.json" ]]; then
        echo "Already recorded: $name"; return
    fi
    echo "Running: $name"
    python3 "$runner" --observer --output "$result_dir/$name" "$@" > "$result_dir/$name.driver.txt" 2>&1
    python3 - "$result_dir/$name.json" <<'PY'
import json,sys
x=json.load(open(sys.argv[1]));print(x['settings']['engine'],x['settings']['batch'],
    'completed-event:',x.get('completion_throughput'),'GPU MiB',x['gpu_peak_sampled_bytes']/2**20,flush=True)
PY
}
deep_batch=$(best_batch deepsets)
part_batch=$(best_batch part)
run_case baseline_async --engine none --events 100
for streams in 1 4 8; do
    run_case deepsets_graph_async_s${streams} --engine graph --model deepsets --events 100 --threads "$streams" --streams "$streams"
    run_case deepsets_torch_async_s${streams} --engine torch --model deepsets --batch "$deep_batch" --events 100 --threads "$streams" --streams "$streams"
done
run_case part_graph_async_s1 --engine graph --model part --events 30
run_case part_torch_async_s1 --engine torch --model part --batch "$part_batch" --events 30
run_case part_graph_async_s4 --engine graph --model part --events 60 --threads 4 --streams 4
run_case part_torch_async_s4 --engine torch --model part --batch "$part_batch" --events 60 --threads 4 --streams 4
run_case part_int8_async_s1 --engine graph --model part --graph /tmp/benchmark_part_int8.nnm --events 30
run_case deepsets_torch_force_ieee_async_s1 --engine torch --model deepsets --batch "$deep_batch" --events 100 --force-ieee
run_case part_torch_force_ieee_async_s1 --engine torch --model part --batch "$part_batch" --events 30 --force-ieee
run_case part_torch_tf32_async_s1 --engine torch --model part --batch "$part_batch" --events 30 --tf32
run_case compiled_qat_async_s1 --engine compiled --model part --events 30
run_case graph_qat_async_s1 --engine graph --model part --graph PhysicsTools/ParticleTransformerAlpaka/data/b_kinadd.nnm --events 30
# DeepSets differences are small: repeat the throughput plateau configurations.
for repeat in 2 3; do
    run_case deepsets_graph_async_s4_repeat${repeat} --engine graph --model deepsets --events 100 --threads 4 --streams 4
    run_case deepsets_torch_async_s1_repeat${repeat} --engine torch --model deepsets --batch "$deep_batch" --events 100
 done
