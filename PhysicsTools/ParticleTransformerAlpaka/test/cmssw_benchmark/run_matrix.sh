#!/usr/bin/env bash
set -euo pipefail
# Run in PhysicsTools/ParticleTransformerAlpaka after cmsenv and scram b.
result_dir=${1:-test/performance_results/cmssw}
runner=test/cmssw_benchmark/run.py
run_case() {
    local name=$1
    shift
    if [[ -f "$result_dir/$name.json" ]]; then
        echo "Already recorded: $name"
        return
    fi
    echo "Running: $name"
    python3 "$runner" --output "$result_dir/$name" "$@" > "$result_dir/$name.driver.txt" 2>&1
    python3 - "$result_dir/$name.json" <<'PY'
import sys,json
x=json.load(open(sys.argv[1])); print(x['settings']['engine'],x['settings']['batch'],
  'ms/event:',x['steady'].get('total'), 'tagger:',x['steady'].get('SoftJetIdSC4'),flush=True)
PY
}
run_case baseline_sync --engine none --events 100 --sync
run_case deepsets_graph_sync --engine graph --model deepsets --events 100 --sync
for batch in 128 512 2048 4096 8192 65536; do
    run_case deepsets_torch_b${batch}_sync --engine torch --model deepsets --batch "$batch" --events 100 --sync
done
run_case part_graph_fp32_sync --engine graph --model part --events 100 --sync
run_case part_graph_int8_sync --engine graph --model part --graph /tmp/benchmark_part_int8.nnm --events 100 --sync
run_case part_graph_qat_sync --engine graph --model part --graph PhysicsTools/ParticleTransformerAlpaka/data/b_kinadd.nnm --events 100 --sync
run_case part_compiled_qat_sync --engine compiled --model part --events 100 --sync
for batch in 128 256 512 1024 4096; do
    run_case part_torch_b${batch}_sync_tune --engine torch --model part --batch "$batch" --events 20 --sync
done
run_case part_torch_full_sync_tune --engine torch --model part --batch 65536 --events 10 --sync
run_case validate_part_graph --engine graph --model part --events 10 --sync --validate
run_case validate_part_int8 --engine graph --model part --graph /tmp/benchmark_part_int8.nnm --events 10 --sync --validate
run_case validate_part_torch --engine torch --model part --batch 256 --events 10 --sync --validate
run_case validate_qat_graph --engine graph --model part --graph PhysicsTools/ParticleTransformerAlpaka/data/b_kinadd.nnm --events 10 --sync --validate
run_case validate_qat_compiled --engine compiled --model part --events 10 --sync --validate
