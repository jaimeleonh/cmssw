#!/usr/bin/env bash
# Run inside the package after `eval "$(scram runtime -sh)"`.
set -euo pipefail
benchmark_dir=${1:-/tmp}
mkdir -p "$benchmark_dir"
torch_path=$(python3 -c 'import torch; print(torch.__path__[0])')
torch_abi=$(python3 -c 'import torch; print(int(torch._C._GLIBCXX_USE_CXX11_ABI))')
nvcc -O3 -std=c++20 -arch=sm_89 -Xcompiler -fPIC -shared \
  -Itest/cudastubs -I../.. test/benchmark_runtime_bridge.cu \
  -o "$benchmark_dir/particle_transformer_benchmark_graph.so"
c++ -O3 -std=c++17 -fPIC -shared -D_GLIBCXX_USE_CXX11_ABI="$torch_abi" \
  -I"$torch_path/include" -I"$torch_path/include/torch/csrc/api/include" \
  test/benchmark_torch_bridge.cc -L"$torch_path/lib" -Wl,-rpath,"$torch_path/lib" \
  -ltorch -ltorch_cpu -lc10 -o "$benchmark_dir/particle_transformer_benchmark_torch.so"
