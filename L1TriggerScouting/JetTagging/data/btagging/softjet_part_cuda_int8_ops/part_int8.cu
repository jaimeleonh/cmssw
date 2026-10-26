
#include <ATen/ATen.h>
#include <c10/cuda/CUDAGuard.h>
#include <c10/cuda/CUDAException.h>
#include <c10/cuda/CUDAStream.h>
#include <cuda_runtime.h>
#include <algorithm>
#include <cstdint>
#include <optional>

void check_inputs(const at::Tensor&, const at::Tensor&, const at::Tensor&,
                  const std::optional<at::Tensor>&);

// One block computes [16 rows, 32 outputs], with K tiles of 32. Each thread
// computes two outputs. Weights [K,N] are loaded coalesced and converted in
// registers. There is never a full dequantized weight allocation.
__global__ void linear_kernel(const float* __restrict__ x,
    const int8_t* __restrict__ w, const float* __restrict__ s,
    const float* __restrict__ bias, float* __restrict__ y,
    int64_t rows, int64_t outputs, int64_t inputs, int64_t row_tiles) {
  __shared__ float a[16][32];
  __shared__ int8_t q[32][32];
  const int tx = threadIdx.x, ty = threadIdx.y, tid = ty * 32 + tx;
  const int64_t tiles = row_tiles * ((outputs + 31) / 32);
  // A grid-stride tile loop avoids CUDA grid.y's 65535 limit at large batches.
  for (int64_t tile = blockIdx.x; tile < tiles; tile += gridDim.x) {
    const int64_t row0 = (tile % row_tiles) * 16;
    const int64_t col0 = (tile / row_tiles) * 32;
    const int64_t col = col0 + tx;
    float v0 = 0.f, v1 = 0.f;
    for (int64_t k0 = 0; k0 < inputs; k0 += 32) {
      for (int i = tid; i < 16 * 32; i += 256) {
        const int r = i / 32, k = i % 32;
        a[r][k] = row0 + r < rows && k0 + k < inputs
                    ? x[(row0 + r) * inputs + k0 + k] : 0.f;
      }
      for (int i = tid; i < 32 * 32; i += 256) {
        const int k = i / 32, c = i % 32;
        q[k][c] = k0 + k < inputs && col0 + c < outputs
                    ? w[(k0 + k) * outputs + col0 + c] : int8_t(0);
      }
      __syncthreads();
      #pragma unroll
      for (int k = 0; k < 32; ++k) {
        const float weight = static_cast<float>(q[k][tx]);
        v0 = fmaf(a[ty][k], weight, v0);
        v1 = fmaf(a[ty + 8][k], weight, v1);
      }
      __syncthreads();
    }
    if (col < outputs) {
      const float scale = s[col], b = bias ? bias[col] : 0.f;
      if (row0 + ty < rows)
        y[(row0 + ty) * outputs + col] = fmaf(v0, scale, b);
      if (row0 + ty + 8 < rows)
        y[(row0 + ty + 8) * outputs + col] = fmaf(v1, scale, b);
    }
    __syncthreads();
  }
}

at::Tensor linear_cuda(const at::Tensor& x, const at::Tensor& w,
    const at::Tensor& s, const std::optional<at::Tensor>& b) {
  check_inputs(x, w, s, b);
  const c10::cuda::CUDAGuard guard(x.device());
  auto xc = x.contiguous(), wc = w.contiguous(), sc = s.contiguous();
  auto bc = b.has_value() ? b->contiguous() : at::Tensor();
  auto shape = x.sizes().vec();
  shape.back() = w.size(1);
  auto result = at::empty(shape, x.options());
  const int64_t rows = x.numel() / w.size(0);
  if (rows == 0) return result;
  const int64_t row_tiles = (rows + 15) / 16;
  const int64_t tiles = row_tiles * ((w.size(1) + 31) / 32);
  const int blocks = static_cast<int>(std::min<int64_t>(tiles, 65535));
  linear_kernel<<<blocks, dim3(32, 8), 0,
      c10::cuda::getCurrentCUDAStream(x.get_device()).stream()>>>(
      xc.data_ptr<float>(), wc.data_ptr<int8_t>(), sc.data_ptr<float>(),
      bc.defined() ? bc.data_ptr<float>() : nullptr, result.data_ptr<float>(),
      rows, w.size(1), w.size(0), row_tiles);
  C10_CUDA_KERNEL_LAUNCH_CHECK();
  return result;
}
