#pragma once
// Minimal CUDA stand-in for the parts of the Alpaka API the interpreter uses.
// The real alpaka CUDA backend maps every one of these one-to-one, so a kernel
// compiled against this shim issues the same instructions as the CMSSW build:
// syncBlockThreads is __syncthreads(), the math calls are the CUDA intrinsics,
// and the index queries are the built-in variables.  It exists so the kernel
// can be benchmarked and profiled without a CMSSW link.
#include <cstdint>
#include <cuda_runtime.h>

#define ALPAKA_FN_ACC __host__ __device__  // host too, so the host-side weight transpose can call the generated accessors
#define ALPAKA_FN_HOST_ACC __host__ __device__

namespace alpaka {
  struct Block {};
  struct Threads {};
  struct Grid {};
  struct Blocks {};

  struct Index1 {
    unsigned int value;
    __device__ unsigned int operator[](unsigned int) const { return value; }
  };

  template <typename TScope, typename TUnit, typename TAcc>
  __device__ inline Index1 getIdx(TAcc const&) {
    if constexpr (std::is_same_v<TScope, Block>)
      return Index1{threadIdx.x};
    else
      return Index1{blockIdx.x};
  }

  template <typename TScope, typename TUnit, typename TAcc>
  __device__ inline Index1 getWorkDiv(TAcc const&) {
    if constexpr (std::is_same_v<TScope, Block>)
      return Index1{blockDim.x};
    else
      return Index1{gridDim.x};
  }

  template <typename TAcc>
  __device__ inline void syncBlockThreads(TAcc const&) {
    __syncthreads();
  }
}  // namespace alpaka

namespace alpaka::math {
  template <class A> __device__ inline float sqrt(A const&, float x) { return ::sqrtf(x); }
  template <class A> __device__ inline float log(A const&, float x) { return ::logf(x); }
  template <class A> __device__ inline float exp(A const&, float x) { return ::expf(x); }
  template <class A> __device__ inline float erf(A const&, float x) { return ::erff(x); }
  template <class A> __device__ inline float tanh(A const&, float x) { return ::tanhf(x); }
  template <class A> __device__ inline float atan2(A const&, float y, float x) { return ::atan2f(y, x); }
  template <class A> __device__ inline float max(A const&, float x, float y) { return ::fmaxf(x, y); }
  template <class A> __device__ inline float min(A const&, float x, float y) { return ::fminf(x, y); }
  template <class A> __device__ inline float abs(A const&, float x) { return ::fabsf(x); }
}  // namespace alpaka::math
