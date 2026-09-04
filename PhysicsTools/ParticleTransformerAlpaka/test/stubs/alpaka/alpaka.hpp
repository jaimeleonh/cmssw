#pragma once
// Minimal host-side stand-in for the parts of the Alpaka API used by the
// ParticleTransformer kernels.  It exists only so that the generated graph can
// be exercised, and its block simulation compared against the scalar
// reference, without a CMSSW or CUDA installation.
#include <array>
#include <cstdint>
#include <cmath>
#define ALPAKA_FN_ACC
#define ALPAKA_FN_HOST_ACC
namespace alpaka {
  struct Block {};
  struct Threads {};
  template <typename TScope, typename TUnit, typename TAcc>
  auto getIdx(TAcc const& acc) {
    return std::array<uint32_t, 1>{acc.lane};
  }
  template <typename TScope, typename TUnit, typename TAcc>
  auto getWorkDiv(TAcc const& acc) {
    return std::array<uint32_t, 1>{acc.lanes};
  }
  template <typename TAcc>
  void syncBlockThreads(TAcc const& acc) {
    acc.sync();
  }
  template <typename T, int TUniqueId, typename TAcc>
  T& declareSharedVar(TAcc const& acc) {
    return *static_cast<T*>(acc.sharedBlock);
  }
}  // namespace alpaka
namespace alpaka::math {
  template <class A> float sqrt(A const&, float x) { return std::sqrt(x); }
  template <class A> float log(A const&, float x) { return std::log(x); }
  template <class A> float exp(A const&, float x) { return std::exp(x); }
  template <class A> float erf(A const&, float x) { return std::erf(x); }
  template <class A> float atan2(A const&, float y, float x) { return std::atan2(y, x); }
  template <class A> float max(A const&, float x, float y) { return x > y ? x : y; }
  template <class A> float min(A const&, float x, float y) { return x < y ? x : y; }
  template <class A> float abs(A const&, float x) { return std::fabs(x); }
}  // namespace alpaka::math
