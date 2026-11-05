// Standalone regression test for the DeepSets fast kernel.
//
// It runs the same jets through a scalar reference transcription of the
// exported INT8 DeepSets graph and through the block-cooperative production
// kernel (simulated with real host threads and a std::barrier), and reports
// the largest disagreement.  The two paths share no code beyond the generated
// tensor descriptors; the reference evaluates every layer on exactly one
// thread, the kernel with the full block.
//
// Run with the minimal Alpaka stub:
//   g++ -std=c++20 -O2 -I src -I src/PhysicsTools/ParticleTransformerAlpaka/test/stubs \
//       -o deepsets_kernel_test \
//       src/PhysicsTools/ParticleTransformerAlpaka/test/compile_deepsets_kernel.cc -pthread
#include <algorithm>
#include <array>
#include <barrier>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <thread>
#include <vector>

#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/DeepSetsFastKernel.h"

struct Acc {
  uint32_t lane = 0;
  uint32_t lanes = 1;
  std::barrier<>* blockBarrier = nullptr;
  void* sharedBlock = nullptr;

  void sync() const {
    if (blockBarrier != nullptr)
      blockBarrier->arrive_and_wait();
  }
};

namespace deepsets = part::deepsets;

template <typename T>
std::vector<T> readBlob(char const* path, std::size_t expected) {
  std::ifstream input(path, std::ios::binary | std::ios::ate);
  if (!input || static_cast<std::size_t>(input.tellg()) != expected * sizeof(T))
    throw std::runtime_error(std::string("model file has the wrong size: ") + path);
  input.seekg(0);
  std::vector<T> values(expected);
  input.read(reinterpret_cast<char*>(values.data()), static_cast<std::streamsize>(values.size() * sizeof(T)));
  if (!input)
    throw std::runtime_error(std::string("model file read failed: ") + path);
  return values;
}

// ---------------------------------------------------------------- scalar reference

void batchNorm(part::generated::BatchNorm const& p, float* x, int rows) {
  for (int i = 0; i < rows * p.size; ++i) {
    int const c = i % p.size;
    x[i] = (x[i] - p.mean[c]) / std::sqrt(p.var[c] + p.eps) * p.gamma[c] + p.beta[c];
  }
}

void runReference(part::generated::ModelView model,
                  float const* features,
                  bool const* mask,
                  float* probabilities) {
  constexpr int N = deepsets::kMaxParticles;
  constexpr int F = part::generated::input_features;
  constexpr int W = part::generated::embed_width;
  float x[N * W], quantScratch[W], state[W];

  auto const input = part::generated::inputLinear(model);
  for (int t = 0; t < N; ++t)
    for (int o = 0; o < W; ++o) {
      float value = input.bias[o];
      for (int k = 0; k < F; ++k)
        value += features[t * F + k] * input.weight[o * F + k];
      x[t * W + o] = value;
    }
  batchNorm(part::generated::phiNorm(model, 0), x, N);
  for (int t = 0; t < N; ++t)
    for (int c = 0; c < W; ++c)
      x[t * W + c] = x[t * W + c] > 0.f ? x[t * W + c] : 0.f;

  auto quantLinear = [&](part::generated::Linear const& l, float const* in, float* out) {
    float const sx = l.activation[0], inv = 1.f / sx;
    for (int o = 0; o < l.out; ++o) {
      int acc = 0;
      for (int k = 0; k < l.in; ++k)
        acc += int(part::quantizeOne(in[k], inv)) * int(l.weight[o * l.in + k]);
      out[o] = float(acc) * (sx * l.scale[o]) + l.bias[o];
    }
  };

  auto const phi = part::generated::phiProject(model);
  for (int t = 0; t < N; ++t) {
    quantLinear(phi, x + t * W, quantScratch);
    for (int c = 0; c < W; ++c)
      x[t * W + c] = quantScratch[c];
  }
  batchNorm(part::generated::phiNorm(model, 1), x, N);
  for (int t = 0; t < N; ++t)
    for (int c = 0; c < W; ++c)
      x[t * W + c] = mask[t] && x[t * W + c] > 0.f ? x[t * W + c] : 0.f;

  for (int c = 0; c < W; ++c) {
    float value = 0.f;
    for (int t = 0; t < N; ++t)
      value += x[t * W + c];
    state[c] = value;
  }

  quantLinear(part::generated::rhoProject(model), state, quantScratch);
  for (int c = 0; c < W; ++c)
    state[c] = quantScratch[c];
  batchNorm(part::generated::rhoNorm(model), state, 1);
  for (int c = 0; c < W; ++c)
    state[c] = state[c] > 0.f ? state[c] : 0.f;

  auto const head = part::generated::classifier(model);
  float logits[2];
  for (int o = 0; o < 2; ++o) {
    float value = head.bias[o];
    for (int k = 0; k < W; ++k)
      value += state[k] * head.weight[o * W + k];
    logits[o] = value;
  }
  float const maximum = std::max(logits[0], logits[1]);
  float const e0 = std::exp(logits[0] - maximum), e1 = std::exp(logits[1] - maximum);
  probabilities[0] = e0 / (e0 + e1);
  probabilities[1] = e1 / (e0 + e1);
}

// ---------------------------------------------------------------- fast kernel

void runFast(part::generated::ModelView model,
             float const* features,
             bool const* mask,
             float* probabilities,
             uint32_t blockThreads) {
  auto shared = std::make_unique<deepsets::DeepSetsShared>();
  std::memset(shared.get(), 0, sizeof(deepsets::DeepSetsShared));
  std::copy_n(features, deepsets::kMaxParticles * part::generated::input_features, shared->features);
  for (int t = 0; t < deepsets::kMaxParticles; ++t)
    shared->mask[t] = mask[t] ? 1.f : 0.f;

  std::barrier blockBarrier(static_cast<std::ptrdiff_t>(blockThreads));
  {
    std::vector<std::jthread> threads;
    threads.reserve(blockThreads);
    for (uint32_t lane = 0; lane < blockThreads; ++lane)
      threads.emplace_back([&, lane] { deepsets::runDeepSets(Acc{lane, blockThreads, &blockBarrier, nullptr}, model, *shared); });
  }
  probabilities[0] = shared->logits[0];
  probabilities[1] = shared->logits[1];
}

// ---------------------------------------------------------------- warp layout

// One warp of real threads evaluates one jet through a WarpAcc, with one
// std::barrier per warp, so warps never wait for each other -- exactly the
// independence the __syncwarp path relies on on the GPU.
void runWarp(part::generated::ModelView model,
             float const* const* features,   // [warps][N * F]
             bool const* const* mask,        // [warps][N]
             float* probabilities,           // [warps * 2]
             uint32_t warps) {
  std::vector<std::unique_ptr<deepsets::DeepSetsShared>> states(warps);
  std::vector<std::unique_ptr<std::barrier<>>> barriers(warps);
  for (uint32_t w = 0; w < warps; ++w) {
    states[w] = std::make_unique<deepsets::DeepSetsShared>();
    std::memset(states[w].get(), 0, sizeof(deepsets::DeepSetsShared));
    std::copy_n(features[w], deepsets::kMaxParticles * part::generated::input_features, states[w]->features);
    for (int t = 0; t < deepsets::kMaxParticles; ++t)
      states[w]->mask[t] = mask[w][t] ? 1.f : 0.f;
    barriers[w] = std::make_unique<std::barrier<>>(static_cast<std::ptrdiff_t>(deepsets::kWarpLanes));
  }
  {
    std::vector<std::jthread> threads;
    threads.reserve(warps * deepsets::kWarpLanes);
    for (uint32_t w = 0; w < warps; ++w)
      for (int lane = 0; lane < deepsets::kWarpLanes; ++lane)
        threads.emplace_back([&, w, lane] {
          Acc acc{0, 1, barriers[w].get(), nullptr};  // sync() covers this warp's lanes
          deepsets::runDeepSets(deepsets::WarpAcc<Acc>{acc, lane}, model, *states[w]);
        });
  }
  for (uint32_t w = 0; w < warps; ++w) {
    probabilities[2 * w] = states[w]->logits[0];
    probabilities[2 * w + 1] = states[w]->logits[1];
  }
}

// ---------------------------------------------------------------- driver

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "usage: compile_deepsets_kernel WEIGHTS PARAMS\n";
    return 2;
  }
  auto weights = readBlob<std::int8_t>(argv[1], part::generated::weights_size);
  auto params = readBlob<float>(argv[2], part::generated::params_size);
  part::generated::ModelView model{weights.data(), params.data()};

  std::cout << "shared bytes per jet: " << sizeof(deepsets::DeepSetsShared) << "\n";

  constexpr int N = deepsets::kMaxParticles;
  constexpr int F = part::generated::input_features;
  std::mt19937 rng(3407);
  std::normal_distribution<float> normal(0.f, .6f);
  constexpr uint32_t threadings[] = {1, 32, 96, 256};

  double worst = 0.;
  int jets = 0;
  for (int multiplicity = 1; multiplicity <= N; ++multiplicity) {
    for (int hole = 0; hole < 2; ++hole) {
      float features[N * F];
      bool mask[N];
      for (float& f : features)
        f = normal(rng);
      for (int t = 0; t < N; ++t)
        mask[t] = t < multiplicity;
      if (hole != 0 && multiplicity > 2)
        mask[1] = false;  // a hole in the middle of the sequence
      float reference[2];
      runReference(model, features, mask, reference);
      for (uint32_t blockThreads : threadings) {
        float fast[2];
        runFast(model, features, mask, fast, blockThreads);
        double const diff = std::max(std::fabs(double(reference[0]) - fast[0]), std::fabs(double(reference[1]) - fast[1]));
        worst = std::max(worst, diff);
        if (multiplicity == N && hole == 0 && blockThreads == 256)
          std::cout << "full jet reference " << reference[0] << " " << reference[1] << " | kernel " << fast[0] << " "
                    << fast[1] << "\n";
      }
      ++jets;
    }
  }
  // Warp layout: up to kWarpJetsPerBlock independent jets, one per warp, each
  // with its own barrier.  Centres the parity argument of the whole design.
  for (uint32_t warps : {1u, 3u, uint32_t(deepsets::kWarpJetsPerBlock)}) {
    std::vector<std::vector<float>> warpFeatures(warps);
    std::vector<float const*> featurePtrs(warps);
    std::vector<bool const*> maskPtrs(warps);
    std::vector<std::array<bool, N>> maskStorage(warps);
    std::vector<float> warpProbs(2 * warps);
    for (uint32_t w = 0; w < warps; ++w) {
      warpFeatures[w].resize(N * F);
      for (float& f : warpFeatures[w])
        f = normal(rng);
      int const multiplicity = 1 + (int(w) * 15) / int(warps);
      for (int t = 0; t < N; ++t)
        maskStorage[w][t] = t < multiplicity;
      featurePtrs[w] = warpFeatures[w].data();
      maskPtrs[w] = maskStorage[w].data();
    }
    runWarp(model, featurePtrs.data(), maskPtrs.data(), warpProbs.data(), warps);
    for (uint32_t w = 0; w < warps; ++w) {
      float reference[2];
      runReference(model, featurePtrs[w], maskPtrs[w], reference);
      double const diff = std::max(std::abs(double(reference[0]) - warpProbs[2 * w]),
                                   std::abs(double(reference[1]) - warpProbs[2 * w + 1]));
      worst = std::max(worst, diff);
    }
  }
  std::cout << jets << " block jets + warp batches; worst |dp| = " << worst << "\n";
  if (worst > 1e-6) {
    std::cerr << "FAILED: disagreement exceeds 1e-6\n";
    return 1;
  }
  std::cout << "OK\n";
  return 0;
}
