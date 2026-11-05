// The same benchmark for the hand-written compiled kernel, so the two can be
// compared on identical jets.  Input generation is copied verbatim from
// bench_graph.cu; keep them in step.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/ModelTranspose.h"
#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/ParticleTransformerFastKernel.h"

#define CHECK(call)                                                                 \
  do {                                                                              \
    cudaError_t status = (call);                                                    \
    if (status != cudaSuccess) {                                                    \
      std::fprintf(stderr, "%s:%d %s\n", __FILE__, __LINE__, cudaGetErrorString(status)); \
      std::exit(1);                                                                 \
    }                                                                               \
  } while (0)

struct Acc {};

__global__ __launch_bounds__(256) void compiledKernel(part::generated::ModelView model,
                                                      float const* features,
                                                      float const* vectors,
                                                      float const* masks,
                                                      float* outputs,
                                                      int jets) {
  __shared__ part::JetShared shared;
  Acc acc;
  int const lane = static_cast<int>(threadIdx.x);
  int const lanes = static_cast<int>(blockDim.x);

  for (int jet = blockIdx.x; jet < jets; jet += gridDim.x) {
    float const* jetFeatures = features + static_cast<long>(jet) * part::kMaxParticles * part::kInputFeatures;
    float const* jetVectors = vectors + static_cast<long>(jet) * part::kMaxParticles * part::kFourVector;
    float const* jetMask = masks + static_cast<long>(jet) * part::kMaxParticles;

    for (int index = lane; index < part::kMaxParticles * part::kInputFeatures; index += lanes)
      shared.na[index] = jetFeatures[index];
    for (int index = lane; index < part::kMaxParticles * part::kFourVector; index += lanes)
      shared.p4[index] = jetVectors[index];
    for (int index = lane; index < part::kMaxParticles; index += lanes)
      shared.mask[index] = jetMask[index];
    __syncthreads();

    part::runParticleTransformer(acc, model, shared);

    if (lane == 0) {
      outputs[static_cast<long>(jet) * 2 + 0] = shared.logits[0];
      outputs[static_cast<long>(jet) * 2 + 1] = shared.logits[1];
    }
    __syncthreads();
  }
}

namespace {
  template <typename T>
  std::vector<T> readBinary(char const* path, std::size_t count) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
      std::fprintf(stderr, "cannot open %s\n", path);
      std::exit(1);
    }
    std::vector<T> values(count);
    input.read(reinterpret_cast<char*>(values.data()), static_cast<std::streamsize>(count * sizeof(T)));
    if (!input) {
      std::fprintf(stderr, "%s is shorter than the %zu elements the model declares\n", path, count);
      std::exit(1);
    }
    return values;
  }

  template <typename T>
  T* upload(T const* host, std::size_t count) {
    T* device = nullptr;
    CHECK(cudaMalloc(&device, count * sizeof(T)));
    CHECK(cudaMemcpy(device, host, count * sizeof(T), cudaMemcpyHostToDevice));
    return device;
  }
}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr,
                 "usage: bench_compiled WEIGHTS PARAMS [JETS] [REPEATS] [--minmult N] [--wrap] [--dump FILE]\n");
    return 2;
  }
  int const jets = argc > 3 ? std::atoi(argv[3]) : 50000;
  int const repeats = argc > 4 ? std::atoi(argv[4]) : 20;
  char const* dump = nullptr;
  // Accepted so the two benchmarks can be swept together and keep generating
  // the same jets.  It does not change this kernel's cost -- it always runs
  // kMaxParticles rows -- which is exactly the asymmetry the sweep exposes.
  int minmult = 1;
  bool wrap = false;
  for (int i = 1; i < argc; ++i)
    if (std::strcmp(argv[i], "--wrap") == 0)
      wrap = true;
  for (int i = 1; i < argc - 1; ++i) {
    if (std::strcmp(argv[i], "--dump") == 0)
      dump = argv[i + 1];
    else if (std::strcmp(argv[i], "--minmult") == 0)
      minmult = std::atoi(argv[i + 1]);
  }

  auto weights = readBinary<std::int8_t>(argv[1], part::generated::weights_size);
  auto params = readBinary<float>(argv[2], part::generated::params_size);
  std::vector<std::int8_t> transposed(part::generated::weights_size, 0);
  auto const covered = part::transposeWeights(weights.data(), params.data(), transposed.data());
  if (covered != part::generated::weights_size) {
    std::fprintf(stderr, "the tensor walk covered %zu of %zu weights\n", covered,
                 part::generated::weights_size);
    return 1;
  }

  int const tokens = part::kMaxParticles;
  int const cols = part::kInputFeatures;

  // Identical to bench_graph.cu, same seed, so the two see the same jets.
  std::mt19937 generator(20260921);
  minmult = std::max(1, std::min(minmult, tokens));
  std::uniform_int_distribution<int> multiplicity(minmult, tokens);
  std::normal_distribution<float> feature(0.f, 1.f);
  std::uniform_real_distribution<float> momentum(-30.f, 30.f);

  std::vector<float> features(static_cast<std::size_t>(jets) * tokens * cols, 0.f);
  std::vector<float> vectors(static_cast<std::size_t>(jets) * tokens * 4, 0.f);
  std::vector<float> masks(static_cast<std::size_t>(jets) * tokens, 0.f);
  for (int jet = 0; jet < jets; ++jet) {
    int const live = multiplicity(generator);
    for (int particle = 0; particle < live; ++particle) {
      masks[static_cast<std::size_t>(jet) * tokens + particle] = 1.f;
      for (int c = 0; c < cols; ++c)
        features[(static_cast<std::size_t>(jet) * tokens + particle) * cols + c] = feature(generator);
      float px = momentum(generator), py = momentum(generator), pz = momentum(generator);
      auto* p4 = &vectors[(static_cast<std::size_t>(jet) * tokens + particle) * 4];
      p4[0] = px; p4[1] = py; p4[2] = pz;
      p4[3] = std::sqrt(px * px + py * py + pz * pz);
    }
    // The model's YAML declares `pad_mode: wrap` for features, points and
    // vectors (only the mask is padded with a constant), so in production the
    // padded slots hold *copies of real particles*, not zeros.  --wrap
    // reproduces that; the zero padding this generator used by default is not
    // what either kernel actually sees.
    if (wrap && live > 0) {
      for (int particle = live; particle < tokens; ++particle) {
        int const source = particle % live;
        for (int c = 0; c < cols; ++c)
          features[(static_cast<std::size_t>(jet) * tokens + particle) * cols + c] =
              features[(static_cast<std::size_t>(jet) * tokens + source) * cols + c];
        for (int c = 0; c < 4; ++c)
          vectors[(static_cast<std::size_t>(jet) * tokens + particle) * 4 + c] =
              vectors[(static_cast<std::size_t>(jet) * tokens + source) * 4 + c];
      }
    }
  }

  auto* deviceWeights = upload(transposed.data(), transposed.size());
  auto* deviceParams = upload(params.data(), params.size());
  auto* deviceFeatures = upload(features.data(), features.size());
  auto* deviceVectors = upload(vectors.data(), vectors.size());
  auto* deviceMasks = upload(masks.data(), masks.size());
  float* deviceOutputs = nullptr;
  CHECK(cudaMalloc(&deviceOutputs, static_cast<std::size_t>(jets) * 2 * sizeof(float)));

  part::generated::ModelView model{deviceWeights, deviceParams};
  std::printf("compiled kernel: shared %zu B\n", sizeof(part::JetShared));

  compiledKernel<<<jets, 256>>>(model, deviceFeatures, deviceVectors, deviceMasks,
                                deviceOutputs, jets);
  CHECK(cudaDeviceSynchronize());
  CHECK(cudaGetLastError());

  cudaEvent_t start, stop;
  CHECK(cudaEventCreate(&start));
  CHECK(cudaEventCreate(&stop));
  float best = 1e30f, total = 0.f;
  for (int repeat = 0; repeat < repeats; ++repeat) {
    CHECK(cudaEventRecord(start));
    compiledKernel<<<jets, 256>>>(model, deviceFeatures, deviceVectors, deviceMasks,
                                  deviceOutputs, jets);
    CHECK(cudaEventRecord(stop));
    CHECK(cudaEventSynchronize(stop));
    float milliseconds = 0.f;
    CHECK(cudaEventElapsedTime(&milliseconds, start, stop));
    best = milliseconds < best ? milliseconds : best;
    total += milliseconds;
  }
  std::printf("%d jets: best %.3f ms, mean %.3f ms over %d launches\n",
              jets, best, total / repeats, repeats);

  if (dump != nullptr) {
    std::vector<float> outputs(static_cast<std::size_t>(jets) * 2);
    CHECK(cudaMemcpy(outputs.data(), deviceOutputs, outputs.size() * sizeof(float),
                     cudaMemcpyDeviceToHost));
    std::ofstream out(dump);
    out.precision(9);
    for (int jet = 0; jet < jets; ++jet)
      out << outputs[static_cast<std::size_t>(jet) * 2] << " "
          << outputs[static_cast<std::size_t>(jet) * 2 + 1] << "\n";
  }
  return 0;
}
