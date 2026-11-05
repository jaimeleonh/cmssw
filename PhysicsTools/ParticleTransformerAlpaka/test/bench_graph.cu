// Time the graph interpreter on a GPU, one block per jet.
//
//   nvcc -O3 -std=c++20 -arch=sm_89 -Itest/cudastubs -I<src-root> \
//        -o bench_graph test/bench_graph.cu
//   ./bench_graph data/b_kinadd.nnm [JETS] [REPEATS] [--dump FILE]
//
// This exists because the interpreter's GPU behaviour cannot be reasoned about
// from the CPU stub, and a full CMSSW build is a slow way to ask a question
// like "did that change help".  It compiles the same headers the producer
// does, against a CUDA shim for the handful of Alpaka calls they use, so it is
// also what `ncu` can be pointed at.  `test/bench_compiled.cu` is its twin for
// the hand-written kernel, generating the same jets from the same seed, which
// is what makes the two directly comparable.
//
// Multiplicities are drawn uniformly from 1 to the model's padded length.
// These are random features, not physical ones: the point here is instruction
// mix and memory behaviour, not numerics.
//
// `--minmult N` raises the floor of that distribution, and it matters more
// than it looks.  The interpreter sizes every op by the *live* token count
// (resolveRows -> shared.tokens, and kRowsPairs is quadratic in it), while the
// compiled kernel in bench_compiled.cu always runs kMaxParticles rows and lets
// the mask zero the result.  So the default [1, 16] hands the interpreter a
// discount its rival does not take, and the two kernels look ~8% apart when at
// equal work they are ~70% apart.  Sweep --minmult before quoting any ratio,
// and say which occupancy the number came from.  bench_compiled.cu takes the
// same flag so the jets stay identical.
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/graph/ModelLoader.h"

#define CHECK(call)                                                                 \
  do {                                                                              \
    cudaError_t status = (call);                                                    \
    if (status != cudaSuccess) {                                                    \
      std::fprintf(stderr, "%s:%d %s\n", __FILE__, __LINE__, cudaGetErrorString(status)); \
      std::exit(1);                                                                 \
    }                                                                               \
  } while (0)

struct Acc {};

__global__ __launch_bounds__(nn::kThreadsPerJet) void graphKernel(nn::ModelView model,
                                                                 float const* features,
                                                                 float const* vectors,
                                                                 float const* masks,
                                                                 float* outputs,
                                                                 int outputStride,
                                                                 int jets) {
  __shared__ nn::GraphShared shared;
  Acc acc;
  int const lane = static_cast<int>(threadIdx.x);
  int const lanes = static_cast<int>(blockDim.x);
  int const tokens = model.maxTokens;
  int const cols = model.features.cols;
  int const components = model.vectors.valid() ? model.vectors.cols : 0;

  for (int jet = blockIdx.x; jet < jets; jet += gridDim.x) {
    float const* jetFeatures = features + static_cast<long>(jet) * tokens * cols;
    float const* jetVectors = components ? vectors + static_cast<long>(jet) * tokens * 4 : nullptr;
    float const* jetMask = masks + static_cast<long>(jet) * tokens;

    for (int index = lane; index < tokens; index += lanes)
      shared.mask[index] = jetMask[index];
    __syncthreads();

    nn::runGraph(acc, model, shared, nn::GraphInputs{jetFeatures, jetVectors});

    for (int index = lane; index < model.outputSize && index < nn::kMaxOutputs; index += lanes)
      outputs[static_cast<long>(jet) * outputStride + index] = shared.output[index];
    __syncthreads();
  }
}

// The producer's warp layout: one jet per warp, eight per block, with the
// inputs staged straight into each warp's arena and runGraph handed null
// pointers -- exactly GraphWarpJetKernel, so what this times and checks is the
// path production takes rather than the LoadInput path graphKernel uses.
__global__ __launch_bounds__(nn::kThreadsPerJet) void graphWarpKernel(nn::ModelView model,
                                                                     float const* features,
                                                                     float const* vectors,
                                                                     float const* masks,
                                                                     float* outputs,
                                                                     int outputStride,
                                                                     int jets) {
  __shared__ nn::GraphWarpBlock block;
  Acc acc;
  int const warp = static_cast<int>(threadIdx.x) / nn::kWarpLanes;
  int const lane = static_cast<int>(threadIdx.x) % nn::kWarpLanes;
  nn::SubGroup<Acc> const group{acc, lane, nn::kWarpLanes};
  nn::GraphWarpState& shared = block.jet[warp];
  int const tokens = model.maxTokens;
  int const cols = model.features.cols;
  int const components = model.vectors.valid() ? model.vectors.cols : 0;
  int const stride = static_cast<int>(gridDim.x) * nn::kWarpJetsPerBlock;

  for (int jet = static_cast<int>(blockIdx.x) * nn::kWarpJetsPerBlock + warp; jet < jets;
       jet += stride) {
    float const* jetFeatures = features + static_cast<long>(jet) * tokens * cols;
    float const* jetVectors = vectors + static_cast<long>(jet) * tokens * 4;
    float const* jetMask = masks + static_cast<long>(jet) * tokens;

    for (int index = lane; index < tokens * cols; index += nn::kWarpLanes) {
      int const particle = index / cols;
      int const c = index - particle * cols;
      shared.arena[model.features.offset + particle * model.features.stride + c] =
          jetFeatures[index];
    }
    for (int index = lane; index < tokens * components; index += nn::kWarpLanes) {
      int const particle = index / components;
      int const c = index - particle * components;
      shared.arena[model.vectors.offset + particle * model.vectors.stride + c] =
          jetVectors[particle * 4 + c];
    }
    for (int index = lane; index < tokens; index += nn::kWarpLanes)
      shared.mask[index] = jetMask[index];
    nn::barrier(group);

    nn::runGraph(group, model, shared, nn::GraphInputs{nullptr, nullptr});

    for (int index = lane; index < model.outputSize && index < nn::kMaxOutputs;
         index += nn::kWarpLanes)
      outputs[static_cast<long>(jet) * outputStride + index] = shared.output[index];
    nn::barrier(group);
  }
}

namespace {
  std::vector<char> readFile(char const* path) {
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input) {
      std::fprintf(stderr, "cannot open %s\n", path);
      std::exit(1);
    }
    auto const bytes = static_cast<std::size_t>(input.tellg());
    input.seekg(0);
    std::vector<char> data(bytes);
    input.read(data.data(), static_cast<std::streamsize>(bytes));
    return data;
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
  if (argc < 2) {
    std::fprintf(stderr,
                 "usage: bench_graph MODEL.nnm [JETS] [REPEATS] [--minmult N] [--wrap] [--warp] [--dump FILE]\n");
    return 2;
  }
  int const jets = argc > 2 ? std::atoi(argv[2]) : 50000;
  int const repeats = argc > 3 ? std::atoi(argv[3]) : 20;
  char const* dump = nullptr;
  int minmult = 1;
  bool wrap = false;
  bool warpLayout = false;
  for (int i = 1; i < argc; ++i)
    if (std::strcmp(argv[i], "--wrap") == 0)
      wrap = true;
  for (int i = 1; i < argc; ++i)
    if (std::strcmp(argv[i], "--warp") == 0)
      warpLayout = true;
  for (int i = 1; i < argc - 1; ++i) {
    if (std::strcmp(argv[i], "--dump") == 0)
      dump = argv[i + 1];
    else if (std::strcmp(argv[i], "--minmult") == 0)
      minmult = std::atoi(argv[i + 1]);
  }

  auto const image = readFile(argv[1]);
  nn::LoadedModel loaded;
  std::string error;
  if (!nn::loadModel(image.data(), image.size(), loaded, error)) {
    std::fprintf(stderr, "the model file was rejected: %s\n", error.c_str());
    return 1;
  }
  auto const& header = *loaded.header;
  int const tokens = static_cast<int>(header.maxTokens);
  int const cols = loaded.features.valid() ? loaded.features.cols : 0;
  int const components = loaded.vectors.valid() ? loaded.vectors.cols : 0;

  // ---------------------------------------------------------------- inputs
  std::mt19937 generator(20260921);
  minmult = std::max(1, std::min(minmult, tokens));
  std::uniform_int_distribution<int> multiplicity(minmult, tokens);
  std::printf("occupancy: multiplicity uniform in [%d, %d]\n", minmult, tokens);
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
        features[(static_cast<std::size_t>(jet) * tokens + particle) * cols + c] =
            feature(generator);
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

  // ---------------------------------------------------------------- upload
  auto* deviceOps = upload(loaded.ops, header.opCount);
  auto* deviceBuffers = upload(loaded.buffers, header.bufferCount);
  auto* deviceWeights = upload(loaded.weights, header.weightsCount ? header.weightsCount : 1);
  auto* deviceParams = upload(loaded.params, header.paramsCount ? header.paramsCount : 1);
  auto* deviceFeatures = upload(features.data(), features.size());
  auto* deviceVectors = upload(vectors.data(), vectors.size());
  auto* deviceMasks = upload(masks.data(), masks.size());

  int const outputStride = static_cast<int>(header.outputSize);
  float* deviceOutputs = nullptr;
  CHECK(cudaMalloc(&deviceOutputs, static_cast<std::size_t>(jets) * outputStride * sizeof(float)));

  auto const model = loaded.view(deviceWeights, deviceParams, deviceOps, deviceBuffers);

  if (warpLayout && !loaded.fitsWarpMode()) {
    std::fprintf(stderr, "--warp: this model's activations reach %d floats / %d int8 bytes, "
                         "past the per-warp arena of %d / %d\n",
                 loaded.floatsUsed, loaded.int8Used, nn::kWarpArenaFloats, nn::kWarpArenaInt8);
    return 1;
  }
  std::printf("model %s: %u ops, %u buffers, %d tokens, %d outputs, arena used %d floats, "
              "layout %s, shared %zu B/block\n",
              argv[1], header.opCount, header.bufferCount, tokens, outputStride,
              loaded.floatsUsed, warpLayout ? "warp" : "block",
              warpLayout ? sizeof(nn::GraphWarpBlock) : sizeof(nn::GraphShared));

  // ---------------------------------------------------------------- run
  int const blocks = warpLayout ? (jets + nn::kWarpJetsPerBlock - 1) / nn::kWarpJetsPerBlock : jets;
  auto launch = [&] {
    if (warpLayout)
      graphWarpKernel<<<blocks, nn::kThreadsPerJet>>>(model, deviceFeatures, deviceVectors,
                                                      deviceMasks, deviceOutputs, outputStride,
                                                      jets);
    else
      graphKernel<<<blocks, nn::kThreadsPerJet>>>(model, deviceFeatures,
                                                  components ? deviceVectors : nullptr,
                                                  deviceMasks, deviceOutputs, outputStride, jets);
  };
  launch();
  CHECK(cudaDeviceSynchronize());
  CHECK(cudaGetLastError());

  cudaEvent_t start, stop;
  CHECK(cudaEventCreate(&start));
  CHECK(cudaEventCreate(&stop));
  float best = 1e30f, total = 0.f;
  for (int repeat = 0; repeat < repeats; ++repeat) {
    CHECK(cudaEventRecord(start));
    launch();
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
    std::vector<float> outputs(static_cast<std::size_t>(jets) * outputStride);
    CHECK(cudaMemcpy(outputs.data(), deviceOutputs, outputs.size() * sizeof(float),
                     cudaMemcpyDeviceToHost));
    std::ofstream out(dump);
    out.precision(9);
    for (int jet = 0; jet < jets; ++jet) {
      for (int index = 0; index < outputStride; ++index)
        out << (index ? " " : "") << outputs[static_cast<std::size_t>(jet) * outputStride + index];
      out << "\n";
    }
    // The inputs too, so the NumPy interpreter can be given exactly these jets.
    std::ofstream inputs(std::string(dump) + ".inputs");
    inputs.precision(9);
    for (int jet = 0; jet < jets; ++jet) {
      inputs << tokens;
      for (int t = 0; t < tokens; ++t)
        inputs << " " << masks[static_cast<std::size_t>(jet) * tokens + t];
      for (int i = 0; i < tokens * cols; ++i)
        inputs << " " << features[static_cast<std::size_t>(jet) * tokens * cols + i];
      if (components)
        for (int i = 0; i < tokens * 4; ++i)
          inputs << " " << vectors[static_cast<std::size_t>(jet) * tokens * 4 + i];
      inputs << "\n";
    }
  }
  return 0;
}
