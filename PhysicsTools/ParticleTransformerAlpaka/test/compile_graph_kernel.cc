// Standalone regression test for the graph interpreter.
//
// It runs the same jets twice: once through the hand-written
// ParticleTransformer kernel reading the generated header and its blobs, and
// once through the generic interpreter reading a .nnm file exported from that
// same quantized model.  The two share no code beyond the Alpaka stub, so
// agreement means the op encoding, the arena layout, the tiled regions and
// every layer implementation are all consistent with the kernel they replace.
//
// Both are driven with real host threads and a std::barrier, so the block
// cooperation is exercised rather than simulated away, at several block sizes.
#include <algorithm>
#include <barrier>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/ParticleTransformerFastKernel.h"
#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/ModelTranspose.h"
#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/graph/ModelLoader.h"

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

namespace {

  std::vector<char> readFile(char const* path) {
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input)
      throw std::runtime_error(std::string("cannot open ") + path);
    auto const bytes = static_cast<std::size_t>(input.tellg());
    input.seekg(0);
    std::vector<char> data(bytes);
    input.read(data.data(), static_cast<std::streamsize>(bytes));
    if (!input)
      throw std::runtime_error(std::string("cannot read ") + path);
    return data;
  }

  template <typename T>
  std::vector<T> readExact(char const* path, std::size_t count) {
    auto const raw = readFile(path);
    if (raw.size() != count * sizeof(T))
      throw std::runtime_error(std::string(path) + " has the wrong size");
    std::vector<T> values(count);
    std::memcpy(values.data(), raw.data(), raw.size());
    return values;
  }

  // The reference: the package's existing block kernel.
  void runHandWritten(part::generated::ModelView model,
                      float const* features,
                      float const* vectors,
                      bool const* mask,
                      float* probabilities,
                      uint32_t blockThreads) {
    auto shared = std::make_unique<part::JetShared>();
    std::memset(shared.get(), 0, sizeof(part::JetShared));
    std::copy_n(features, part::kMaxParticles * part::kInputFeatures, shared->na);
    std::copy_n(vectors, part::kMaxParticles * part::kFourVector, shared->p4);
    for (int particle = 0; particle < part::kMaxParticles; ++particle)
      shared->mask[particle] = mask[particle] ? 1.f : 0.f;

    std::barrier blockBarrier(static_cast<std::ptrdiff_t>(blockThreads));
    {
      std::vector<std::jthread> threads;
      threads.reserve(blockThreads);
      for (uint32_t lane = 0; lane < blockThreads; ++lane)
        threads.emplace_back([&, lane] {
          part::runParticleTransformer(Acc{lane, blockThreads, &blockBarrier, nullptr}, model,
                                       *shared);
        });
    }
    probabilities[0] = shared->logits[0];
    probabilities[1] = shared->logits[1];
  }

  // The interpreter, over the exported file.
  void runInterpreter(nn::ModelView const& model,
                      float const* features,
                      float const* vectors,
                      bool const* mask,
                      float* probabilities,
                      uint32_t blockThreads) {
    auto shared = std::make_unique<nn::GraphShared>();
    std::memset(shared.get(), 0, sizeof(nn::GraphShared));
    for (int particle = 0; particle < model.maxTokens; ++particle)
      shared->mask[particle] = mask[particle] ? 1.f : 0.f;

    nn::GraphInputs const inputs{features, vectors};
    std::barrier blockBarrier(static_cast<std::ptrdiff_t>(blockThreads));
    {
      std::vector<std::jthread> threads;
      threads.reserve(blockThreads);
      for (uint32_t lane = 0; lane < blockThreads; ++lane)
        threads.emplace_back([&, lane] {
          nn::runGraph(Acc{lane, blockThreads, &blockBarrier, nullptr}, model, *shared, inputs);
        });
    }
    probabilities[0] = shared->output[0];
    probabilities[1] = shared->output[1];
  }

}  // namespace

int main(int argc, char** argv) {
  if (argc != 5) {
    std::cerr << "usage: compile_graph_kernel MODEL.nnm WEIGHTS PARAMS PRINTED_JET\n";
    return 2;
  }

  auto const image = readFile(argv[1]);
  nn::LoadedModel loaded;
  std::string error;
  if (!nn::loadModel(image.data(), image.size(), loaded, error)) {
    std::cerr << "the model file was rejected: " << error << "\n";
    return 1;
  }
  auto const model = loaded.view(loaded.weights, loaded.params, loaded.ops, loaded.buffers);

  std::cout << "graph shared bytes per jet: " << sizeof(nn::GraphShared) << "\n"
            << "ops " << loaded.header->opCount << ", buffers " << loaded.header->bufferCount
            << ", arena " << loaded.header->arenaFloats << " floats + "
            << loaded.header->arenaInt8 << " int8 bytes\n";

  // The hand-written kernel, over the same quantized model in its own format.
  auto weights = readExact<std::int8_t>(argv[2], part::generated::weights_size);
  auto params = readExact<float>(argv[3], part::generated::params_size);
  std::vector<std::int8_t> transposed(part::generated::weights_size, 0);
  part::transposeWeights(weights.data(), params.data(), transposed.data());
  part::generated::ModelView const handWritten{transposed.data(), params.data()};

  // The printed jet.
  float features[part::kMaxParticles * part::kInputFeatures] = {};
  float vectors[part::kMaxParticles * part::kFourVector] = {};
  bool mask[part::kMaxParticles] = {};
  {
    std::ifstream jet(argv[4]);
    std::string tag;
    int particle = 0;
    while (jet >> tag && particle < part::kMaxParticles) {
      double row[21];
      for (double& value : row)
        jet >> value;
      for (int f = 0; f < part::kInputFeatures; ++f)
        features[particle * part::kInputFeatures + f] = static_cast<float>(row[f + 2]);
      for (int c = 0; c < part::kFourVector; ++c)
        vectors[particle * part::kFourVector + c] = static_cast<float>(row[c + 17]);
      mask[particle] = true;
      ++particle;
    }
    if (particle == 0) {
      std::cerr << "no particles were read from " << argv[4] << "\n";
      return 1;
    }
  }

  float reference[2];
  float interpreted[2];
  runHandWritten(handWritten, features, vectors, mask, reference, part::kThreadsPerJet);
  runInterpreter(model, features, vectors, mask, interpreted, part::kThreadsPerJet);
  std::cout << std::setprecision(9) << "printed jet, hand-written kernel " << reference[0] << " "
            << reference[1] << "\nprinted jet, graph interpreter     " << interpreted[0] << " "
            << interpreted[1] << "\n";

  // Randomised jets at several block sizes: every multiplicity, half of them
  // with a hole punched in the mask, so the masked-key and partial-tile paths
  // are exercised rather than assumed.
  std::mt19937 generator(20240917u);
  std::uniform_real_distribution<float> feature(-2.f, 2.f);
  std::uniform_real_distribution<float> momentum(-30.f, 30.f);
  float worst = std::abs(interpreted[0] - reference[0]);
  int comparisons = 1;

  for (int multiplicity = 1; multiplicity <= part::kMaxParticles; ++multiplicity) {
    for (int variant = 0; variant < 2; ++variant) {
      float randomFeatures[part::kMaxParticles * part::kInputFeatures] = {};
      float randomVectors[part::kMaxParticles * part::kFourVector] = {};
      bool randomMask[part::kMaxParticles] = {};
      for (int particle = 0; particle < multiplicity; ++particle) {
        randomMask[particle] = true;
        for (int f = 0; f < part::kInputFeatures; ++f)
          randomFeatures[particle * part::kInputFeatures + f] = feature(generator);
        float const px = momentum(generator);
        float const py = momentum(generator);
        float const pz = momentum(generator);
        randomVectors[particle * part::kFourVector + 0] = px;
        randomVectors[particle * part::kFourVector + 1] = py;
        randomVectors[particle * part::kFourVector + 2] = pz;
        // Massless, which is the convention the input transform uses.
        randomVectors[particle * part::kFourVector + 3] =
            std::sqrt(px * px + py * py + pz * pz);
      }
      // A hole in the middle of the mask: `active` is the last unmasked index
      // plus one, so the padded tail and an interior gap are different cases.
      if (variant == 1 && multiplicity > 2)
        randomMask[multiplicity / 2] = false;

      float expected[2];
      runHandWritten(handWritten, randomFeatures, randomVectors, randomMask, expected,
                     part::kThreadsPerJet);
      for (uint32_t threads : {1u, 32u, 96u, 256u}) {
        float actual[2];
        runInterpreter(model, randomFeatures, randomVectors, randomMask, actual, threads);
        worst = std::max(worst, std::abs(actual[0] - expected[0]));
        worst = std::max(worst, std::abs(actual[1] - expected[1]));
        ++comparisons;
      }
    }
  }

  std::cout << "compared " << comparisons << " evaluations against the hand-written kernel\n"
            << "largest difference: " << worst << "\n";

  // Quantized inference is not bitwise stable across reduction orders: a few
  // 1e-7 before rounding occasionally moves an activation across an INT8
  // boundary, and each flip shifts one input by a whole quantization step.
  // This is the same tolerance compile_kernel_int8.cc uses against the scalar
  // reference, and for the same reason.
  constexpr float kTolerance = 3.e-3f;
  if (!(worst < kTolerance)) {
    std::cerr << "FAILED: the interpreter and the hand-written kernel differ by " << worst
              << ", above the " << kTolerance << " quantization tolerance\n";
    return 1;
  }
  std::cout << "OK\n";
  return 0;
}
