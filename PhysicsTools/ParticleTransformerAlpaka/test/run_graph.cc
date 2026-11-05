// Run any .nnm through the C++ interpreter and print the outputs.
//
// This is the tool that lets a converted model be checked against the NumPy
// interpreter and against the framework it came from, without a GPU or a CMSSW
// build.  It is also the quickest way to see what a model file actually does
// when a conversion looks wrong.
//
//   run_graph MODEL.nnm INPUTS.txt [THREADS]
//
// INPUTS.txt holds one jet per line:
//
//   <particles> <mask_0..mask_{n-1}> <features row-major> [<four-vectors row-major>]
//
// with `particles` the padded sequence length the model declares.  The
// four-vector block is optional and only read when the model has an input for
// it.  Each line produces one line of outputs.
#include <barrier>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

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
    return data;
  }

}  // namespace

int main(int argc, char** argv) {
  bool prestaged = false;
  std::vector<char*> args;
  for (int i = 0; i < argc; ++i) {
    if (std::strcmp(argv[i], "--prestaged") == 0)
      prestaged = true;
    else
      args.push_back(argv[i]);
  }
  argc = static_cast<int>(args.size());
  argv = args.data();

  if (argc < 3 || argc > 4) {
    std::cerr << "usage: run_graph MODEL.nnm INPUTS.txt [THREADS] [--prestaged]\n";
    return 2;
  }
  auto const threads = argc == 4 ? static_cast<uint32_t>(std::stoul(argv[3])) : 256u;

  auto const image = readFile(argv[1]);
  nn::LoadedModel loaded;
  std::string error;
  if (!nn::loadModel(image.data(), image.size(), loaded, error)) {
    std::cerr << "the model file was rejected: " << error << "\n";
    return 1;
  }
  auto const model = loaded.view(loaded.weights, loaded.params, loaded.ops, loaded.buffers);
  auto const& header = *loaded.header;

  int const tokens = static_cast<int>(header.maxTokens);
  int const features = loaded.features.valid() ? loaded.features.cols : 0;
  int const components = loaded.vectors.valid() ? loaded.vectors.cols : 0;

  std::ifstream inputs(argv[2]);
  if (!inputs) {
    std::cerr << "cannot open " << argv[2] << "\n";
    return 1;
  }

  std::string line;
  int jet = 0;
  while (std::getline(inputs, line)) {
    if (line.empty())
      continue;
    std::istringstream stream(line);
    int declared = 0;
    stream >> declared;
    if (declared != tokens) {
      std::cerr << "line " << jet << " declares " << declared << " particles but the model was "
                << "packed for " << tokens << "\n";
      return 1;
    }

    std::vector<float> mask(tokens, 0.f);
    std::vector<float> featureValues(static_cast<std::size_t>(tokens) * features, 0.f);
    std::vector<float> vectorValues(static_cast<std::size_t>(tokens) * (components ? components : 1),
                                    0.f);
    for (auto& value : mask)
      stream >> value;
    for (auto& value : featureValues)
      stream >> value;
    if (components != 0)
      for (auto& value : vectorValues)
        stream >> value;
    if (!stream) {
      std::cerr << "line " << jet << " has too few numbers\n";
      return 1;
    }

    auto shared = std::make_unique<nn::GraphShared>();
    std::memset(shared.get(), 0, sizeof(nn::GraphShared));
    for (int particle = 0; particle < tokens && particle < nn::kMaxTokens; ++particle)
      shared->mask[particle] = mask[particle];

    // --prestaged reproduces what GraphInferenceProducer does: write the inputs
    // straight into their arena buffers and hand runGraph null pointers, rather
    // than letting LoadInput fetch them.  This is not a variation for its own
    // sake -- it is the only path the producer takes, and until it was covered
    // here every test and benchmark passed real pointers, which is how a
    // runGraph that zeroed the arena *after* the caller had staged its inputs
    // survived: it silently reduced the network to a function of the token
    // count and no oracle in this package could see it.
    nn::GraphInputs jetInputs{featureValues.data(),
                              components != 0 ? vectorValues.data() : nullptr};
    if (prestaged) {
      if (model.features.valid())
        for (int row = 0; row < model.features.rows; ++row)
          for (int col = 0; col < model.features.cols; ++col)
            shared->arena[model.features.offset + row * model.features.stride + col] =
                featureValues[static_cast<std::size_t>(row) * features + col];
      if (components != 0 && model.vectors.valid())
        for (int row = 0; row < model.vectors.rows; ++row)
          for (int col = 0; col < model.vectors.cols; ++col)
            shared->arena[model.vectors.offset + row * model.vectors.stride + col] =
                vectorValues[static_cast<std::size_t>(row) * components + col];
      jetInputs = nn::GraphInputs{nullptr, nullptr};
    }
    std::barrier blockBarrier(static_cast<std::ptrdiff_t>(threads));
    {
      std::vector<std::jthread> workers;
      workers.reserve(threads);
      for (uint32_t lane = 0; lane < threads; ++lane)
        workers.emplace_back([&, lane] {
          nn::runGraph(Acc{lane, threads, &blockBarrier, nullptr}, model, *shared, jetInputs);
        });
    }

    std::cout << std::setprecision(9);
    for (uint32_t index = 0; index < header.outputSize; ++index)
      std::cout << (index ? " " : "") << shared->output[index];
    std::cout << "\n";
    ++jet;
  }
  return 0;
}
