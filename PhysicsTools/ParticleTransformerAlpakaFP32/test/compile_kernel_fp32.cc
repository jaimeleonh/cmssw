// Standalone regression test.
//
// It runs the same jets through the scalar reference graph and through the
// block-cooperative production kernel (simulated with real host threads and a
// std::barrier), and reports the largest disagreement.  The two paths share no
// code beyond the generated tensor descriptors: the reference reads the
// original row-major weights, the fast kernel reads the transposed image.
#include <algorithm>
#include <barrier>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <thread>
#include <vector>

#include "PhysicsTools/ParticleTransformerAlpakaFP32/interface/alpaka/ParticleTransformerKernel.h"
#include "PhysicsTools/ParticleTransformerAlpakaFP32/interface/alpaka/ParticleTransformerFastKernel.h"
#include "PhysicsTools/ParticleTransformerAlpakaFP32/interface/alpaka/ModelTranspose.h"

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

std::vector<float> readFloats(char const* path, std::size_t expected) {
  std::ifstream input(path, std::ios::binary | std::ios::ate);
  if (!input || static_cast<std::size_t>(input.tellg()) != expected * sizeof(float))
    throw std::runtime_error("model file has the wrong size");
  input.seekg(0);
  std::vector<float> values(expected);
  input.read(reinterpret_cast<char*>(values.data()), static_cast<std::streamsize>(values.size() * sizeof(float)));
  if (!input)
    throw std::runtime_error("model file read failed");
  return values;
}

void runReference(partfp32::generated::ModelView model,
                  float const* features,
                  float const* vectors,
                  bool const* mask,
                  float* probabilities) {
  partfp32::Kernel{}(Acc{}, model, features, vectors, mask, probabilities);
}

void runFast(partfp32::generated::ModelView model,
             float const* features,
             float const* vectors,
             bool const* mask,
             float* probabilities,
             uint32_t blockThreads) {
  auto shared = std::make_unique<partfp32::JetShared>();
  std::memset(shared.get(), 0, sizeof(partfp32::JetShared));
  std::copy_n(features, partfp32::kMaxParticles * partfp32::kInputFeatures, shared->na);
  std::copy_n(vectors, partfp32::kMaxParticles * partfp32::kFourVector, shared->p4);
  for (int particle = 0; particle < partfp32::kMaxParticles; ++particle)
    shared->mask[particle] = mask[particle] ? 1.f : 0.f;

  std::barrier blockBarrier(static_cast<std::ptrdiff_t>(blockThreads));
  {
    std::vector<std::jthread> threads;
    threads.reserve(blockThreads);
    for (uint32_t lane = 0; lane < blockThreads; ++lane)
      threads.emplace_back([&, lane] {
        partfp32::runParticleTransformer(Acc{lane, blockThreads, &blockBarrier, nullptr}, model, *shared);
      });
  }
  probabilities[0] = shared->logits[0];
  probabilities[1] = shared->logits[1];
}

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "usage: compile_kernel_fp32 WEIGHTS PARAMS\n";
    return 2;
  }
  auto weights = readFloats(argv[1], partfp32::generated::weights_size);
  auto params = readFloats(argv[2], partfp32::generated::params_size);

  std::vector<float> transposed(partfp32::generated::weights_size, 0.f);
  auto const covered = partfp32::transposeWeights(weights.data(), transposed.data());
  if (covered != partfp32::generated::weights_size) {
    std::cerr << "the tensor walk covered " << covered << " of " << partfp32::generated::weights_size
              << " weights\n";
    return 1;
  }
  std::cout << "shared bytes per jet: " << sizeof(partfp32::JetShared) << "\n";
  if (sizeof(partfp32::JetShared) > 48u * 1024u) {
    std::cerr << "the shared footprint exceeds the 48 kB static limit\n";
    return 1;
  }

  partfp32::generated::ModelView reference{weights.data(), params.data()};
  partfp32::generated::ModelView fast{transposed.data(), params.data()};

  float features[16 * 15] = {};
  float vectors[16 * 4] = {};
  bool mask[16] = {};
  float referenceProbabilities[2];
  float fastProbabilities[2];
  constexpr float rows[14][21] = {
      {0.013672f,-0.013212f,1.202409f,1.105805f,2.427650f,2.430056f,-0.723950f,0.013672f,-0.013212f,0.104687f,0.062500f,1,1,0,0,0,0,-28.023132f,-12.039688f,18.866461f,35.863541f},
      {-0.008145f,0.026058f,0.728588f,0.624072f,1.953830f,1.948323f,-0.690794f,-0.008145f,0.026058f,-0.045313f,0.062500f,1,1,0,0,0,0,-13.990071f,-6.672922f,9.192499f,18.020878f},
      {0.052942f,-0.021938f,0.717206f,0.635448f,1.942447f,1.959699f,-0.570771f,0.052942f,-0.021938f,0.054688f,0.062500f,1,1,0,0,0,0,-14.063565f,-5.897343f,10.144864f,18.316133f},
      {-0.112864f,-0.026302f,0.549435f,0.410465f,1.774676f,1.734716f,-0.336446f,-0.112864f,-0.026302f,-0.045313f,0.093750f,1,1,0,0,0,0,-11.086554f,-4.592202f,5.692147f,13.281587f},
      {0.039852f,0.056602f,0.220432f,0.133641f,1.445674f,1.457892f,-0.523106f,0.039852f,0.056602f,0.054688f,0.109375f,1,1,0,0,0,0,-6.667628f,-3.434055f,4.871789f,8.943395f},
      {-0.095411f,-0.109205f,0.092807f,-0.040840f,1.318049f,1.283411f,-0.219945f,-0.095411f,-0.109205f,0.204687f,0.203125f,1,1,0,0,0,0,-5.952474f,-1.905402f,3.085850f,6.970292f},
      {0.096575f,0.165684f,-0.219594f,-0.283980f,1.005648f,1.040271f,-0.032896f,0.096575f,0.165684f,-0.145313f,0.125000f,1,1,0,0,0,0,-3.335543f,-2.207748f,2.873172f,4.924949f},
      {-0.021235f,0.122051f,-0.481879f,-0.591025f,0.743362f,0.733226f,-0.304461f,-0.021235f,0.122051f,-0.145313f,0.187500f,1,1,0,0,0,0,-2.357210f,-1.416355f,1.589214f,3.176177f},
      {-0.191404f,0.100235f,-0.548597f,-0.709327f,0.676645f,0.614924f,0.064246f,-0.191404f,0.100235f,-0.145313f,0.203125f,1,1,0,0,0,0,-2.170497f,-1.240541f,0.971980f,2.682302f},
      {-0.169588f,-0.017575f,-0.622349f,-0.777402f,0.602893f,0.546849f,-0.118017f,-0.169588f,-0.017575f,0.054688f,0.203125f,1,1,0,0,0,0,-2.071136f,-0.879145f,0.927662f,2.433733f},
      {0.083485f,0.095871f,-0.622349f,-0.692041f,0.602893f,0.632210f,-0.291495f,0.083485f,0.095871f,2.304688f,0,0.949219f,0,0,1,0,0,-1.958300f,-1.107953f,1.580034f,2.749365f},
      {-0.012508f,-0.187745f,-0.906174f,-1.012244f,0.319067f,0.312007f,-0.047357f,-0.012508f,-0.187745f,2.304688f,0,0.292969f,0,1,0,0,0,-1.460069f,-0.343801f,0.881996f,1.740091f},
      {0.275471f,-0.109205f,-1.033800f,-1.018251f,0.191442f,0.306001f,0.385311f,0.275471f,-0.109205f,2.304688f,0,0.183594f,0,1,0,0,0,-1.190495f,-0.381080f,1.189073f,1.725223f},
      {0.240565f,-0.078661f,-1.033800f,-1.034864f,0.191442f,0.289387f,0.212395f,0.240565f,-0.078661f,2.304688f,0,0.839844f,0,0,1,0,0,-1.178302f,-0.417259f,1.129563f,1.684759f},
  };
  for (int particle = 0; particle < 14; ++particle) {
    mask[particle] = true;
    for (int feature = 0; feature < 15; ++feature)
      features[particle * 15 + feature] = rows[particle][feature + 2];
    for (int component = 0; component < 4; ++component)
      vectors[particle * 4 + component] = rows[particle][component + 17];
  }

  runReference(reference, features, vectors, mask, referenceProbabilities);
  runFast(fast, features, vectors, mask, fastProbabilities, partfp32::kThreadsPerJet);
  std::cout << std::setprecision(9) << "printed jet reference " << referenceProbabilities[0] << " "
            << referenceProbabilities[1] << "\n"
            << "printed jet fast      " << fastProbabilities[0] << " " << fastProbabilities[1] << "\n";

  bool valid = std::abs(referenceProbabilities[0] - 0.25344047f) < 2.e-5f &&
               std::abs(fastProbabilities[0] - referenceProbabilities[0]) < 2.e-5f &&
               std::abs(fastProbabilities[1] - referenceProbabilities[1]) < 2.e-5f;
  float worst = std::max(std::abs(fastProbabilities[0] - referenceProbabilities[0]),
                         std::abs(fastProbabilities[1] - referenceProbabilities[1]));

  // Randomised jets, every multiplicity, both block shapes, and holes in the
  // mask so that the "skip the padded tail" logic is exercised properly.
  std::mt19937 generator(12345u);
  std::uniform_real_distribution<float> featureDistribution(-2.f, 2.f);
  std::uniform_real_distribution<float> momentumDistribution(-20.f, 20.f);
  for (int testCase = 0; testCase < 34; ++testCase) {
    int const validParticles = 1 + testCase % 16;
    bool const punchHole = testCase >= 16;
    std::fill_n(features, 16 * 15, 0.f);
    std::fill_n(vectors, 16 * 4, 0.f);
    std::fill_n(mask, 16, false);
    for (int particle = 0; particle < validParticles; ++particle) {
      for (int feature = 0; feature < 15; ++feature)
        features[particle * 15 + feature] = featureDistribution(generator);
      float const px = momentumDistribution(generator);
      float const py = momentumDistribution(generator);
      float const pz = momentumDistribution(generator);
      vectors[particle * 4] = px;
      vectors[particle * 4 + 1] = py;
      vectors[particle * 4 + 2] = pz;
      vectors[particle * 4 + 3] = std::sqrt(px * px + py * py + pz * pz) + 1.f;
      mask[particle] = true;
    }
    if (punchHole && validParticles > 2)
      mask[validParticles / 2] = false;

    runReference(reference, features, vectors, mask, referenceProbabilities);
    for (uint32_t blockThreads : {1u, 32u, 96u, 256u}) {
      runFast(fast, features, vectors, mask, fastProbabilities, blockThreads);
      float const difference = std::max(std::abs(fastProbabilities[0] - referenceProbabilities[0]),
                                        std::abs(fastProbabilities[1] - referenceProbabilities[1]));
      worst = std::max(worst, difference);
      if (!(difference < 2.e-5f)) {
        std::cerr << "case " << testCase << " with " << validParticles << " particles and "
                  << blockThreads << " threads: " << referenceProbabilities[0] << " vs "
                  << fastProbabilities[0] << "\n";
        valid = false;
      }
    }
  }
  std::cout << "largest difference from the scalar reference: " << worst << "\n";
  return valid ? 0 : 1;
}
