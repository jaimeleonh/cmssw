#ifndef PhysicsTools_ParticleTransformerAlpaka_ModelConstants_h
#define PhysicsTools_ParticleTransformerAlpaka_ModelConstants_h

#include <cstddef>

namespace part {

  // Fixed graph dimensions of the exported Weaver ParticleTransformer.
  inline constexpr int kMaxParticles = 16;   // padded input sequence length
  inline constexpr int kInputFeatures = 15;  // pf_features from b_kinadd.yaml
  inline constexpr int kFourVector = 4;      // px, py, pz, e
  inline constexpr int kEmbedDim = 128;
  inline constexpr int kHeads = 8;
  inline constexpr int kHeadDim = kEmbedDim / kHeads;  // 16
  inline constexpr int kQkvDim = 3 * kEmbedDim;        // 384
  inline constexpr int kFfnDim = 512;
  inline constexpr int kPairInputs = 4;
  inline constexpr int kPairHidden = 64;
  inline constexpr int kParticleBlocks = 8;
  inline constexpr int kClassBlocks = 2;
  inline constexpr int kClasses = 2;

  // Upper-triangular pair storage.  The pairwise features are symmetric under
  // i <-> j, so only i <= j is evaluated and stored.
  inline constexpr int kMaxPairs = kMaxParticles * (kMaxParticles + 1) / 2;  // 136

  // Cooperative kernel tiling.
  inline constexpr int kThreadsPerJet = 256;  // block size of the inference kernel
  // Upper bounds of the register tile searched by tokenBlockFor/columnBlockFor.
  inline constexpr int kTokenBlockMax = 8;   // register-blocked rows of a GEMM
  inline constexpr int kColumnBlockMax = 4;  // register-blocked columns of a GEMM
  // Granularity the token count is padded to, so a tile never straddles the
  // end of a buffer.
  inline constexpr int kTokenBlock = 4;
  inline constexpr int kFfnTokenTile = 8;     // tokens whose 512-wide hidden state is live
  inline constexpr int kPairTile = 32;        // pairs evaluated per pair-embedding pass
  inline constexpr int kReduceThreads = 128;  // threads contributing to a LayerNorm reduction

  // Scratch shared by the transformer blocks, the class blocks, the pair
  // embedding and the feed-forward tiles.  Sized by the largest consumer,
  // which is the 16-token QKV projection.
  inline constexpr int kBigScratch = kMaxParticles * kQkvDim;  // 6144
  // Only one feed-forward tile of normalised tokens is live at a time.
  inline constexpr int kNormScratch = kFfnTokenTile * kEmbedDim;
  // INT8 staging for the activations that feed the integer GEMMs, in bytes.
  inline constexpr int kQuantScratch = kFfnTokenTile * kFfnDim;

  static_assert(kBigScratch >= kFfnTokenTile * (kFfnDim + kEmbedDim));
  static_assert(kBigScratch >= (kMaxParticles + 1) * 2 * kEmbedDim + 3 * kEmbedDim + kHeads * (kMaxParticles + 1));
  static_assert(kBigScratch >= 2 * kPairTile * kPairHidden);
  static_assert(kBigScratch >= kMaxParticles * kEmbedDim + kFfnTokenTile * kFfnDim);
  static_assert(kNormScratch >= kMaxParticles * kInputFeatures);
  // The INT8 staging area must hold the widest quantized activation tile: the
  // feed-forward hidden state, one 8-token tile of 512 channels.
  static_assert(kQuantScratch >= kFfnTokenTile * kFfnDim);
  static_assert(kQuantScratch >= kMaxParticles * kEmbedDim);
  static_assert(kQuantScratch % 16 == 0);
  static_assert(kFfnTokenTile % kTokenBlock == 0);
  static_assert(kFfnTokenTile % kTokenBlockMax == 0);
  static_assert(kPairTile % kTokenBlockMax == 0);

}  // namespace part

#endif
