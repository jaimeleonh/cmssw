#ifndef PhysicsTools_ParticleTransformerAlpaka_graph_GraphShared_h
#define PhysicsTools_ParticleTransformerAlpaka_graph_GraphShared_h

// The per-jet block-shared state of the graph interpreter.
//
// Block shared memory is sized when the kernel is compiled, so the interpreter
// cannot grow an arena to fit whatever model it is handed.  Instead it
// publishes a fixed arena here, and the exporter's allocator packs the model's
// activations into it and refuses -- at export time, with the offending buffer
// named -- any model that does not fit.  A model file therefore either runs or
// never gets written; there is no runtime allocation failure to handle.
//
// The budget is the 47 KiB Alpaka CPU block-shared arena, which is also
// comfortably inside the 48 kB CUDA static limit.  Raising kArenaFloats means
// rebuilding CMSSW *and* re-exporting with a matching --arena-floats, so the
// two are checked against each other when a model is loaded.

#include <cstdint>

namespace nn {

  // Activation arenas.  The float arena holds everything the graph computes;
  // the int8 arena is staging for the activations that feed the dp4a GEMMs.
  // Overridable so a build can be sized to a particular model.  The default is
  // what the shipped ParT needs; a smaller network reserving the same 47 KiB
  // gets far fewer blocks per SM than its arithmetic deserves, because shared
  // memory -- not registers or occupancy -- is what limits residency here.
#ifndef NN_ARENA_FLOATS
#define NN_ARENA_FLOATS 10560
#endif
#ifndef NN_ARENA_INT8
#define NN_ARENA_INT8 4096
#endif
  inline constexpr int kArenaFloats = NN_ARENA_FLOATS;  // 42,240 B by default
  inline constexpr int kArenaInt8 = NN_ARENA_INT8;      //  4,096 B by default

  // Threads per jet.  The interpreter is written for any block size -- every
  // loop is a grid-stride over lanes -- but this is what the producer launches.
  // Overridable for the same reason as the arena: the right block size depends
  // on the model.  A wide network wants 256 lanes; a small one spends most of
  // its time at the block-wide barrier that ends every op, and fewer lanes make
  // that barrier cheaper.
#ifndef NN_THREADS_PER_JET
#define NN_THREADS_PER_JET 256
#endif
  inline constexpr int kThreadsPerJet = NN_THREADS_PER_JET;

  // Lanes that take part in a LayerNorm reduction.
  inline constexpr int kReduceThreads = 128;

  // Upper bounds the interpreter imposes on a model.  ModelLoader checks a
  // file against each of them and refuses one that exceeds any, so the kernel
  // itself never has to bounds-check.
  inline constexpr int kMaxTokens = 32;      // padded sequence length
  inline constexpr int kMaxReduceRows = 32;  // rows one normalization may cover
  inline constexpr int kMaxOutputs = 16;     // values a model may publish

  // The per-jet state, parameterised by its capacities so the same interpreter
  // can run either one jet per block with the full arena, or several jets per
  // block with a small arena each (see warp mode below).  The capacities are
  // carried as members so the interpreter reads them from the type it is
  // handed rather than from the global constants.
  template <int Floats, int Int8, int Reduce>
  struct alignas(16) GraphState {
    static constexpr int kFloats = Floats;
    static constexpr int kInt8 = Int8;
    static constexpr int kReduce = Reduce;

    alignas(16) float arena[Floats];
    alignas(16) std::int8_t quant[Int8];

    // Two distinct partial arrays, so no barrier is needed between reading the
    // mean partials and writing the variance partials.
    float reduceA[Reduce];
    float reduceB[Reduce];
    float mean[kMaxReduceRows];
    float rstd[kMaxReduceRows];

    float mask[kMaxTokens];
    float output[kMaxOutputs];

    int tokens;    // padded active token count
    int active;    // last unmasked particle plus one
    int pairs;     // active * (active + 1) / 2
    int tileBase;  // first row of the tiled region currently executing
    int tileRows;  // rows in that tile
  };

  // One jet per block: the whole arena.  This is what every model that needs
  // most of the arena -- the ParT and tau networks -- runs with.
  using GraphShared = GraphState<kArenaFloats, kArenaInt8, kReduceThreads>;

  // 47 KiB is the Alpaka CPU back-end's default block-shared arena and keeps
  // the struct inside the 48 kB CUDA static limit.
  static_assert(sizeof(GraphShared) <= 47u * 1024u,
                "the per-jet shared state grew past its budget; lower kArenaFloats");

  // ---------------------------------------------------------------- warp mode
  //
  // One jet per warp, several jets per block.  A small network spends most of
  // its time at the barrier that ends every op: with one jet per 256-lane
  // block there are only a few hundred elements of work between barriers, and
  // the block-wide barrier is what the SM stalls on.  Running each jet in its
  // own warp replaces that with a warp-level sync, and puts eight jets' worth
  // of independent work in every block.
  //
  // It only applies to a model whose activations fit the small arena; the
  // ParT and tau networks use nearly all of the full one and stay one jet per
  // block.  ModelLoader measures what a model actually uses, so no re-export
  // is needed to run a small model this way.
  inline constexpr int kWarpLanes = 32;
  static_assert(kThreadsPerJet % kWarpLanes == 0, "a block must hold a whole number of warps");
  inline constexpr int kWarpJetsPerBlock = kThreadsPerJet / kWarpLanes;
  // Overridable like the arena above.  Besides sizing a build to a model, it
  // is how the warp path gets tested on the large networks: with one warp per
  // block (NN_THREADS_PER_JET=32) the full arena fits, and ParT's whole op set
  // -- attention, LayerNorm, the pair embedding -- runs through SubGroup.
#ifndef NN_WARP_ARENA_FLOATS
#define NN_WARP_ARENA_FLOATS 1152
#endif
#ifndef NN_WARP_ARENA_INT8
#define NN_WARP_ARENA_INT8 128
#endif
  inline constexpr int kWarpArenaFloats = NN_WARP_ARENA_FLOATS;  // 4,608 B by default
  inline constexpr int kWarpArenaInt8 = NN_WARP_ARENA_INT8;

  using GraphWarpState = GraphState<kWarpArenaFloats, kWarpArenaInt8, kWarpLanes>;

  struct GraphWarpBlock {
    GraphWarpState jet[kWarpJetsPerBlock];
  };

  static_assert(sizeof(GraphWarpBlock) <= 47u * 1024u,
                "the warp-mode block state grew past the shared-memory budget; lower "
                "kWarpArenaFloats or run fewer jets per block");

}  // namespace nn

#endif
