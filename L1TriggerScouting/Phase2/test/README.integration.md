# CMSSW 20.0.0-pre1 integration

This area starts from `jaimeleonh/cmssw:20_0_0_pre1_alpaka` at
`ad4a0136d7bde3a0cf0cd2196234df6907a8ac4c` and ports the clustering changes in
`jonjanik/cmssw:P2_jet_clustering` at
`2384008cf5505abf3bd55136368ea506f9a74a20`.
The clustering-specific delta is the Phase2 subtree difference from
`de721132de59cdee2eb3e7bfb12405abd6381982` to that source tip.
The release uses `SCRAM_ARCH=el9_amd64_gcc13`.

## Conflicts and resolutions

The three-way application reported conflicts in these five files, relative to
`L1TriggerScouting/Phase2/`:

| File | Resolution |
| --- | --- |
| `plugins/alpaka/L1TScPhase2SCJets.cc` | Use the source's SeededCone / SCNMSW / LinkTree dispatcher and validation, with the target's `BxLookupDevice`, queue-first collection constructors, and optional `minConstituents` cut. |
| `plugins/alpaka/L1TScPhase2SCJetsKernels.h` | Adapt source signatures to the target product types and retain the constituent-cut setting. |
| `plugins/alpaka/L1TScPhase2SCJetsKernels.dev.cc` | Port the source launches and finalizers to the target's SoABlocks BX and association-map views. Retain the target's constituent filter and use the shared finalizer when filtering SCNMSW outputs. |
| `plugins/alpaka/radixSort64.h` | The versions differed only in whitespace; retain the target version. |
| `test/runScoutingPhase2PF_cfg.py` | Use the source's algorithm names, defaults, rounds and synchronization options. Retain the target's vertex unpacker and TorchScript tagger behind `runJetTagging=True`, with corrected model-path handling and a local FastTimerService output filename. |

Additional compatibility work:

- Add the source's `int is_seed` column to
  `DataFormats/L1ScoutingSoA/interface/ClustersSoA.h`. LinkTree uses it to identify
  roots; the other algorithms populate it as well.
- Adapt the new `JetOutput.h`, `SeededCone.h`, `SCNMSW.h` and `LinkTree.h` to the
  target's SoABlocks accessors without replacing its BX or association-map schemas.
- Avoid a zero-block CUDA launch when the shared finalizer receives no jets.
  The regression test exercises this with an empty event.

No merge conflicts remain. The changes are staged locally, without a commit.

## Running

From the release's `src` directory:

```sh
export SCRAM_ARCH=el9_amd64_gcc13
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval "$(scram runtime -sh)"
cmsRun L1TriggerScouting/Phase2/test/testJetClustering_cfg.py backend=serial_sync
cmsRun L1TriggerScouting/Phase2/test/testJetClustering_cfg.py backend=cuda_async
```

The PF job accepts `jetAlgo=SeededCone`, `jetAlgo=SCNMSW` or `jetAlgo=LinkTree`.
See `README.jet-clustering.md` for the algorithm and throughput options.
`minConstituents` defaults to zero and can be used with all three algorithms.
`runJetTagging` defaults to false so clustering jobs do not require a model;
use `runJetTagging=True jetTagModel=<CMSSW-relative model path>` to enable the
retained tagger. Model weights and custom-operator setup are external inputs;
jet-tagger inference was not part of this clustering validation.

## Validation

- Full `scram b -j 6` build, including serial and CUDA clustering and jet-tagging
  plugins and data-format dictionaries.
- Synthetic CPU and CUDA checks over 21 combinations of algorithms, SCNMSW
  exponents/rounds, the SeededCone jet cap, and constituent cuts. Inputs include
  a phi branch cut, tied pT, a BX with 132 candidates, an empty BX and an empty
  event. Checks cover BX IDs, counts, associations, constituent ordering and
  scalar jet centroids.
- Full RAW unpacking-to-jets runs for all three algorithms on serial and CUDA,
  each over two orbits (7,128 bunch crossings) from the existing run-42 PUPPI
  input, including OrbitFlatTable ROOT output.
- The CPU/GPU ROOT comparison found identical jet constituent memberships and
  exact pT agreement for all three algorithms. SCNMSW and LinkTree eta/phi agreed
  exactly; SeededCone differed by at most 2.3842e-7 in eta and 1.1921e-7 in phi.

Validation artifacts and scripts are in the release's `validation/` directory.
Tests were run outside the execution sandbox: even an empty CMSSW job crashed
at startup inside it, while the same environment outside exposed two NVIDIA L4s.
The target branch still emits its pre-existing `Invalid tool PerfTools/Perfetto`
warning during builds; the build succeeds.
