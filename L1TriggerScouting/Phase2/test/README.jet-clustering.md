# Phase-2 jet clustering

The Alpaka producer supports three algorithms: `SeededCone`, `SCNMSW`, and
`LinkTree`. Their kernels are in the correspondingly named headers in
`plugins/alpaka/`; `JetOutput.h` contains shared output construction, and
`L1TScPhase2SCJetsKernels.dev.cc` handles allocation and launches.

All three use scalar constituent-pT sums and pT-weighted eta/phi offsets, with
wrapped phi differences. SCNMSW recomputes jet properties from the final assigned
constituents. SeededCone accumulates relative to its seed, LinkTree relative to
its root, and SCNMSW relative to its protojet axis. These references matter for
constituents spanning a phi branch cut; this is not four-vector recombination.

## Configuration

Use `runScoutingPhase2PF_cfg.py` with the usual input options and one of:

| Selection | Options |
| --- | --- |
| CPU FastJet AK4 | `run=ak4` |
| CPU seeded cone | `run=sc4` |
| GPU seeded cone, at most 16 jets/BX | `run=sc4Alpaka jetAlgo=SeededCone njets=16` |
| GPU seeded cone, study's unconstrained setting | `run=sc4Alpaka jetAlgo=SeededCone njets=1024` |
| SCNMS4 plot configuration | `run=sc4Alpaka jetAlgo=SCNMSW alphaSeed=0 maxRounds=1` |
| SCNMSW4 plot configuration | `run=sc4Alpaka jetAlgo=SCNMSW alphaSeed=2 maxRounds=1` |
| LinkTree | `run=sc4Alpaka jetAlgo=LinkTree` |

Use `backend=cuda_async numThreads=1 numFwkStreams=1` for the GPU study setup.
Default radii are `jetR=0.4`, `RSeed=0.4`, `RCen=0.4`, `RClu=0.4`, and
`RLink=0.36`. `minSeedPt` affects LinkTree and the CPU seeded cone.

SCNMSW retains configurable rounds: `maxRounds=2` runs up to two rounds;
`maxRounds=None` (or `-1`) runs until no unassigned candidates remain. The
one-round implementation uses a compact shared-memory path for BXs with up to
128 candidates and an indexed fallback for larger BXs. Alpha zero has a
distance-only assignment specialization; other exponents use the weighted
assignment. There is no separate SCNMS algorithm or alpha-one specialization.

The former experimental algorithm names and tuning switches were removed.
Historical job commands must be translated to the options above. CPU AK4 is
independent of these Alpaka implementations.

## Throughput

After `cmsenv`, run the five GPU plot cases over the complete TTbar file with:

```sh
python3 L1TriggerScouting/Phase2/test/throughputPhase2.py \
  --cases SC4_16 SC4_unconstrained SCNMS4 SCNMSW4 LinkTree \
  --orbits -1 --threads 1 --streams 1 --repeats 5 \
  --synchronize False --output-dir throughput-async
```

Run separately with `--synchronize True` and another output directory to measure
completed jet-module work. The synchronization flag adds module-boundary waits
at both unpacking and clustering; use the asynchronous run for end-to-end
throughput and the synchronized run for jet-module throughput. Existing waits
needed to allocate variable-sized outputs remain in both modes. Module timing
includes output compaction and constituent-map construction, not just clustering
kernels. The helper defaults to synchronized timing, refuses to overwrite jobs,
and verifies processed orbit/BX counts. Add `AK4` to `--cases` for the CPU reference.

The finalized study remains in
`throughput-study/single-stream-alpha0-20261104/`; cleanup validation and smoke
runs are in `throughput-study/algorithm-cleanup-20261104/`. Smoke timings do not
replace the study's repeated measurements.

## CMSSW 20.0.0-pre1 port

See `README.integration.md` for the target/source revisions, conflict resolutions
and validation. This port also supports the target's `minConstituents` cut and
optional `runJetTagging=True` workflow. Clustering and throughput jobs default to
`runJetTagging=False`. The historical study directories above belong to the
source branch's study environment and are not included in this release area.
