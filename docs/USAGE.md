# Using the pipeline

Operational reference: sources, modules, restart, outputs and the meaning of
each report field. The [README](../README.md) has the short version.

## One pipeline

```text
surface / strata + exact target and containing rational envelope
           |
candidate source (no mandatory automaton)
  +-- reciprocal Newton spectra: Python or bundled C++
  +-- verified residuals from a previous proof run
  +-- imported spectra: explicitly incomplete unless coverage is rebuilt
  +-- signed folding paths in a generated or saved automaton
           |
applicable optional modules: traces, Schur, Perron, mod-2, peripheral
orbits, Lefschetz, deck / invariant-anti-invariant coupling, cycle graphs
           |
unresolved necessary data / certified exclusions
           |
complete lower-bound certificate + separately certified attaining map
```

All surface families use `ttauto_pipeline.engine.run`. There are no dated
campaign imports, absolute workspace paths, server credentials or known-minimum
lookups in that engine. Sources feed one event loop, journal, resume mechanism
and stage accounting. The large-disc and genus-three arithmetic kernels are
ported into the package, not executed through historical server scripts.

## Install and build

Python 3.10+ is required. The Python computation uses only the standard library.
Native builds require CMake 3.16+ and a C++17 compiler. All required algorithm
and numerical-library sources are bundled; Git is not needed to build a source ZIP.
GCC and Clang are supported; Windows is tested with LLVM-MinGW and Ninja.

```sh
python -m pip install -e .
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

On Windows, add `-G Ninja` and supply `CMAKE_C_COMPILER`, `CMAKE_CXX_COMPILER`
and `CMAKE_MAKE_PROGRAM` when the toolchain is not on PATH. Put its runtime DLL
directory on PATH when running the original C++ examples/tests. The native
automaton exporter is linked statically under MinGW.

All generated native files go under `build/`. The local jlt source includes the
tested CSparse portability and allocation fixes directly. The Python ribbon
builder and geometric certificates can run without a native build.

## Run, switch modules, resume, verify

```sh
ttauto modules
ttauto run experiments/proofs/disc3-minimum.json --output runs/d3-minimum
ttauto run experiments/proofs/genus3-bulk21.json --output runs/g3-lower
ttauto verify runs/g3-lower --replay-search
ttauto run experiments/proofs/disc12.json --output runs/d12 --events 1000
ttauto run experiments/proofs/disc12.json --output runs/d12 --resume --events 100000

ttauto run experiments/disc3.json --output runs/disc3
ttauto run experiments/disc3-native.json --output runs/disc3-native
ttauto verify runs/disc3 --replay-search

ttauto run experiments/closed-g2.json --output runs/g2 --events 10000
ttauto run experiments/closed-g2.json --output runs/g2 --resume --events 1000000
ttauto summary runs/g2

ttauto run experiments/disc3.json --output runs/no-cycle --disable cycle
ttauto run experiments/closed-g3-witness.json --output runs/g3-witness
```

`python -m ttauto_pipeline` is equivalent to `ttauto`. With no installation, set
`PYTHONPATH=python`; this also avoids package downloads when working offline.
The `native` source runs `build/bin/ttauto_pipeline_export`, compiled from this
repository. The old `upstream` source label is accepted as an alias for `native`.
The event budget is a **total** budget, including work already
committed. Increase it on resume. Changing a bound, pruning policy, mathematical
input, executable, automaton or Python source invalidates the checkpoint.
Use another output directory for a different experiment.

`disc3.json` certifies its two surviving maps using the exact homology action on
the branched torus cover, with no external process. `disc3-audit.json` selects
audit-only mode and deliberately leaves those maps unresolved.

The built-in geometry currently certifies the D3 four-one-prong pillowcase case
and the supported closed-surface canonical carriers. Other punctured-disc
survivors remain **unresolved**: their search and optional pruning still run,
but this version does not implement a general disc mapping-class classifier.
Removing the former classifier does not transfer its conclusions to these routines.

## Strata and campaigns

```sh
ttauto catalog disc 6 --one-prong --upper 2 --output configs/disc6
ttauto catalog disc 12 --one-prong --upper 313/250 --output configs/disc12
ttauto catalog disc 22 --upper 1131/1000 --output configs/disc22-all-markings
ttauto catalog closed 3 --upper 3/2 --output configs/genus3
ttauto campaign configs/disc6 --output runs/disc6 --events 10000
ttauto campaign configs/disc6 --output runs/disc6 --resume --events 1000000
```

Catalog entries default to built-in certificates, with explicit unresolved
outcomes outside their supported scope. `--one-prong`
fixes all original interior punctures to one prong. Without it, puncture prongs
vary too. Catalogs satisfy the Euler relation but do not assert stratum existence.
Closed catalogs use auxiliary regular inner marks; a fixed regular mark does not
cover every map on the closed surface without an additional theorem.

Use `experiments/proofs/disc*.json` for the migrated large-disc proof route.
Those specifications enumerate homological spectra without building a graph.
They use the bundled `ttauto_trace_prefix` binary; change the source to
`{"kind":"newton"}` to use arbitrary-precision Python integers instead.
Large examples carry explicit event/model/cycle budgets and may remain incomplete.
The old top-level disc preview files still select explicit ribbon enumeration.

## Outputs and interpretation

Each run contains `config.json`, `run.sqlite3` and `result.json`. Folding runs
also contain `automaton.json`; spectral runs contain `candidate_space.json`.
The report gives actual counts at every stage: tested, passed, rejected and
inapplicable, along with geometric outcomes and rational root intervals.
Counts use first rejection, so they are additive within a stage's input cohort.
Paths are counted with multiplicity; these are **not** conjugacy-class counts.

`enumeration_completed` means exhaustion of the chosen source. For spectral
proofs, `candidate_superset_complete` additionally certifies coverage by the
implemented finite containing enumeration. An arbitrary imported list cannot
set this flag. `lower_bound_proved` also requires every candidate/stratum
occurrence to be excluded, under the recorded mathematical theorem dependencies.
`attained` means a particular signed map passed the geometric criterion under
the recorded theorem dependencies. Unsupported gate patterns and computation
failures never become mathematical exclusions.

Spectral proofs set `minimum_proved` only when a complete lower bound matches
a separately revalidated, in-scope geometric witness at exactly the target root.
Explicit fold-search reports retain `minimum_proved=false` without a separate
whole-surface coverage proof. No historical Hall classification is accepted as
a new geometric certificate. Large-disc lower proofs can run fully without a
disc mapping-class classifier; a new Hall-free attaining certificate is still
required to promote them to minimum claims.

See [the proof workflow and mathematical contracts](PROOF_PIPELINE.md) for
source schemas, residual continuation, native limits and the corrected deck action.
The old full distributed campaigns have not been rerun merely by porting their code.

See [the mathematical contracts](PIPELINE.md),
[configuration examples](../experiments), [tests](../python/tests/test_pipeline.py),
[validation results](VALIDATION.md), and [publication instructions](RELEASING.md).
