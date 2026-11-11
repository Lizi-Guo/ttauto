# ttauto: unified computational proof pipeline

Proofs for punctured discs and closed orientable surfaces from finite spectral
supersets or train-track automata, with exact arithmetic, optional obstructions,
signed geometric witnesses and
transactional restart. All algorithm implementations are included as ordinary
source files in this repository, including the C++ train-track core, numerical
utilities and CSparse. Geometry uses this project's Python certificate routines.

The C++ core is adapted from
[Jean-Luc Thiffeault and Erwan Lanneau's ttauto](https://github.com/jeanluct/ttauto),
revision `ca1b2ed030946cd1ad5dd3abe05a85a52f59b0e7`. That link records provenance;
no separate ttauto installation is used. See [the source map](docs/SOURCES.md)
and `LOCAL_PROVENANCE.json` for attribution and local changes.

See [recorded results by configuration](#recorded-results-by-configuration)
for the numerical results, historical counts and their current input specifications.

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

See [the proof workflow and mathematical contracts](docs/PROOF_PIPELINE.md) for
source schemas, residual continuation, native limits and the corrected deck action.
The old full distributed campaigns have not been rerun merely by porting their code.

See [the mathematical contracts](docs/PIPELINE.md),
[configuration examples](experiments), [tests](python/tests/test_pipeline.py),
[validation results](docs/VALIDATION.md), and [publication instructions](docs/RELEASING.md).

## Recorded results by configuration

This register records the project results available through **2026-11-03**, mapped
to the current configuration format on **2026-11-10**. It contains **64 entries**:
26 one-pronged disc strata, five genus-two strata, 27 genus-three strata and six
unrestricted punctured discs. Historical records establish 56 minima; eight
genus-three stratum minima remain open. The D26 result is included. Nonorientable
surface N22 results are excluded; D22 below denotes a punctured disc.
The mapping classes considered here preserve the orientation of the surface.

The [machine-readable register](docs/HISTORICAL_RESULTS.json) stores every row's
exact polynomial, saved rational root interval, historical proof route, source
identifiers and a `current_query_configuration`. Polynomial coefficients are in
descending degree order; the specified root is the largest real root. Values
below are rounded to 12 decimal places. Source paths and SHA256 hashes identify
the original project records; the full historical archives are not bundled.

**These are recorded historical results, not 64 completed runs of this version.**
No mathematical search or proof replay was performed to assemble this register.
The configurations express the same lower-bound questions in the current schema;
they do not import historical completion flags, automata, coverage theorems or
attaining certificates. In particular, historical Hall-dependent attainment
records remain identified as such and are not accepted by the current engine.

### Configuration profiles

All 64 saved queries use `schema: "ttauto-proof-v1"` and the portable
`source: {"kind": "newton"}`. No automaton is required to start them.

| Profile | Surface and stratum selection | Enabled pruning modules |
| --- | --- | --- |
| `disc-stratum` | `kind: disc`, `punctures: n`, `one_pronged: true`; explicit `n`, `interior_prongs`, `outer_prongs` | `trace`, `schur`, `perron`, `mod2`, `peripheral`, `lefschetz`, `deck`, `cycle` |
| `closed-stratum`, H | `kind: closed`, `genus: 2` or `3`; explicit `type: H`, `orders` | `trace`, `schur`, `perron`, `lefschetz` |
| `closed-stratum`, Q | Same, with `type: Q` | `trace`, `schur`, `perron`, `lefschetz`, `base_coupling` |
| `disc-all` | `kind: disc`, `punctures: n`, `one_pronged: false`, `strata: "all"` | Same modules as `disc-stratum` |

These are explicit module choices for new queries, not reconstructed settings
for every historical job. Each query uses the row's recorded lower-bound
polynomial as `target_polynomial` and a rational `upper` strictly above its root.
For an open row, the target is its lower bound, not its attained upper bound.
The starter limits are 1,000 events and 1,000,000 orbit models, with 200,000 cycle
nodes for discs or 500,000 base-catalog nodes for Q strata. The trace horizon is
32, or 192 for large discs. These limits are not historical exhaustion thresholds
and do not promise completion. Disabled modules are explicitly `false`.

To extract and validate one query without running a search, execute this Python
from the repository root, then the validation command:

```python
import json
from pathlib import Path

index = json.loads(Path("docs/HISTORICAL_RESULTS.json").read_text(encoding="utf-8"))
row = next(r for r in index["records"] if r["id"] == "G3-Q-4-2-1-1")
Path("recorded-query.json").write_text(
    json.dumps(row["current_query_configuration"], indent=2) + "\n",
    encoding="utf-8",
)
```

```sh
ttauto validate recorded-query.json
```

For subsequent computation, the bundled C++ spectral source can replace `newton`;
see [the large-disc example](experiments/proofs/disc12.json) for the
`newton_native` executable and batch settings. Historical folding searches belong
to `ttauto-experiment-v1`, using a supplied or generated automaton and separately
justified coverage. They are not recreated by these spectral query files.
In the proof schema, `target_polynomial` asks for exclusion **strictly below** the
target root. In the folding schema it is an **equality filter**; omit it when
searching all smaller folding candidates.

Homology defines the spectral candidate space; it is not a `homology` switch in
the proof schema. The folding schema does have a `homology` filter. The current
`deck` and `base_coupling` modules use the actual deck/invariant decomposition:
the second lift acts by `JA`, not by a universal substitution `A -> -A`.
Historical counts below retain their original scopes and are not recomputed
counts for the corrected implementation. Stronger historical cycle-growth,
sector-coverage and special geometric exclusions are not all supplied by the
generic modules. A surviving necessary spectrum still needs further exclusion
or a geometric realization certificate.

### One-pronged disc strata: n = 3, 4, 5, 6, 7

Profile: `disc-stratum`. All original interior punctures are one-pronged.
`I` lists interior singularity prong numbers and `O` is the outer prong number;
`s` is the zero-based historical stratum index. Equal-prong singularities and
punctures may be permuted. All 26 rows have a historically proved minimum.

| Record / configuration ID | n | I | O | Minimum |
| --- | --- | --- | --- | --- |
| `D3-s0` | 3 | `[]` | 1 | 2.618033988750 |
| `D4-s0` | 4 | `[]` | 2 | 2.618033988750 |
| `D4-s1` | 4 | `[3]` | 1 | 2.296630262887 |
| `D5-s0` | 5 | `[]` | 3 | 1.722083805739 |
| `D5-s1` | 5 | `[3]` | 2 | 1.722083805739 |
| `D5-s2` | 5 | `[4]` | 1 | 2.153721375542 |
| `D5-s3` | 5 | `[3, 3]` | 1 | 2.015357181280 |
| `D6-s0` | 6 | `[]` | 4 | 1.883203505914 |
| `D6-s1` | 6 | `[3]` | 3 | 1.839286755214 |
| `D6-s2` | 6 | `[4]` | 2 | 1.883203505914 |
| `D6-s3` | 6 | `[5]` | 1 | 2.081018996625 |
| `D6-s4` | 6 | `[3, 3]` | 2 | 2.081018996625 |
| `D6-s5` | 6 | `[3, 4]` | 1 | 1.883203505914 |
| `D6-s6` | 6 | `[3, 3, 3]` | 1 | 2.171126688697 |
| `D7-s0` | 7 | `[]` | 5 | 1.556030191323 |
| `D7-s1` | 7 | `[3]` | 4 | 1.465571231877 |
| `D7-s2` | 7 | `[4]` | 3 | 1.465571231877 |
| `D7-s3` | 7 | `[5]` | 2 | 1.556030191323 |
| `D7-s4` | 7 | `[6]` | 1 | 2.042490533941 |
| `D7-s5` | 7 | `[3, 3]` | 3 | 1.610941783098 |
| `D7-s6` | 7 | `[3, 4]` | 2 | 2.475412101314 |
| `D7-s7` | 7 | `[3, 5]` | 1 | 1.809789333466 |
| `D7-s8` | 7 | `[4, 4]` | 1 | 1.754877666247 |
| `D7-s9` | 7 | `[3, 3, 3]` | 2 | 1.610941783098 |
| `D7-s10` | 7 | `[3, 3, 4]` | 1 | 2.042490533941 |
| `D7-s11` | 7 | `[3, 3, 3, 3]` | 1 | 2.025975350174 |

The historical family combined homological/orbit obstructions, reductions of
markings and finite folding searches for the residual cases, followed by
separate attainment checks. It was not one identical sequence of filters for
every row. Representative saved counts are:

| Case | Recorded residual reduction or search | Final outcome |
| --- | --- | --- |
| `D6-s3` | One lower residual; 20,562 folding prefixes | Residual excluded |
| `D6-s6` | 39 spectra after Lefschetz; mod-2 removes 7, strict-root filtering removes 5, leaving 27; 25 search groups, 15,700,700 prefixes, 341,299 closed paths | Three matches, all reducible; no smaller realization |
| `D7-s6` | Last three stages: 137,369,974 / 26,857,999 / 31,454,172 prefixes; each uses 99 reflection representatives | Zero matches in each stage |
| `D7-s11` | Phases 6-24 exhaust 21 residual spectra; final phase: 46 representatives and 2,809,458 prefixes | Zero matches in final phase; all residuals closed |

The first `D7-s6` prefix count includes donor trees and is not a count of globally
distinct paths. The difficult upper witnesses `D6-s6`, `D7-s6` and `D7-s11` have
independent geometric certificates. Earlier Hall-based witnesses retain their
historical dependency. Sources: `disc6-final`, `small-disc-26-completion`,
`small-disc-26-addendum`, `archived-S01` in the register.

### Closed orientable surfaces: genus 2

Profile: `closed-stratum`. `H` lists Abelian zero orders and `Q` lists zero orders
of non-square quadratic differentials. These denote the whole indicated stratum,
including its components.
All five minima are historically proved. The old proof uses the hyperelliptic
quotient and complete disc automata; the new spectral configurations are alternate
lower-bound questions and do not replace that coverage argument.

| Record / configuration ID | Stratum | Minimum | Historical quotient cases |
| --- | --- | --- | --- |
| `G2-H-2` | H(2) | 1.722083805739 | `(5,0)` |
| `G2-H-1-1` | H(1,1) | 1.883203505914 | `(6,0)` |
| `G2-Q-2-2` | Q(2,2) | 2.618033988750 | `(3,0)`, `(4,0)` |
| `G2-Q-2-1-1` | Q(2,1,1) | 1.722083805739 | `(5,0)` + regular branch point |
| `G2-Q-1-1-1-1` | Q(1,1,1,1) | 1.839286755214 | `(5,3)`, `(6,1)`, `(6,4)` |

The quotient's `(n,s)` refers to the disc indexing above. `Q(2,1,1)` also requires
the regular branch-point condition. Historical homology certified the upper
witnesses, and Lefschetz arguments supplied regular periodic-point information;
the lower proof did not begin with an enumeration of homological polynomials.

| Quotient `(n,s)` | Start representatives | Prefixes | Closures | Primitive checks |
| --- | ---: | ---: | ---: | ---: |
| `(3,0)` | 1 | 13 | 12 | 2 |
| `(4,0)` | 4 | 1,310 | 1,436 | 74 |
| `(5,3)` | 9 | 10,239 | 1,072 | 228 |
| `(6,0)` | 49 | 641,630 | 58,155 | 3,948 |
| `(6,1)` | 138 | 6,572,375 | 309,981 | 12,836 |
| `(6,4)` | 110 | 7,628,379 | 191,074 | 14,940 |

The `(5,0)` case reuses an earlier certificate on an 11-state automaton. Columns
count different objects, so their differences are not rejection counts. Source:
`genus2-proof`. The example [closed-g2.json](experiments/closed-g2.json) searches
an auxiliary regular-fixed-point sector; it alone does not certify these five
whole-stratum results.

### Closed orientable surfaces: genus 3

Profile: `closed-stratum`. Nineteen strata have historically proved minima:

| Record / configuration ID | Stratum | Minimum | Historical lower route |
| --- | --- | --- | --- |
| `G3-H-1-1-1-1` | H(1,1,1,1) | 1.401268367940 | L |
| `G3-Q-5-1-1-1` | Q(5,1,1,1) | 1.401268367940 | A |
| `G3-H-2-1-1` | H(2,1,1) | 1.465571231877 | A |
| `G3-Q-6-1-1` | Q(6,1,1) | 1.465571231877 | A |
| `G3-H-4` | H(4) | 1.506135679554 | M + S |
| `G3-H-3-1` | H(3,1) | 1.506135679554 | A |
| `G3-Q-3-3-2` | Q(3,3,2) | 1.506135679554 | A |
| `G3-Q-5-3` | Q(5,3) | 1.582347183680 | A |
| `G3-H-2-2` | H(2,2) | 1.635573129922 | M + S |
| `G3-Q-4-3-1` | Q(4,3,1) | 1.635573129922 | A |
| `G3-Q-4-2-1-1` | Q(4,2,1,1) | 1.635573129922 | S (joint homology) |
| `G3-Q-3-2-2-1` | Q(3,2,2,1) | 1.635573129922 | A |
| `G3-Q-4-4` | Q(4,4) | 1.722083805739 | S (rank four) |
| `G3-Q-4-2-2` | Q(4,2,2) | 1.722083805739 | S + independent A |
| `G3-Q-6-2` | Q(6,2) | 1.883203505914 | A |
| `G3-Q-5-2-1` | Q(5,2,1) | 1.883203505914 | P |
| `G3-Q-2-2-2-2` | Q(2,2,2,2) | 1.883203505914 | M (two sectors + low matrices) |
| `G3-Q-7-1` | Q(7,1) | 1.987793166846 | S + A (spectral gap) |
| `G3-Q-8` | Q(8) | 2.081018996625 | M + P (fixed / period two) |

Here **A** means complete stratum automaton search; **P** complete projected
automaton search; **S** spectral/homological exclusion; **M** mixed periodic-point
or sector arguments; **L** a published lower bound. Each minimum also has a
separate historical attaining construction. `H(4)` and `H(2,2)` are whole-stratum
minima; their hyperelliptic components have the larger minima 1.556030191323 and
1.839286755214, respectively.

The remaining eight rows record proved lower bounds and attained upper bounds;
their exact stratum minima are open:

| Record / configuration ID | Stratum | Lower bound | Attained upper bound |
| --- | --- | --- | --- |
| `G3-Q-4-1-1-1-1` | Q(4,1,1,1,1) | 1.401268367940 | 1.610941783098 |
| `G3-Q-3-3-1-1` | Q(3,3,1,1) | 1.401268367940 | 1.809789333466 |
| `G3-Q-3-2-1-1-1` | Q(3,2,1,1,1) | 1.401268367940 | 2.618033988750 |
| `G3-Q-3-1-1-1-1-1` | Q(3,1,1,1,1,1) | 1.401268367940 | 2.618033988750 |
| `G3-Q-2-2-2-1-1` | Q(2,2,2,1,1) | 1.401268367940 | 1.722083805739 |
| `G3-Q-2-2-1-1-1-1` | Q(2,2,1,1,1,1) | 1.401268367940 | 2.042490533941 |
| `G3-Q-2-1-1-1-1-1-1` | Q(2,1,1,1,1,1,1) | 1.401268367940 | 1.610941783098 |
| `G3-Q-1-1-1-1-1-1-1-1` | Q(1,1,1,1,1,1,1,1) | 1.363878022887 | 1.839286755214 |

Historical counts illustrating the different proof routes:

| Case | Saved counts | Scope |
| --- | --- | --- |
| `H(3,1)` | 936 states, 6,048 arcs; 12,474,712 prefixes, 233,804 closures, 4,504 primitive checks, 85 spectra | Complete stratum automaton |
| `H(2,1,1)` | 14,784 states, 119,760 arcs; 190,150,702 prefixes, 15,881,744 closures, 145,406 primitive checks | Complete stratum automaton |
| `Q(5,1,1,1)` | 8,871,800 states, 74,664,240 arcs; 3,736,081 states excluded by coordinate tests; 23,672 primitive closures, 214 spectra | Complete stratum automaton |
| `Q(5,2,1)` | 173,025 original states, 879,210 arcs; 96,939 projected states, 386 ranges, 3,526 raw matrices | Complete projected search |
| `Q(4,2,1,1)` | 1,107 coefficient triples -> 141 disk-bounded spectra -> 16 positive expanding spectra; 432 singularity actions -> 54 joint sequences | Joint invariant/anti-invariant exclusion leaves no smaller candidate; no full automaton |
| `Q(2,2,2,2)` | First sector: 75,194 states, 1,175 ranges, 4,674 final replayed matrices | Sector count only; completion also uses the other sector, low-matrix exclusions and an attained V4 cover |

The whole genus-three minimum is still reduced to two possible values:

$$
\delta_3\in\{\beta,c\},\qquad
\beta\approx1.363878022887,\quad c\approx1.401268367940,
$$

where beta is the largest real root of `x^10 - x^8 - 2*x^5 - x^2 + 1`
and c of `x^6 - x^4 - x^3 - x^2 + 1`. Twenty-six strata have lower bound c.
Only `Q(1,1,1,1,1,1,1,1)` can still realize beta; its saved remaining joint model
is `[1,0]` at the recorded horizon 16. Attainment of c is established, while
attainment of beta is unresolved.

The historical bulk lower-bound proof for 21 strata used finite spectra, not a
whole genus-three automaton: rank 4 used 81 coefficient possibilities; rank 6
used 15,093 -> 75 small-house spectra -> three positive candidates in the
containing interval; rank 8 used 8,719 -> 120 -> 25 -> 20 below c. Boundary,
root-rotation, dominance and joint orbit obstructions complete those exclusions;
the rank-eight joint stage has 75 occurrences and zero survivors. These are
different rank cohorts, not one subtractable chain.

The principal-stratum P12 residual was separately closed by
4,959 unsigned objects -> 2,332 signed pairs -> 21,748 word families ->
19,704 coarse keys -> 39,408 key rotations -> 1,226,440 cyclic orders ->
zero coarse surfaces. This specialized geometric search is not a generic
current `cycle` run. Sources include `genus3-bulk21-proof`, `genus3-two-values`,
`genus3-current-state`, `archived-S07`, `archived-S08` and the individual
genus-three proof records listed in the register. None of these counts asserts
that a complete automaton for all genus-three maps was built.

### Unrestricted punctured discs: n = 10, 12, 14, 18, 22, 26

Profile: `disc-all`. These are whole-disc minima, allowing the original marked
points to have different prong types. Every row has a historically proved
minimum. The old lower proofs enumerated homological spectra, applied exact
algebraic and orbit obstructions and closed residuals with arithmetic or cycle
graph certificates. They did not construct complete train-track automata.

| Record / configuration ID | Minimum | Polynomial defining the minimum |
| --- | --- | --- |
| `D10` | 1.343719995611 | `x^9 - 2*x^5 - 2*x^4 + 1` |
| `D12` | 1.251413186052 | `x^12 - 2*x^7 - 2*x^5 + 1` |
| `D14` | 1.222573998205 | `x^14 - 2*x^9 - 2*x^5 + 1` |
| `D18` | 1.164318329674 | `x^18 - 2*x^11 - 2*x^7 + 1` |
| `D22` | 1.130709350883 | `x^22 - 2*x^13 - 2*x^9 + 1` |
| `D26` | 1.108686289468 | `x^26 - 2*x^15 - 2*x^11 + 1` |

Let `rho_n = sigma_1 ... sigma_(n-1)`. The historical attaining braids are
`(sigma_1^2 sigma_2 ... sigma_9)^2 sigma_9^-2` for D10, and
`rho_n^k sigma_(n-1)^-2` with `(n,k) = (12,5), (14,3), (18,5), (22,5), (26,7)`.
Their historical classification depends on Hall; a new Hall-free upper
certificate is not supplied by the spectral configurations.

| Disc | Historical strata | Retained entries | After algebraic | After topological | Final unresolved |
| --- | --- | --- | --- | --- | --- |
| `D10` | 120 | 358,430 | - | 3 | 0 |
| `D12` | 97 | 156,012 | 17 | 0 | 0 |
| `D14` | 195 | 85,639 | 938 | 2 | 0 |
| `D18` | 684 | 220,107 | 1,590 | 7 | 0 |
| `D22` | 2,087 | 1,588,303 | 2,383 | 19 | 0 |
| `D26` | 5,763 | 15,809,615 | - | 7 | 0 |

Here "retained" means native retained polynomial entries; "algebraic" counts
spectra; "topological" and "final" count spectrum/stratum occurrences. A dash
means the consulted summary does not separately report that count. Differences
between these columns must not be called stage rejection counts. The old D10
enumeration covers one-pronged cores `m = 3,...,10`. D12-D26 first exclude `m<n`
using a published lower bound and enumerate the `m=n` strata shown above. Current
`disc-all` configurations instead cover all `m=3,...,n`, so their stratum and
candidate counts will differ.

The last historical exclusions are recorded as follows:

| Disc | Closure of the topological residuals |
| --- | --- |
| D10 | Three occurrences: one negative periodic-orbit count, one mod-2 kernel determinant obstruction (`16` does not divide `8`), one cycle search with 419 nodes and 10 matches; zero connected realizations |
| D12 | All 17 algebraic spectra excluded by topological tests |
| D14 | Two occurrences: one opposite-root obstruction; one cycle search with 1,936 nodes and 10 matches, zero connected realizations |
| D18 | Seven occurrences: four opposite-root exclusions; the other three spectra give six real transition polynomials, excluded by six cycle searches totaling 11,963 nodes |
| D22 | Nineteen occurrences: 13 opposite-root exclusions; six residual branches closed by 192 certificates |
| D26 | Seven occurrences: five already excluded; the last two give three real transition polynomials. Cycle searches visit 47,267 / 63,712 / 239,734 nodes and produce 272 / 272 / 1,926 matches, all with zero connected realizations |

D26 also uses stronger cycle-pair, twin and induced-growth obstructions than the
current basic cycle kernel. The saved D18 and D22 first-rejection algebraic counts
are preserved in `large_disc_counts` in the register. Sources: `D10-proof`,
`D12-D14-proof`, `D18-proof`, `D22-final`, `disc-all5-final`,
`disc-all5-archive`, `D26-last-two-archive` and the algebraic summary snapshots.

### Completed checks of the current implementation

These are the **already saved 0.2.0 validation runs**, distinct from the historical
tables and from the unexecuted configurations in the register. See
[the validation record](docs/PROOF_VALIDATION_0.2.0.json).

| Saved run / configuration | Recorded accounting | Certified scope |
| --- | --- | --- |
| [Genus-three bulk 21](experiments/proofs/genus3-bulk21.json) | 11,715 spectra; trace rejects 11,228, Schur rejects 46, domain rejects 413, Perron rejects 11, leaving 17 spectra; 46 stratum occurrences, Lefschetz rejects 45 and base coupling rejects 1 | Complete lower bound c for these 21 strata; no attaining witness, `minimum_proved=false` |
| [D3 minimum](experiments/proofs/disc3-minimum.json) | Seven spectra excluded from the strict-below-target domain; separate small signed-map witness | `minimum_proved=true` |
| [D12 preview](experiments/proofs/disc12.json) | 50 committed events, including 47 processed spectra; paused at budget | Incomplete; `lower_bound_proved=false` |

The documentation update only validates configurations and cross-checks the saved
numbers. It does not rerun these searches or extend their certified scope.

## License and attribution

GPL-3.0-or-later, preserving upstream authorship and notices. `extern/jlt` and
its CSparse dependency retain their own licenses.
`dependencies.lock.json` records the included native source revisions.
