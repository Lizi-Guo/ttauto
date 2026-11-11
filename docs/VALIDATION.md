# Validation records

## 0.2.0 unified spectral and folding proofs (2026-11-10)

The full CMake build succeeds, including the bundled exact Newton enumerator.
All 20 CTest targets pass. After adding an independent native-versus-Python
candidate-set check, the Python CTest target was rerun: all 115 tests pass.
The GitHub workflow has been updated but has not been run remotely.

| Fresh run | Spectral candidates | Stratum occurrences reaching topology | Unresolved | Lower bound | Attainment |
|---|---:|---:|---:|---|---|
| g3 bulk21, target root of x^6-x^4-x^3-x^2+1 | 11,715 | 46 | 0 | proved for these 21 strata | not supplied |
| D3, target root of x^2-3x+1 | 7 | 0 | 0 | proved | built-in witness verified |
| D12 native preview, 50 committed events | 47 | 0 | 0 observed | incomplete | not supplied |

All three runs pass independent record verification and full source replay.
The g3 run constructs no automaton: trace rejects 11,228 spectra, Schur 46,
the target-domain check 413, and Perron 11. The remaining 17 spectra produce
46 stratum occurrences: Lefschetz rejects 45 and invariant/anti-invariant
base coupling rejects the last one. This establishes the common lower bound
approximately 1.40126836794 on the specified 21 strata, not a global genus-three
minimum or attainment in each stratum.

The D12 preview is explicitly `paused_budget`; an empty observed residual set
does not establish coverage. Its three native batches visit 485 nodes, with
375 parity-admissible completions, 300 trace prunes and 75 retained spectra;
28 retained spectra remain queued. Large-disc historical campaigns were not
fully rerun, and their Hall-dependent upper classifications are not imported
as independent certificates.

The lift regression uses an actual Q(8) cover: the first four traces of A are
`[0, 4, 6, 20]`, and those of the other lift JA are `[-2, 4, -14, 20]`.
Replacing JA by -A would incorrectly give `[0, 4, -6, 20]`. Tests also cover
negative-return orbit inversion, coupled branch permutations, incomplete
budgets, residual coverage, evidence tampering, native frontier resume and
the D14 abstract cycle control (1,936 nodes, 3,008 capacity prunes, 10 matching
graphs and no connected survivor).

Machine-readable fresh counts are in `PROOF_VALIDATION_0.2.0.json`.
The wheel is built from a fresh temporary source tree with an exact module
inventory check, so removed classifier files cannot reappear from build caches.

## 0.1.2 built-in geometry, without Hall (2026-11-10)

The source archive contains no Hall implementation, adapter or process client.
It builds in a fresh directory with Git absent from PATH. All 20 main CTest
targets pass, including 92 Python tests. The four core C++ executables and both
short examples pass. A complete D3 search and replay also pass with every
external process launch forbidden by the test harness.

| Fresh run | States | Arcs | Closed-path candidates | Attained | Unresolved | Best attained value |
|---|---:|---:|---:|---:|---:|---:|
| D3, U=8/3, built-in geometry | 1 | 2 | 20 | 2 | 0 | 2.618033988750 |
| D3, U=8/3, audit only | 1 | 2 | 20 | 0 | 2 | -- |
| g2, I2/O6, U=7/4 | 76 | 262 | 28,614 | 10 | 0 | 1.722083805739 |
| g3 Q(8) witness subgraph, depth 5 | 5 | 5 | 5 | 5 | 0 | 2.081018996625 |
| D4, U=3, depth 4, norm-only pruning | 4 | 14 | 362 | 0 | 38 | -- |

All candidate records in these runs pass independent replay. The D3, g3 and
bounded D4 runs also pass full search re-enumeration. The g2 run reproduces
408,954 proposed prefixes, 236,986 norm rejections and 8,970 cycle rejections.
The g3 fixture and D4 depth cap do not exhaust a whole-stratum search.

D4's 38 primitive, in-bound survivors deliberately remain unresolved; 240
candidates fail primitivity and 84 more fail the growth bound. No general
punctured-disc classifier is claimed. The supported D3 and closed-surface
certificates are described in `PIPELINE.md`. Higher-dimensional homology alone
does not replace a geometric proof. The large-disc campaigns were not rerun.

External geometry configurations are rejected. Replay checks every field in the
geometric evidence; inserting even an extra purported classifier report is
detected. Earlier Hall-based validations below are historical records only.

## 0.1.1 standalone source validation (2026-11-10)

A source archive was extracted into a separate tree with no `.git` or
`.gitmodules` files. CMake configured and built the complete project with Git
absent from PATH. All numerical and train-track implementations came from
the extracted source tree; no dependency checkout or download was used.

- All 20 main CTest targets pass, including all 90 Python tests with native integration.
- The four core C++ executables and both short examples pass in the same environment.
- All 17 portable jlt test targets pass, including CSparse and matrix arithmetic.
  The POSIX-shell command tests and Git-metadata tests were excluded from this
  Windows/no-Git check. Optional LAPACK, Matlab-library and Boost-timer tests were
  not built; those libraries are not used by the pipeline.
- A fresh D3 run with `source.kind = "native"` has 1 state, 2 arcs, 20 closed-path
  candidates and 2 attained paths, with value 2.618033988750. All 20 candidate
  records and the complete search replay pass verification.
- The source tree has no Git submodule entries. Source archives include the local
  algorithm implementations, dependency license files and per-file checksums.

The bundled jlt tests exposed a pre-existing template compilation error in
matrix compound addition/subtraction. Separate mutable and const iterator
declarations fix it; the existing arithmetic tests now pass.

The larger historical campaigns were not rerun for this packaging change.
The results below remain the dated 0.1.0 baseline, not new full reruns.

## 0.1.0 baseline (2026-11-08)

Tested locally on Windows with Python 3.12, Clang/LLVM-MinGW 20.1 and CMake 3.31.6.
The GitHub workflow additionally defines an Ubuntu build; that remote workflow
has not been run as part of this local preparation.

## Automated checks

- 89 Python tests pass, including native export/Hall integration.
- All 20 CTest targets pass (19 C++ tests and the Python suite).
- The four upstream core executables and both short examples pass, with test assertions enabled.
- AddressSanitizer checks pass for gate analysis, graph construction, bad-word
  pruning and native search.
- Each optional pruning module is compared against disabled pruning on a
  complete small search. The retained signed witnesses agree.
- Restart equals uninterrupted execution; changed policy/source identity is
  refused. Damaged records, modified proof claims and mismatched matrices are rejected.
- A simulated backend timeout rolls back work and is successfully resumed;
  it is not counted as an exclusion.
- Rooted colored canonicalization agrees with explicit symmetry-orbit enumeration
  on a finite independent control family.
- An allowlisted source archive was extracted to a separate directory. Both native
  adapters build from that snapshot, all 89 Python tests pass there, and a fresh
  D3 upstream run passes complete search replay without importing the original workspace.

## End-to-end checks

The counts below are path counts, not conjugacy classes. Bounds are inclusive.

| Experiment | States | Arcs | Closed paths with terminal symmetry | Cycle rejections | Primitive survivors after cycle | Geometrically attained paths | Best attained value |
|---|---:|---:|---:|---:|---:|---:|---:|
| D3, U=8/3, Python ribbon source | 1 | 2 | 20 | 6 | 2 | 2 | 2.618033988750 |
| D3, U=8/3, upstream C++ source | 1 | 2 | 20 | 6 | 2 | 2 | 2.618033988750 |
| g2, I2/O6, U=7/4 | 76 | 262 | 28,614 | 8,970 | 10 | 10 | 1.722083805739 |
| g3 Q(8) witness subgraph, depth 5 | 5 | 5 | 5 | 0 | 5 | 5 | 2.081018996625 |

The g2 search examines 408,954 proposed prefixes; the norm module rejects
236,986. All ten primitive survivors pass growth, homology, Lefschetz, Perron
and actual-stratum checks. Their real polynomial is
`x^4 - x^3 - x^2 - x + 1`.

Every candidate record in these runs was replayed, including all 28,614 g2
records. Both D3 runs and the g3 witness-subgraph run also passed full search
re-enumeration. The g3 witness polynomial is
`x^6 - x^5 - 3*x^4 + 3*x^2 + x - 1`; its reconstructed carrier has one
ten-prong singularity and a two-gate single-edge junction.

The generic builder separately generated the complete prescribed g3 I2/O10
state inventory: **8,198 states and 29,822 arcs**. This was a state-construction
check, not an exhaustive g3 minimum search.

The same catalog generator produces 1, 2, 4, 7 and 12 one-prong disc strata
for n=3 through 7, and 45, 97, 195, 684, 2,087 and 5,763 for
n=10, 12, 14, 18, 22 and 26. It produces 27 formal genus-three H/Q prong types.
These are inventory counts, not completed searches.

The D3 and g2 bounded runs finish their stated automata. The g3 witness fixture
is deliberately a subgraph with a short depth cap, and its completion flag for
the mathematical bound remains false. None of these reports claims a global
minimum solely from an automaton run.

## Upstream defects corrected

1. The N=3 seed builder allocated two singularity entries but wrote a third
   sentinel. AddressSanitizer identified the heap-buffer overflow; the buffer
   now includes that sentinel. A dedicated boundary regression was added.
2. The native result-map insertion decremented `lower_bound` even when it was
   `begin()`. The insertion now uses the valid lower-bound hint directly.
3. In 0.1.0, a generated CSparse compatibility overlay deduced its integer pointer
   type and used the C library's matching deallocators, leaving the submodule
   unchanged. Version 0.1.1 embeds these same fixes directly in the local header.

## Scope of this validation

The old distributed D10/D12/D14/D18/D22/D26 computations and all g3 sectors
were not rerun. Their historical proof status is not altered. Larger-disc
configurations are budgeted previews; normal-form and sector-coverage arguments
remain separate mathematical obligations.
