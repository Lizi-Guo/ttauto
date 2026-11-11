# Included algorithm sources

The implementation is self-contained at the source level. Building or running
it does not contact `jeanluct/ttauto` or `jeanluct/jlt`.
The URLs and original revisions in `dependencies.lock.json` document provenance.

| Algorithm or component | Included source |
|---|---|
| Train-track states, multigons, numbering, elementary folds and signed maps | `include/traintracks/`, `lib/traintracks/` |
| Folding automata, transition matrices and original path search | `include/ttauto/` |
| Local C++ automaton entrypoint for the unified pipeline | `native/export.cpp` |
| General ribbon state enumeration and signed-map checks | `python/ttauto_pipeline/ribbon_states.py`, `ribbon.py`, `carrying.py`, `topology.py` |
| Shared search and optional pruning modules | `python/ttauto_pipeline/automaton.py`, `engine.py`, `pruning.py` |
| Exact arithmetic and certificate primitives | `python/ttauto_cert/` |
| Matrix, polynomial and numerical utilities | `extern/jlt/jlt/` |
| Sparse matrix primitives | `extern/jlt/extern/CSparse/` |
| Independent D3 and closed-surface geometric certificates | `python/ttauto_pipeline/geometry.py`, `gates.py` |
| Automaton-free proof inventory and resumable coefficient sources | `proof_inventory.py`, `spectral_source.py`, `native_spectral.py`, `native/trace_prefix.cpp` |
| Migrated spectral, deck, base-coupling and abstract cycle kernels | `python/ttauto_cert/spectral.py`, `cover_orbits.py`, `cycle_obstruction.py` |
| Proof modules, claims, evidence replay and residual continuation | `proof_pruning.py`, `proof_backend.py`, `proof_verify.py`, `residual.py` |

The train-track and automaton core was copied from ttauto revision
`ca1b2ed030946cd1ad5dd3abe05a85a52f59b0e7`. The N=3 builder sentinel allocation
and the native search's result-map insertion were corrected locally. The
unified pipeline builds this copy directly into its executable. It does not
link to a separately installed ttauto library.

jlt was copied from revision `ea4e9558f46dd12487be7e3e85d6a3a4eaed8895`.
Its CSparse wrapper has two local fixes: pointer types follow CSparse's index
type on Windows, and C-allocated structures use their matching C deallocators.
These fixes are in the checked-in header, not injected by a build-time patch.
The matrix compound-addition and subtraction loops also declare mutable and
const iterators separately, fixing a template compilation error in jlt's tests.
No automatic update or download is performed.

The bundled CSparse keeps `Source/`, `Include/`, `Demo/`, `Config/`, its build
files and `Doc/` (the licence texts and change log). Its `MATLAB/` interface,
`Matrix/` sample inputs and `Tcov/` coverage harness are not included: nothing
here builds or reads them. CSparse's own `README.txt` still describes the full
upstream layout; take those three directories from upstream CSparse if needed.

Version 0.1.2 removes Hall Trains, its native adapter, exact-decision process
client and graph-format exporter. There is no replacement external classifier.
The built-in certificates have a stated scope; unsupported disc maps remain
unresolved. Some small historical JSON graphs remain as structural test inputs.
Their reported mapping-class types are not used to certify any run.

`source.kind = "native"` means the local C++ fold graph. `"ribbon"` means the
local Python ribbon builder. Both feed the same signed-path enumeration,
pruning, geometric classification and result journal. `"automaton"` loads a
saved graph. The earlier label `"upstream"` remains a configuration alias for
`"native"`; it has no network behavior.

Proof configurations use `newton`, `newton_native`, `polynomial_catalog` or
`residual` sources and do not construct an automaton. Their mathematical kernels
were ported from the project's large-disc and genus-three proof jobs; original
source hashes are recorded in `PROOF_MIGRATION_SOURCES.json`. Historical server
orchestration and external-classifier integrations are not imported.

Original authorship notices and license files are retained. The main project
uses GPL-3.0-or-later; jlt retains MIT, and CSparse retains its own
notices and license. `docs/UPSTREAM_README.md` is an archived
upstream document; use the root README for current build instructions.
