# Mathematical and execution contracts

This page describes the explicit **folding-path backend**. The common engine
also accepts automaton-free proof sources; see [PROOF_PIPELINE.md](PROOF_PIPELINE.md).

## Inputs and enumeration

The surface is orientable. A disc has genus zero and one distinguished outer
face. Closed experiments have genus at least two. Inner and outer data are
**prong counts**, not zero orders; a regular point has two prongs. Unmarked
filled one-prong points are rejected. Inner polygons have independent puncture
or filled colors, which canonicalization and terminal symmetries preserve.

For a prescribed type with genus g and I/O inner/outer faces, the real dimension
is m = 2g - 2 + I + O. The builder enumerates connected labelled ribbon pairings
with positive fans and the requested complementary prongs, then quotients by
oriented colored isomorphism using rooted polygon walks. It never pools all
pairings or all factorial polygon permutations. This is still expensive for
large m. A decoration budget is an operational stop, not an exclusion.

Each arc includes its complete signed word images and vertex permutation. Its
inverse, endpoint continuity, peripheral transport, and ordered noncrossing band
carrying are checked. The independent carrying checker uses chord nonalternation.
Every closed path is tested with all color-preserving oriented terminal symmetries.
The engine does not merge paths by their unsigned matrix or replace a map by a power.

## Finite bound

For an irreducible nonnegative integer matrix of size m with Perron root at most
the exact rational U, use the Ham--Song bound

    sum(M) <= floor(U^m) + m - 1 = B.

A left elementary fold is permutation-plus-one; it raises the total entry sum
of the prefix by at least one. A terminal permutation preserves the sum.
Hence at most B-m folds need be considered. This depth limit remains present
when optional norm pruning is disabled. A user depth smaller than B-m is marked
truncated even after its finite DFS finishes.

## Optional modules

| Module | Stage | Rejection means | Scope |
|---|---|---|---|
| norm | prefix | sum exceeds B | only expanding permutation-plus-one folds |
| cycle | closed matrix | exact selected-cycle growth exceeds U | never a prefix with future row permutations |
| homology | signed map | positive orienting-lift growth cannot equal the real growth | prescribed formal carrier/prongs |
| mod2 | signed map | no permitted colored branch permutation has the lift's characteristic polynomial modulo 2 | sphere double covers only |
| lefschetz | signed map | no singular/regular periodic-orbit pattern survives the chosen horizon | both lifts when the orientation cover is connected |
| perron | signed map | repeated dominant root or proven rotational companions | positive orienting lift only |

Disabling any of these does not disable final primitivity, exact growth or geometry.
Finite cycle-discovery limits and Lefschetz horizons weaken tests; a successful
necessary test never asserts realizability. Inapplicable tests return `skip`.
Errors propagate to a failed, restartable run. They never return `reject`.
The formal-carrier tests must not be applied to an unrelated polynomial, a
quotient map after an unverified collapse, or a guessed cover.

The orientation cover is built from the **full** ribbon graph. Crossing a band
flips the local tangent orientation label. Its filled cellular H1 is computed
using an exact integral quotient. A connected double cover doubles odd prongs
and splits even prongs; that identity and the Euler genus are checked. For an
orientable carrier the base homology is normalized by its tangent-action sign.
For a connected cover, the deck matrix J is computed from the actual lifted
graph. The two actual lift actions are A and JA; Lefschetz does not replace
JA with -A on a higher-genus base. The trace API separates the supplied action
from the foliation-orientation sign.

## Geometry

`geometry.kind = "builtin"` uses only project-owned Python certificate code.
`"audit"` checks the signed surface-map pipeline but never claims attainment.
There are no external classifier calls; executable and timeout geometry options
are rejected. All generated configurations use built-in certificates by default.

For D3 with three one-prong punctures and a one-prong distinguished outer mark,
the actual signed map lifts to the torus branched over these four points. The
integral H1 matrix A has determinant 1. The four-punctured-sphere correspondence
and torus classification give pseudo-Anosov type precisely when |trace(A)| > 2.
Its stretch factor is the larger root of x^2 - |trace(A)|x + 1. The code compares
that root exactly with the candidate's root and checks the requested orientability.
This use of homology is confined to the four-marked sphere; it does not classify
maps with additional punctures or higher genus. See Farb--Margalit,
[*A Primer on Mapping Class Groups*, §§2.2.5 and 13.1](https://pagine.dm.unipi.it/~a019210/Farb%20Magalit_Primer%20on%20Teichmuller%20theory.pdf).

Other disc types return `unresolved`, never a rejection or an attained claim
based solely on a primitive matrix or a finite set of necessary conditions.
The former Hall backend and its protocol/export code have been removed.

For closed surfaces, inner faces are filled to form a spine. The derivative
pair-orbit algorithm computes gates and the finite closure of taken turns. The
canonical infinitesimal graph is reconstructed for cyclic or one-edge-missing
cyclic local patterns. In particular, two gates are joined by **one edge**, not
an artificial bigon. Other patterns remain unresolved.

The reconstructed embedded carrier is checked again for signed paths, vertex
and face transport, genus, switch intertwining and reachability of every
infinitesimal band. Primitivity and growth > 1 give positive PF data; switch
intertwining and the periodic infinitesimal block give the switch equations.
Attainment uses the efficient filling carrier construction of Bestvina--Handel
(1995), sections 3.3--3.4. Actual complementary prongs and tangent orientability
are computed from the reconstruction before matching the requested stratum.
This criterion is not inferred from irreducibility of closed homology alone.

## Restart, accounting, verification

SQLite commits candidate records, counters and the explicit DFS frontier in one
transaction. A crash replays only the uncommitted batch. Computation failure rolls
back that batch. Restart identity binds mathematical configuration, source code,
native executables and the complete automaton. Bounds and module policy cannot
change during resume. Increasing the event budget is allowed.

Candidate records and checkpoints carry SHA-256 digests. `verify` checks these
bindings, independently reconstructs every path's signed map and real matrix,
and reruns each decision and geometric certificate, comparing the entire geometry
record without backend-report exceptions. `--replay-search` also
re-enumerates the bounded run and compares counters and the attained minimum
among the enumerated witnesses. It is deliberately more expensive.

An event is either a proposed fold prefix or a closed path with a terminal
symmetry. Module counts record the first-rejection funnel, not hypothetical
rejections by every module. Search completion and whole-surface minimum proofs
are distinct fields. No unverified coverage flag from an input file upgrades
the report to a global minimum theorem.

## Extending pruning

Implement `name`, `stage`, `contract`, and `check(Context) -> Decision` and add the
class to `REGISTRY`, the strict configuration keys and CLI options. A rejection
must include an exact necessary-condition argument; a limit or unsupported case
returns `skip` or an operational error. Add a complete small on/off comparison,
a boundary-equality case and a mutation control. Minima-only dominance/bad-word
pruning must not be plugged into full geometric-witness enumeration without a
separate proven objective contract.
