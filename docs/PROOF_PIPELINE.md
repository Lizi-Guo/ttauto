# Unified proof workflow (0.2.0)

`engine.run` is one transactional event loop with two families of candidate
sources. A train-track automaton is an optional representation, not a prerequisite
for a lower bound. No Hall source, process, report flag or imported classification
is used to certify a result.

## Spectral configuration

```json
{
  "schema": "ttauto-proof-v1",
  "surface": {"kind": "closed", "genus": 3},
  "strata": [{"type": "Q", "orders": [5, 1, 1, 1]}],
  "upper": "701/500",
  "target_polynomial": [1, 0, -1, -1, -1, 0, 1],
  "source": {"kind": "newton"},
  "pruning": {"trace_horizon": 32},
  "limits": {"events": 100000, "models": 500000}
}
```

Orders are differential zero orders here, unlike prong counts in explicit
automaton specifications. `strata:"all"` enumerates all formal H/Q partitions.
Every connected component and every permutation of equal-order zeros is retained.
Formal empty strata only enlarge the containing search space.

`upper` is a rational containing spectral radius. When a target polynomial is
given, its largest real root must lie strictly in `(1, upper)` and the proof
excludes dilatations *strictly below* that root. Without a target polynomial,
the rational upper itself is the proposed lower bound; no attainment is implied.

For a disc specify `{"kind":"disc","punctures":12}`. The inventory fills
non-one-pronged interior punctures and includes every one-pronged Euler stratum
for every `3 <= m <= 12`, preserving the distinguished outer mark. This is a
lower-bound superset; arbitrary extra markings are not asserted realizable.
`one_pronged:true` restricts the question to exactly n one-pronged punctures.
Explicit disc strata use `n`, `interior_prongs` (unmarked singularities >=3),
and `outer_prongs`. No published lower-bound table is silently substituted for
the omitted smaller-m searches.

## Candidate sources and completion

* `newton`: exact Python integer/Fraction recursion. For degree 2h,
  `|Tr(A^k)| <= h(U^k + U^-k)` bounds the first h coefficients through Newton
  identities. The reciprocal tail is determined. All coefficient prefixes are
  retained in a resumable DFS cursor. This is a finite containing set, including
  spectra which do not correspond to any surface map.
* `newton_native`: bundled `native/trace_prefix.cpp`, built as
  `build/bin/ttauto_trace_prefix`. Configure `executable`, `batch_nodes` (default
  10,000), and source `trace_horizon` (default 32). A batch returns spectra AND
  every deferred prefix. Optional `trace` and sphere `mod2` filters can run early.
  The native stage uses all uncolored branch permutations; later colored checks
  refine it. Early filtering changes intermediate counts, not the covered maps.
  The native implementation supports rank <=48, horizon <=192 and trace bounds
  <=10^15 with checked int64 coefficients and safe int128 intermediates. Range
  errors fail the run; choose Python integers for larger ranges. The pruning
  `trace_horizon` can be 192 while native coefficient enumeration uses 48.
* `polynomial_catalog`: JSON list, or object with `polynomials`, at `path`.
  Imported spectra are rechecked. Historical `completed`/`minimum_proved` flags
  are ignored. This source never independently establishes catalogue coverage.
* `residual`: `path` to a previous spectral run. The parent evidence is copied
  with a transactional SQLite backup, reverified, and retained in
  `parent_evidence/`. Only unresolved polynomial/stratum occurrences are emitted.
  Surface, strata and target must agree. Modules and model budgets may change.
  A complete parent containing search transfers its coverage; a partial parent
  never becomes complete just because its *current* residual list is empty.

All sources feed the same `engine.run` loop and SQLite journal. Event limits
are total committed budgets. Use `--resume --events N` to increase the event
limit. Mathematical policies and module budgets are identity-bound; use a new
residual run to strengthen exclusions or increase exhausted model/cycle budgets.

```sh
ttauto import-catalog old/polynomials.txt --output spectra.json
ttauto import-catalog old/accepted.json --field survivors --output spectra.json
ttauto run residual-config.json --output runs/stronger-residuals
ttauto verify runs/stronger-residuals --replay-search
```

## Modules and their contracts

All named modules are optional and configurable with `pruning` or the CLI's
`--enable` / `--disable`. Default spectral modules are enabled except `cycle`.
The first-half Newton bounds define the containing source, and the strict target
interval defines the mathematical question; those are not optional heuristics.

| Module | Domain | Necessary condition |
|---|---|---|
| trace | reciprocal spectrum | extended exact reciprocal-pair trace bounds |
| schur | reciprocal spectrum | exact strict Schur disk test at the containing radius |
| perron | positive expanding spectrum | no repeated dominant root, dominant negative root, or rotating dominant factor |
| mod2 | sphere orienting cover | branch permutations respect singularity colors and fixed outer point |
| peripheral | one-pronged disc model | `chi_M = chi_H chi_B`, B a signed permutation of even-prong peripheral directions; trace and formal orbit lower bounds |
| lefschetz | stratum and spectrum | nonnegative integral ordinary periodic-orbit counts for allowed singular rotations |
| deck | disc lifts | singular deck transport, regular orbit pairing and, if mod2 enabled, the same branch permutation in both tests |
| base_coupling | closed nonorientable foliation | correlated singular indices, invariant base spectrum and anti-invariant spectrum |
| cycle | disc residual real polynomial | exhaustive containing cycle disjointness graphs; a primitive matrix requires connected cycle intersection graph |

`deck` includes its necessary lift-index calculation, even when the standalone
Lefschetz stage is disabled. `base_coupling` likewise includes ordinary orbit
integrality as part of its coupling condition. Their contracts explicitly say
what they test; these dependencies never silently enable another configured stage.

The abstract cycle search enumerates all signed peripheral factors first.
Its pruning uses coefficient matching, simple-cycle counts bounded by primitive
orbit counts, and total disjoint length bounded by matrix dimension. It does
not use an unproved monotonicity of formal logarithmic coefficients. Finding
a connected abstract graph leaves an unresolved candidate, not a realization.
The `cycle` module on an explicit folding matrix is separately the existing
short-cycle growth bound; the CLI reports each module's domain.

`limits.models`, `limits.cycle_nodes`, and `limits.base_catalog_nodes` are
operational bounds. Exhaustion returns an unresolved/skip outcome, not reject.
Inapplicable modules also skip. Prefix events, spectra, stratum occurrences and
abstract graph nodes have different units and are counted separately. Native
statistics include coefficient completions retained/pruned inside each batch;
they must not be called train-track states or closed paths.

## Correct lifts and genus-three migration

For an actual connected orienting cover the code constructs the deck map J from
the lifted signed graph and computes its filled H1 action. It checks `J^2=I`,
`JA=AJ`, and the invariant rank `2g`. Lefschetz uses the actual polynomials of
`A` and `JA`. Foliation reversal changes local indices, not the supplied traces.
The older `compatible(p, negative=True)` helper is retained only for contexts
where p is the positive-normalized action and the other action really is `-A`.

For a closed genus-g nonorientable foliation with b odd zeros, invariant H1 has
rank 2g, anti-invariant H1 rank `2g-2+b`, and the positive stretch eigenvalue lies
in the anti-invariant part. The two lift actions restrict as `(A+, A-)` and
`(A+, -A-)`, not `(A+, A-)` and `(-A+, -A-)`.

The migrated closed-Q kernel keeps correlated singular vectors `(S(k), Z(k), E(k))`.
With anti-invariant traces t- and base traces t+ it computes

    R(k) = t-(k) - S(k)
    R_negative(k) = (2 - t+(k) - Z(k) + t-(k) + E(k))/2.

Ordinary exact-period counts use all divisors; negative-return counts use only
divisors d with k/d odd, and must lie between zero and the ordinary orbit count.
The base catalogue contains every reciprocal spectrum within the envelope,
including nonexpanding spectra. No base is omitted because it lacks a positive
expanding eigenvalue. Shared +/- stretch roots and exact disk obstructions can
exclude base spectra before the full correlated coupling.

The supplied bulk21 configuration ports the 2026-10-10 lower-bound argument.
Other genus-three strata remain subject to their actual residual outcomes;
this port does not assert that a global genus-three minimum has been solved.
The source mapping and hashes are in `PROOF_MIGRATION_SOURCES.json`.
The project's written arguments are retained in
[`proofs/disc_weight_homology.txt`](proofs/disc_weight_homology.txt) and
[`proofs/genus3_bulk_homology_lower.tex`](proofs/genus3_bulk_homology_lower.tex).
Their dated job names and counts describe the original proofs; this version's
new replay counts are recorded separately in `PROOF_VALIDATION_0.2.0.json`.

## Upper witnesses and report claims

A proof may optionally specify a `fold_path` witness with `path` to a saved
automaton, `start`, `folds`, and `terminal`. This loads only that supplied graph
and replays its signed maps with the built-in geometric validator. Its actual
surface/stratum and algebraic stretch root must match the claimed target.

`lower_bound_proved` requires complete containing coverage and no unresolved
occurrences. `minimum_proved` additionally requires an in-scope independently
checked attaining witness with the exact target root. Spectral survivors never
become attained maps. The existing geometry supports the D3 pillowcase and
specified closed canonical carriers. Historical large-disc upper certificates
whose classification depended on Hall are NOT promoted to Hall-free certificates.
A general replacement classifier is not included in this migration.

Every proof is conditional on its recorded topological theorem dependencies;
arithmetic replay is not a proof-assistant formalization. `verify` reconstructs
every decision, count and claim. `--replay-search` additionally repeats the full
containing enumeration and checks its frontier/checkpoint digest. Native replay
requires the content-bound native executable; the Python route needs no native
tools. Old 0.1.2 runs remain historical artifacts and cannot resume across this
mathematical code change.
