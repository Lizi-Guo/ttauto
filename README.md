# ttauto

**Exact computation of minimal stretch factors of pseudo-Anosov maps.**

This version of ttauto is developed jointly by **Lizi Guo, Jean-Luc Thiffeault,
and Erwan Lanneau**.

A pseudo-Anosov homeomorphism $f$ of a surface $S$ stretches its invariant
foliations by a factor $\lambda(f)>1$, its *stretch factor* (or dilatation).
On a fixed surface the set of stretch factors is discrete, so there is a least
one,

$$\delta(S)=\min\lbrace\lambda(f) : f \text{ pseudo-Anosov on } S\rbrace,$$

and the same is true inside each stratum, that is, for a prescribed list of
singularities of the invariant foliation. Determining these minima is a finite
but large computation. This repository contains the programs that carry it
out with exact arithmetic, and the results obtained with them for punctured
discs and for closed surfaces of genus two and three.

It combines C++ train-track automata with an automaton-free method for lower
bounds and certificates for the maps that attain them.

## Results

Throughout, mapping classes preserve orientation, and every $\lambda$ below
is an algebraic integer given exactly by a polynomial. Decimals are rounded
to six places; [the register](docs/HISTORICAL_RESULTS.json) has the exact
polynomial and a rational isolating interval for each entry, and
[docs/RESULTS.md](docs/RESULTS.md) gives the proof route, counts and sources.

### Punctured discs

$D_n$ is the disc with $n$ punctures; its mapping classes are braids on $n$
strands. Here the punctures may have any number of prongs.

**Theorem 1.** The minimal stretch factor of a pseudo-Anosov braid on $n$
strands, for $n = 10, 12, 14, 18, 22, 26$, is as follows.

| Surface | Minimum $\lambda$ | $\lambda$ is the largest real root of | Attained by the braid |
| --- | --- | --- | --- |
| $D_{10}$ | 1.343720 | $x^{9}-2x^{5}-2x^{4}+1$ | $(\sigma_1^2\sigma_2\cdots\sigma_9)^2\sigma_9^{-2}$ |
| $D_{12}$ | 1.251413 | $x^{12}-2x^{7}-2x^{5}+1$ | $\rho_{12}^{5}\sigma_{11}^{-2}$ |
| $D_{14}$ | 1.222574 | $x^{14}-2x^{9}-2x^{5}+1$ | $\rho_{14}^{3}\sigma_{13}^{-2}$ |
| $D_{18}$ | 1.164318 | $x^{18}-2x^{11}-2x^{7}+1$ | $\rho_{18}^{5}\sigma_{17}^{-2}$ |
| $D_{22}$ | 1.130709 | $x^{22}-2x^{13}-2x^{9}+1$ | $\rho_{22}^{5}\sigma_{21}^{-2}$ |
| $D_{26}$ | 1.108686 | $x^{26}-2x^{15}-2x^{11}+1$ | $\rho_{26}^{7}\sigma_{25}^{-2}$ |

Here $\rho_n=\sigma_1\sigma_2\cdots\sigma_{n-1}$.

**Corollary.** The minimum stretch factor $\delta(D_n)$ is now known for every
$n\ge3$. The six cases above complete the determination when combined with
the previously known small cases and the general theorem of
[Tsang and Zeng](https://arxiv.org/html/2412.01648v2).

### Strata of punctured discs, $3\le n\le 7$

Now every puncture is one-pronged, the foliation has the listed interior
singularities (by number of prongs), and the boundary has the listed number
of prongs.

**Theorem 2.** The minimal stretch factor in each of the 26 strata of
$D_3,\dots,D_7$ is as follows.

| $n$ | Interior singularities | Boundary prongs | Minimum $\lambda$ | $\lambda$ is the largest real root of |
| --- | --- | --- | --- | --- |
| 3 | none | 1 | 2.618034 | $x^{2}-3x+1$ |
| 4 | none | 2 | 2.618034 | $x^{2}-3x+1$ |
| 4 | 3 | 1 | 2.296630 | $x^{4}-2x^{3}-2x+1$ |
| 5 | none | 3 | 1.722084 | $x^{4}-x^{3}-x^{2}-x+1$ |
| 5 | 3 | 2 | 1.722084 | $x^{4}-x^{3}-x^{2}-x+1$ |
| 5 | 4 | 1 | 2.153721 | $x^{4}-3x^{3}+3x^{2}-3x+1$ |
| 5 | 3, 3 | 1 | 2.015357 | $x^{6}-x^{5}-4x^{3}-x+1$ |
| 6 | none | 4 | 1.883204 | $x^{4}-2x^{3}+x^{2}-2x+1$ |
| 6 | 3 | 3 | 1.839287 | $x^{3}-x^{2}-x-1$ |
| 6 | 4 | 2 | 1.883204 | $x^{4}-2x^{3}+x^{2}-2x+1$ |
| 6 | 5 | 1 | 2.081019 | $x^{4}-x^{3}-2x^{2}-x+1$ |
| 6 | 3, 3 | 2 | 2.081019 | $x^{4}-x^{3}-2x^{2}-x+1$ |
| 6 | 3, 4 | 1 | 1.883204 | $x^{4}-2x^{3}+x^{2}-2x+1$ |
| 6 | 3, 3, 3 | 1 | 2.171127 | $x^{8}-2x^{7}+x^{6}-4x^{5}+4x^{4}-4x^{3}+x^{2}-2x+1$ |
| 7 | none | 5 | 1.556030 | $x^{6}-x^{5}-x^{4}+x^{3}-x^{2}-x+1$ |
| 7 | 3 | 4 | 1.465571 | $x^{7}-2x^{4}-2x^{3}+1$ |
| 7 | 4 | 3 | 1.465571 | $x^{7}-2x^{4}-2x^{3}+1$ |
| 7 | 5 | 2 | 1.556030 | $x^{7}-2x^{5}-2x^{2}+1$ |
| 7 | 6 | 1 | 2.042491 | $x^{7}-2x^{6}-2x+1$ |
| 7 | 3, 3 | 3 | 1.610942 | $x^{8}-x^{7}-2x^{5}+2x^{4}-2x^{3}-x+1$ |
| 7 | 3, 4 | 2 | 2.475412 | $x^{6}-3x^{5}+3x^{4}-5x^{3}+3x^{2}-3x+1$ |
| 7 | 3, 5 | 1 | 1.809789 | $x^{8}-x^{7}-2x^{5}-2x^{3}-x+1$ |
| 7 | 4, 4 | 1 | 1.754878 | $x^{8}-x^{7}-4x^{4}-x+1$ |
| 7 | 3, 3, 3 | 2 | 1.610942 | $x^{9}-x^{7}-2x^{6}-2x^{3}-x^{2}+1$ |
| 7 | 3, 3, 4 | 1 | 2.042491 | $x^{9}-2x^{8}+x^{7}-2x^{6}-2x^{3}+x^{2}-2x+1$ |
| 7 | 3, 3, 3, 3 | 1 | 2.025975 | $x^{8}-x^{6}-3x^{5}-5x^{4}-3x^{3}-x^{2}+1$ |

### Closed surfaces of genus two

$H(k_1,\dots,k_m)$ denotes the pseudo-Anosov maps whose invariant foliations
are orientable with singularities of $2k_i+2$ prongs, and $Q(k_1,\dots,k_m)$
those whose foliations are non-orientable with singularities of $k_i+2$
prongs. A stratum means the whole stratum, all components included.

**Theorem 3.** The minimal stretch factor in each stratum of genus two is as
follows.

| Stratum | Minimum $\lambda$ | $\lambda$ is the largest real root of |
| --- | --- | --- |
| $H(2)$ | 1.722084 | $x^{4}-x^{3}-x^{2}-x+1$ |
| $H(1,1)$ | 1.883204 | $x^{4}-2x^{3}+x^{2}-2x+1$ |
| $Q(2,2)$ | 2.618034 | $x^{2}-3x+1$ |
| $Q(2,1,1)$ | 1.722084 | $x^{4}-x^{3}-x^{2}-x+1$ |
| $Q(1,1,1,1)$ | 1.839287 | $x^{3}-x^{2}-x-1$ |

### Closed surfaces of genus three

**Theorem 4.** The minimal stretch factor in each of the following 19 strata
of genus three is as follows.

| Stratum | Minimum $\lambda$ | $\lambda$ is the largest real root of | Lower-bound route |
| --- | --- | --- | --- |
| $H(1,1,1,1)$ | 1.401268 | $x^{6}-x^{4}-x^{3}-x^{2}+1$ | L |
| $Q(5,1,1,1)$ | 1.401268 | $x^{6}-x^{4}-x^{3}-x^{2}+1$ | A |
| $H(2,1,1)$ | 1.465571 | $x^{3}-x^{2}-1$ | A |
| $Q(6,1,1)$ | 1.465571 | $x^{3}-x^{2}-1$ | A |
| $H(4)$ | 1.506136 | $x^{6}-x^{5}-x^{3}-x+1$ | M + S |
| $H(3,1)$ | 1.506136 | $x^{6}-x^{5}-x^{3}-x+1$ | A |
| $Q(3,3,2)$ | 1.506136 | $x^{6}-x^{5}-x^{3}-x+1$ | A |
| $Q(5,3)$ | 1.582347 | $x^{6}-x^{4}-2x^{3}-x^{2}+1$ | A |
| $H(2,2)$ | 1.635573 | $x^{6}-2x^{5}+2x^{4}-3x^{3}+2x^{2}-2x+1$ | M + S |
| $Q(4,3,1)$ | 1.635573 | $x^{6}-2x^{5}+2x^{4}-3x^{3}+2x^{2}-2x+1$ | A |
| $Q(4,2,1,1)$ | 1.635573 | $x^{6}-2x^{5}+2x^{4}-3x^{3}+2x^{2}-2x+1$ | S (joint homology) |
| $Q(3,2,2,1)$ | 1.635573 | $x^{6}-2x^{5}+2x^{4}-3x^{3}+2x^{2}-2x+1$ | A |
| $Q(4,4)$ | 1.722084 | $x^{4}-x^{3}-x^{2}-x+1$ | S (rank four) |
| $Q(4,2,2)$ | 1.722084 | $x^{4}-x^{3}-x^{2}-x+1$ | S + independent A |
| $Q(6,2)$ | 1.883204 | $x^{4}-2x^{3}+x^{2}-2x+1$ | A |
| $Q(5,2,1)$ | 1.883204 | $x^{4}-2x^{3}+x^{2}-2x+1$ | P |
| $Q(2,2,2,2)$ | 1.883204 | $x^{4}-2x^{3}+x^{2}-2x+1$ | M (two sectors + low matrices) |
| $Q(7,1)$ | 1.987793 | $x^{6}-2x^{4}-3x^{3}-2x^{2}+1$ | S + A (spectral gap) |
| $Q(8)$ | 2.081019 | $x^{4}-x^{3}-2x^{2}-x+1$ | M + P (fixed / period two) |

The last column records how the lower bound was obtained: **A** complete
search of the stratum's automaton; **P** complete search of a projected
automaton; **S** spectral or homological exclusion; **M** mixed
periodic-point or sector arguments; **L** a published lower bound. $H(4)$ and
$H(2,2)$ are whole-stratum minima; their hyperelliptic components have the
larger minima 1.556030 and 1.839287.

For the remaining eight strata the exact minimum is open. A lower bound is
proved and an upper bound is attained:

| Stratum | Lower bound | Attained upper bound |
| --- | --- | --- |
| $Q(4,1,1,1,1)$ | 1.401268 | 1.610942 |
| $Q(3,3,1,1)$ | 1.401268 | 1.809789 |
| $Q(3,2,1,1,1)$ | 1.401268 | 2.618034 |
| $Q(3,1,1,1,1,1)$ | 1.401268 | 2.618034 |
| $Q(2,2,2,1,1)$ | 1.401268 | 1.722084 |
| $Q(2,2,1,1,1,1)$ | 1.401268 | 2.042491 |
| $Q(2,1,1,1,1,1,1)$ | 1.401268 | 1.610942 |
| $Q(1,1,1,1,1,1,1,1)$ | 1.363878 | 1.839287 |

**Theorem 5.** The minimal stretch factor $\delta_3$ of a closed surface of
genus three satisfies

$$\delta_3\in\lbrace\beta,c\rbrace,\qquad \beta\approx1.363878,\quad c\approx1.401268,$$

where $\beta$ is the largest real root of $x^{10}-x^{8}-2x^{5}-x^{2}+1$ and $c$
that of $x^{6}-x^{4}-x^{3}-x^{2}+1$. Twenty-six strata have lower bound $c$, and
$c$ is attained. Only $Q(1,1,1,1,1,1,1,1)$ can still realize $\beta$;
whether it does is open.

## Method

**Lower bounds.** Let $f$ be pseudo-Anosov with stretch factor $\lambda$ in
a given stratum. On the orientation double cover of its invariant foliation,
the lift of $f$ that preserves the orientation of the foliation acts on
homology (on the anti-invariant part, when the foliation downstairs is not
orientable) by an integer matrix $A$ of known rank $2h$. Its characteristic
polynomial is monic and reciprocal, and $\lambda$ is a simple root of it
that strictly dominates all others. So for a bound $U$ there are only
finitely many candidates: Newton's identities and
$|\operatorname{tr}A^k|\le h\,(U^k+U^{-k})$ confine the coefficients to
finite ranges, which are enumerated completely. Each candidate is then tested
against exact necessary conditions:

| Module | Condition |
| --- | --- |
| `trace`, `schur` | every root lies in the disc of radius $U$ |
| `perron` | the largest real root is simple and strictly dominant |
| `mod2` | on a disc, the polynomial mod 2 is that of a permutation of the branch points |
| `peripheral` | on a disc, compatibility with the transition matrix of a train track |
| `lefschetz` | the Lefschetz numbers of all iterates give non-negative integer counts of periodic orbits |
| `deck`, `base_coupling` | the same for the second lift, coupled through the deck involution |
| `cycle` | on a disc, a primitive matrix with this polynomial has a connected cycle graph |

If no candidate with largest root below the target survives, the target is a
lower bound for the stratum. A candidate that survives needs a further
argument; the general one is an exhaustive search of the stratum's
train-track automaton, which the Ham–Song norm bound makes finite.

**Upper bounds.** A minimum is attained by an explicit map, given as a closed
path of folds in a train-track automaton. That the map is pseudo-Anosov in
the stated stratum is certified from the path itself with the
Bestvina–Handel criterion.

**Exactness.** Every decision uses integers and rationals only; real roots
are compared with Sturm sequences. A resource limit never counts as an
exclusion: a candidate that exhausts its budget stays unresolved. Each proof
depends on the topological theorems it records, and replaying its arithmetic
is not a formalization in a proof assistant. The contracts of the modules are
in [docs/PROOF_PIPELINE.md](docs/PROOF_PIPELINE.md) and
[docs/PIPELINE.md](docs/PIPELINE.md), and the written arguments in
[docs/proofs](docs/proofs).

## Installation

Python 3.10 or later, with no packages beyond the standard library. The
native programs need CMake 3.16 or later and a C++17 compiler (GCC or Clang).
All sources are in the repository; nothing is downloaded.

```sh
python -m pip install -e .
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Without installing, set `PYTHONPATH=python` and use `python -m ttauto_pipeline`
in place of `ttauto`.

## Usage

A computation is described by a JSON file: the surface, the strata, a rational
bound `upper`, and optionally the polynomial whose largest root is the claimed
minimum.

```json
{
  "schema": "ttauto-proof-v1",
  "surface": {"kind": "closed", "genus": 3},
  "strata": [{"type": "Q", "orders": [5, 1, 1, 1]}],
  "upper": "701/500",
  "target_polynomial": [1, 0, -1, -1, -1, 0, 1]
}
```

```sh
ttauto run experiments/proofs/genus3-bulk21.json --output runs/g3   # run a proof
ttauto verify runs/g3 --replay-search                               # replay it
ttauto run experiments/proofs/disc12.json --output runs/d12 --events 1000
ttauto run experiments/proofs/disc12.json --output runs/d12 --resume --events 100000
ttauto modules                                                      # list the tests
```

A run writes `result.json`, which states `lower_bound_proved`, and
`minimum_proved` when a certified attaining map is supplied as well. Runs are
transactional and resume exactly where they stopped. The second kind of
computation, `ttauto-experiment-v1`, searches the closed paths of a
train-track automaton; `ttauto catalog` and `ttauto campaign` generate and
run one configuration per stratum. [docs/USAGE.md](docs/USAGE.md) is the
full reference.

## Configurations and estimated runtime

The choices below are inferred from the recorded proof routes. The
[result register](docs/HISTORICAL_RESULTS.json) supplies the exact target and
surface data for every row; [the timing register](docs/RUNTIME_ESTIMATES.json)
maps all 64 results to these choices and records the timing evidence and
assumptions. No local search was rerun to make these estimates. A configuration
for a lower-bound question does not by itself supply the historical coverage
argument or an attaining map.

Use these module profiles, with unlisted optional modules disabled:

| Profile | Schema and source | Enabled modules |
| --- | --- | --- |
| **D**: disc spectra | `ttauto-proof-v1`; `newton_native` for large discs, `newton` for small arithmetic checks | `trace`, `schur`, `perron`, `mod2`, `peripheral`, `lefschetz`, `deck`, `cycle` |
| **H**: orientable closed stratum | `ttauto-proof-v1`; `newton` | `trace`, `schur`, `perron`, `lefschetz` |
| **Q**: nonorientable foliation on a closed surface | `ttauto-proof-v1`; `newton` | H modules plus `base_coupling` |
| **F**: residual folding search | `ttauto-experiment-v1`; `native` for standard disc automata, or `automaton` for a supplied certified graph | `norm`, `cycle`, `homology`, `lefschetz`, `perron`; also `mod2` for discs; `geometry.kind: builtin` |

For H/Q use trace horizon 32. Large-disc examples use native source horizon 48,
pruning horizon 192 and `batch_nodes: 10000`. The saved common genus-three
configuration already selects the applicable modules; inapplicable modules skip.
In a spectral configuration, `target_polynomial` defines a strict lower-bound
question. In F it selects equality with a root: omit it when searching for
all smaller maps, and omit a manual `limits.depth` for an exhaustive bounded
search. A fixed regular mark alone does not cover a whole closed stratum.

### Which configuration produces which part of a result?

`R(id)` below means the `current_query_configuration` of that ID in the result
register. H/Q labels and disc indices are those in [docs/RESULTS.md](docs/RESULTS.md).
The timing column states the phase being estimated, not a guarantee that the
current programs will finish the entire historical proof.

| Result | Configuration inferred from the records | Time and scope |
| --- | --- | --- |
| $D_3$ minimum | [`disc3-minimum.json`](experiments/proofs/disc3-minimum.json), including its `fold_path` witness | Plan **1–5 s** on one process. Saved CI lower-bound and small folding steps each took **0–1 s** at one-second resolution; the complete combined configuration has no separate saved timer. |
| The other 25 one-pronged disc strata, $4\le n\le7$ | `R(Dn-sj)`, D; use F for the historical residuals, followed by the recorded attaining map | Per-stratum totals vary. Native estimates for the difficult recorded residuals are given below; there is no measured total for all 26 strata. |
| All five genus-two minima | `R(G2-H-...)` or `R(G2-Q-...)` gives the arithmetic question. The recorded complete route uses F on the hyperelliptic quotient disc cases, then the branched-cover certificate | **30 s–25 min** for the six counted native quotient searches together; additional graph, cover and witness checks are not included. [`closed-g2.json`](experiments/closed-g2.json) covers only its regular-fixed-point sector. |
| Common bound $c$ on 21 genus-three strata; in particular the lower half of the $H(1,1,1,1)$ and $Q(5,1,1,1)$ minima | [`genus3-bulk21.json`](experiments/proofs/genus3-bulk21.json); supply each attaining map separately | **4–6 s observed**, then **9–11 s observed** for full replay in saved CI. No large automaton is needed for this lower bound. |
| Genus-three $Q(4,2,1,1)$ minimum | `R(G3-Q-4-2-1-1)`, Q, then its known attaining map; retain the joint invariant/anti-invariant argument | **1–30 min planning estimate** for the arithmetic stage, with low confidence: 1,107 historical coefficient triples and 54 joint sequences. This is not a measured current run. |
| Genus-three $H(4),H(2,2),Q(4,4),Q(4,2,2)$ | Corresponding H/Q query plus the recorded periodic-point, rank-four or sector argument and attaining map | Arithmetic is the intended first step. The available summaries do not give comparable complete timings for these mixed proofs; **no defensible whole-proof ETA** is assigned. |
| The other 12 proved genus-three strata | Corresponding H/Q query, then F on the required complete/projected graph or sectors if arithmetic leaves residuals; separate attainment | Native workload estimates for $H(3,1)$ and $H(2,1,1)$ appear below. The large projected and two-sector cases need their original inputs; **no single runtime estimate applies to this group**. For $Q(8)$, [`closed-g3-witness.json`](experiments/closed-g3-witness.json) checks only the upper witness. |
| Global genus-three two-value reduction | [`genus3-all.json`](experiments/proofs/genus3-all.json) starts the spectral question; the recorded principal-stratum and P12 exclusions are additional steps | The **4–6 s** bulk computation covers only 21 strata. The complete two-value reduction has no measured current runtime; the remaining beta realization question is open. |
| Eight open genus-three stratum intervals | Corresponding `R(G3-Q-...)`, Q, targeting the recorded **lower** polynomial; certify the recorded upper map separately | There is **no completion-time estimate for an exact minimum**. Neither an event budget nor an empty partial residual list resolves an open case. |
| $D_{10},D_{12},D_{14},D_{18},D_{22},D_{26}$ | D with [`disc10.json`](experiments/proofs/disc10.json), [`disc12.json`](experiments/proofs/disc12.json), [`disc14.json`](experiments/proofs/disc14.json), [`disc18.json`](experiments/proofs/disc18.json), [`disc22.json`](experiments/proofs/disc22.json), [`disc26.json`](experiments/proofs/disc26.json) | Native enumeration scenarios below; orbit/cycle exclusions and attainment require additional time and evidence. |
| Other puncture numbers in the all-$n$ corollary | Apply the cited general theorem and known small cases | No new enumeration is required. |

The genus-two quotient cases are `(5,0)` for H(2), `(6,0)` for H(1,1),
`(3,0)` and `(4,0)` for Q(2,2), `(5,0)` with the regular branch-point condition
for Q(2,1,1), and `(5,3)`, `(6,1)`, `(6,4)` for Q(1,1,1,1).
For a small disc row `(n,I,O,j)`, adapt
[`disc3-native.json`](experiments/disc3-native.json): set `punctures=n`,
`inner_prongs=[1]*n+I`, `inner_colors=["puncture"]*n+["filled"]*len(I)`,
`outer_prongs=[O]`, `source.stratum_index=j`, and the row's rational `upper`.
For a closed F search, the historical complete/projected graph must be supplied;
a generic auxiliary-mark catalog is not a substitute for its coverage proof.

To extract any recorded spectral query without a search:

```python
import json
from pathlib import Path

records = json.loads(Path("docs/HISTORICAL_RESULTS.json").read_text())["records"]
row = next(r for r in records if r["id"] == "G3-Q-4-2-1-1")
Path("query.json").write_text(json.dumps(row["current_query_configuration"], indent=2))
```

```sh
ttauto validate query.json
ttauto run query.json --output runs/query --events 100000
ttauto run query.json --output runs/query --resume --events 1000000
```

Event limits are total budgets, not time limits. The register's 1,000-event
defaults and the large-disc examples are starting budgets. Increase them on
resume, or omit `limits.events` for an uncapped source enumeration. Exhausted
`models`, `cycle_nodes` or `base_catalog_nodes` limits leave unresolved cases;
use a new residual run to change those limits. Whole minima additionally need
matching, in-scope geometric witnesses. Some historical coverage arguments and
stronger cycle obstructions still need separate integration, so these inferred
recipes are not claims that every result is already obtainable with one command.

### Native search estimates

The following are **enumeration-only scenarios**, excluding graph construction,
Python/SQLite overhead, orbit tests, geometric certification and replay. The
historical algorithms and current module choices can visit different trees.
For folding, assume **10,000–500,000 prefixes/s per worker**: the saved H(3,1)
search processed 12,474,712 prefixes in 25–66 s on local storage, but took 937 s
in a paired shared-storage control. The wider assumed range allows for different
matrix dimensions and per-prefix work; it is not a runtime bound.

| Historical folding workload | Prefixes | Estimated native time, one worker |
| --- | ---: | --- |
| D6-s3 residual | 20,562 | 0.04–2.1 s |
| D6-s6 residual groups | 15,700,700 | 31 s–26 min |
| D7-s6 last three stages | 195,682,145 | 6.5 min–5.4 h |
| D7-s11 final phase only | 2,809,458 | 6 s–4.7 min |
| Six genus-two quotient searches | 14,853,946 | 30 s–25 min |
| Genus-three H(2,1,1) search | 190,150,702 | 6.3 min–5.3 h |

These phase counts are not complete timings for the surrounding proofs. As a
separate example of **rediscovering** a known attaining map, the archived
Q(5,1,1,1) automaton had 5,135,719 unexcluded start states. Its observed sample
rate of 0.0762 s/start suggests about **109 worker-hours**, or **1–5 days with
two workers** allowing a factor-two variation, before graph construction and
audits. The historical memory admission was 128 GiB. Using the known witness
and the bulk lower bound avoids that rediscovery search.

For large-disc coefficients, the initial distributed waves averaged about
0.52–0.85 million nodes/s per worker. The scenario below assumes
**0.2–1 million nodes/s per worker**, applied to the final historical node counts.
The current engine uses one worker per run. The 64-worker column is an **ideal
externally partitioned workload**, not a supported `workers` configuration flag
or a measured speedup; a valid partition and coverage accounting are required.

| Disc | Historical native nodes | One-worker scenario | Ideal 64-worker scenario |
| --- | ---: | --- | --- |
| D10 | Not recorded in the consulted summary | 12–60 h (low-confidence proxy) | 11–56 min (same proxy) |
| D12 | 18,890,600,014 | 5.2–26 h | 5–25 min |
| D14 | 31,881,873,712 | 8.9–44 h | 8–42 min |
| D18 | 1,119,625,304,756 | 13–65 days | 4.9–24 h |
| D22 | 14,387,847,199,303 | 167–833 days | 2.6–13 days |
| D26 | 151,539,175,623,868 | 4.8–24 years | 27–137 days |

The D10 proxy scales the D12 estimate by retained spectra, `358430/156012`;
its 12,431,171,156 recorded coefficient completions are a different counter.
These estimates are sensitive to pruning, hardware and batching. The historical
distributed enumeration windows were about **10 min, 14 min, 3 h, 9.8 h and
5.2 days** for D12 through D26, using changing allocations and continuations;
they are not measurements on a constant 64-worker machine. Also, those campaigns
used a theorem to remove smaller one-pronged cores, whereas the current disc
examples enumerate all `m=3,...,n`. The scenario is a workload scale, not a
prediction that the current serial implementation completes within that time.

Finally, the saved final auditor invocations took **19.5 min, 4.1 min, 1.9 min,
3.5 h and 11.2 h**, respectively. They reused existing evidence, including prior
replay checkpoints in D18/D22/D26, and must not be reported as fresh search times.
The independent saved CI runs used for the small examples are
[run 42953696873](https://github.com/Lizi-Guo/ttauto/actions/runs/42953696873) and
[run 42953012802](https://github.com/Lizi-Guo/ttauto/actions/runs/42953012802).

## Repository layout

| Path | Contents |
| --- | --- |
| `include/`, `lib/` | train tracks, folds and folding automata (C++) |
| `examples/`, `tests/`, `testsuite/` | programs and tests of the C++ core |
| `python/ttauto_pipeline/` | the engine, its command line, and the two back ends |
| `python/ttauto_cert/` | exact arithmetic and the necessary conditions |
| `native/` | coefficient enumeration and automaton export (C++) |
| `experiments/` | configurations, including `experiments/proofs/` |
| `docs/` | contracts, written arguments, results register, validation records |
| `doc/` | the mathematics and structure of the C++ core |
| `extern/jlt/` | bundled numerical library and CSparse |

[docs/SOURCES.md](docs/SOURCES.md) maps each algorithm to its file.

## References

1. M. Bestvina and M. Handel, *Train-tracks for surface homeomorphisms*,
   Topology **34** (1995), 109–140.
2. J.-Y. Ham and W. T. Song, *The minimum dilatation of pseudo-Anosov
   5-braids*, Experiment. Math. **16** (2007), 167–179.
3. E. Lanneau and J.-L. Thiffeault, *On the minimum dilatation of braids on
   punctured discs*, Geom. Dedicata **152** (2011), 165–182.
4. E. Lanneau and J.-L. Thiffeault, *On the minimum dilatation of
   pseudo-Anosov homeomorphisms on surfaces of small genus*, Ann. Inst.
   Fourier **61** (2011), 105–144.
5. E. Hironaka and C. C. Tsang, Groups Geom. Dyn. **19** (2025),
   Propositions 3.12 and 3.21, as cited in
   [docs/proofs](docs/proofs/disc_weight_homology.txt).
6. B. Farb and D. Margalit, *A Primer on Mapping Class Groups*, Princeton
   University Press, 2012.
7. C. C. Tsang and X. Zeng, *Minimum dilatations of pseudo-Anosov braids*,
   [arXiv:2412.01648v2](https://arxiv.org/abs/2412.01648v2), 2025,
   Corollary 1.2 and Section 1.1.

## Citation

```bibtex
@misc{ttauto-pipeline,
  author       = {Guo, Lizi and Thiffeault, Jean-Luc and Lanneau, Erwan},
  title        = {ttauto: exact computation of minimal stretch factors of pseudo-Anosov maps},
  year         = {2026},
  howpublished = {\url{https://github.com/Lizi-Guo/ttauto}}
}
```

## License

GPL-3.0-or-later. The bundled jlt and CSparse components retain their respective
licenses and copyright notices.
