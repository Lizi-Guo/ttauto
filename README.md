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

## Configurations and computation time

The configurations below follow our recorded proof routes. `R(id)` denotes
`current_query_configuration` in the [result register](docs/HISTORICAL_RESULTS.json),
which contains the exact parameters and pruning settings for all 64 entries.

| Results | Configuration / route |
| --- | --- |
| D3 minimum | [`disc3-minimum.json`](experiments/proofs/disc3-minimum.json) |
| One-pronged disc strata, n = 4-7 | `R(Dn-sj)`; supplement the spectral exclusions with native folding searches for residuals |
| Five genus-two minima | `R(G2-H-...)` / `R(G2-Q-...)`; quotient-disc folding searches and branched-cover certificates |
| Common lower bound on 21 genus-three strata | [`genus3-bulk21.json`](experiments/proofs/genus3-bulk21.json) |
| Genus-three Q(4,2,1,1) | `R(G3-Q-4-2-1-1)`; joint invariant/anti-invariant homology |
| Other genus-three strata | Corresponding `R(G3-...)`, with the recorded sector or folding arguments and attaining map |
| Global genus-three reduction | [`genus3-all.json`](experiments/proofs/genus3-all.json) starts the spectral stage; the principal-stratum exclusions are additional |
| D10, D12, D14, D18, D22, D26 | `experiments/proofs/discN.json`, with N equal to the puncture count |

These are configuration mappings; the complete proofs also use the coverage
arguments and attaining certificates recorded in [docs/RESULTS.md](docs/RESULTS.md).

Our recorded computation times were approximately:

| Computation | Time |
| --- | --- |
| D3 lower-bound and small folding checks | 1 second each |
| Genus-three common lower bound, 21 strata | 4-6 seconds |
| Genus-three H(3,1), native search on local storage | 25-66 seconds |
| D12 coefficient enumeration | 10 minutes |
| D14 coefficient enumeration | 14 minutes |
| D18 coefficient enumeration | 3 hours |
| D22 coefficient enumeration | 10 hours |
| D26 coefficient enumeration | 5 days |

The large-disc times are elapsed windows of our distributed enumeration runs,
including resumptions and excluding later proof audits. The other rows time
only the indicated stages. See the [timing records](docs/RUNTIME_ESTIMATES.json).

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
