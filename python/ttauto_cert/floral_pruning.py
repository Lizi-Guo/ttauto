"""Parameterized necessary conditions from the n26 draft, with explicit scope.

No n26 minimum or n26-specific case catalogue is used. These routines require
the floral setup described in their contracts; arbitrary homology matrices or
general train-track prefixes do not supply that setup automatically.
"""
from fractions import Fraction
from collections import Counter
from .cycle_pruning import _matrix


def _integers(*values, minimum=1):
    if any(type(v) is not int or v < minimum for v in values):
        raise ValueError('integer geometric lengths/counts required')


def floral_completion_budget(lower, upper, *, zeta, filaments, entry_sum_limit):
    """Minimum entry sum within an interval box subject to floral parity.

    The caller supplies a fixed, geometrically justified canonical permutation
    zeta (column j -> row zeta[j]) and filament set. A-zeta is nonnegative;
    filament rows are even; each petal column of G has even sum, where filament
    rows of A-zeta are divided by two. A partial column with odd current sum
    is repairable if its interval has room. No current odd prefix is rejected
    merely for being odd. Output is a necessary condition, not realization.
    """
    lower, upper = _matrix(lower), _matrix(upper)
    size = len(lower)
    _integers(entry_sum_limit, minimum=0)
    if (len(upper) != size or any(type(v) is not int for v in zeta) or
            sorted(zeta) != list(range(size))):
        raise ValueError('wrong permutation or interval dimensions')
    filaments = frozenset(filaments)
    if any(type(v) is not int or not 0 <= v < size for v in filaments):
        raise ValueError('filament index out of range')
    if any((j in filaments) != (zeta[j] in filaments) for j in range(size)):
        raise ValueError('canonical permutation must preserve filament/petal types')
    if any(lower[i][j] > upper[i][j] for i in range(size) for j in range(size)):
        raise ValueError('empty matrix-entry interval')
    total, repairs = 0, []
    for j in range(size):
        units_sum, repair_costs = 0, []
        for i in range(size):
            canonical = int(zeta[j] == i)
            step = 2 if i in filaments else 1
            units = max(0, (lower[i][j]-canonical+step-1)//step)
            value = canonical+step*units
            if value > upper[i][j]:
                return {'excluded': True, 'reason': 'canonical_or_filament_parity', 'cell': [i, j]}
            total += value
            units_sum += units
            if value+step <= upper[i][j]:
                repair_costs.append(step)
        if j not in filaments and units_sum % 2:
            if not repair_costs:
                return {'excluded': True, 'reason': 'petal_exit_parity', 'column': j}
            cost = min(repair_costs)
            total += cost
            repairs.append({'column': j, 'minimum_extra_entries': cost})
    return {'excluded': total > entry_sum_limit, 'minimum_entry_sum': total,
            'entry_sum_limit': entry_sum_limit, 'parity_repairs': repairs,
            'reason': 'entry_sum_budget' if total > entry_sum_limit else 'not_excluded',
            'realizability_claimed': False}


def secondary_euler_budget(n, p, q, petals, s, h):
    """Draft (5), for TWO DISTINCT secondary orbits in the floral setup.

    p,q are primary period bounds, petals is the actual total petal count,
    s,h are embedded secondary lengths; ordinary secondary periods divide
    p,q respectively. This includes a regular first orbit with s <= p.
    Selected or summed petal curves cannot replace the actual petal count.
    """
    _integers(n, p, q, s, h)
    _integers(petals, minimum=0)
    excess = [max(0, s-3*p), max(0, h-3*q)]
    used = p+q+petals+sum(excess)
    return {'excluded': used > n, 'required_budget': used, 'available_budget': n,
            'secondary_excess': excess, 'reason': 'secondary_euler_budget' if used > n else 'not_excluded'}


def full_real_polynomial_obstruction(coefficients, *, complete_standard_real_matrix):
    """Unimodularity/reciprocity only for a FULL standard real matrix.

    The standard geometric model and based-fold factorization are caller
    preconditions. Never apply this to a prefix, relaxed clique subgraph,
    eigenvalue factor or signed homology factor.
    """
    if complete_standard_real_matrix is not True:
        raise ValueError('the full standard real transition polynomial is required')
    p = list(coefficients)
    if len(p) < 2 or p[0] != 1 or any(type(v) is not int for v in p):
        raise ValueError('a monic integer polynomial is required')
    if abs(p[-1]) != 1:
        return {'excluded': True, 'reason': 'full_matrix_not_unimodular'}
    if p != [p[-1]*v for v in reversed(p)]:
        return {'excluded': True, 'reason': 'full_matrix_not_reciprocal_up_to_sign'}
    return {'excluded': False, 'reason': 'not_excluded', 'realizability_claimed': False}


def case_ii_intersecting(n, m, p, q, upper, *, exit_parity=True):
    """Intersecting mixed/petal primary curves, generalized draft Appendix A.

    The caller has established Case II: m filament vertices; mixed length p;
    pure-petal primary length q; p+q+r0 <= n, q <= r0, m <= n. These are
    geometric preconditions, not inferred from a candidate polynomial.
    If p+2q=n, r0=q is forced. If U**q<2, any additional petal-only edge is
    impossible below U; floral exit parity then forces a second mixed pair
    with length <= m+q. Both tests use strict negativity, keeping equality.
    """
    _integers(n, m, p, q)
    if isinstance(upper, float):
        raise TypeError('use an exact rational upper bound')
    upper = Fraction(upper)
    if upper <= 1 or type(exit_parity) is not bool:
        raise ValueError('invalid upper bound or parity switch')
    if m > n or p+2*q > n:
        return {'excluded': True, 'reason': 'integer_primary_budget'}
    t = 1/upper
    base = (1-t**m)*(1-t**q)-2*t**p
    if base < 0:
        return {'excluded': True, 'reason': 'mixed_pair_growth', 'rational_value': str(base)}
    if exit_parity and p+2*q == n and upper**q < 2:
        strengthened = base-2*t**(m+q)
        if strengthened < 0:
            return {'excluded': True, 'reason': 'saturated_petal_exit_parity',
                    'forced_petals': q, 'extra_mixed_pair_max_length': m+q,
                    'rational_value': str(strengthened)}
    return {'excluded': False, 'reason': 'not_excluded', 'realizability_claimed': False}


def case_ii_integer_sieve(n, upper, *, exit_parity=True):
    """Exhaust this Case II integer relaxation, using monotone m-block pruning.

    Both growth expressions increase with m, so excluding m=n excludes the
    whole interval 1 <= m <= n. Other m values are inspected only when the
    block bound fails. This is a branch certificate, not a disc minimum.
    """
    _integers(n)
    counts, survivors, block_witnesses = Counter(), [], []
    for q in range(1, (n-1)//2+1):
        for p in range(1, n-2*q+1):
            counts['integer_triples_covered'] += n
            largest = case_ii_intersecting(n, n, p, q, upper, exit_parity=exit_parity)
            counts['growth_evaluations'] += 1
            if largest['excluded']:
                counts['whole_m_blocks_excluded'] += 1
                counts['triples_excluded'] += n
                block_witnesses.append({'p': p, 'q': q, 'm_interval': [1, n], **largest})
                continue
            survivors.append([n, p, q])
            for m in range(1, n):
                result = case_ii_intersecting(n, m, p, q, upper, exit_parity=exit_parity)
                counts['growth_evaluations'] += 1
                if result['excluded']:
                    counts['triples_excluded'] += 1
                else:
                    survivors.append([m, p, q])
    return {'scope': 'Case II intersecting mixed/petal primary curves only',
            'n': n, 'upper': str(Fraction(upper)), 'completed_integer_relaxation': True,
            'exit_parity_used': exit_parity, 'counts': dict(counts),
            'survivors': sorted(survivors), 'block_witnesses': block_witnesses,
            'global_minimum_proved': False, 'survivors_claimed_realizable': False}
