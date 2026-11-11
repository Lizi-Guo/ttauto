# SPDX-License-Identifier: GPL-3.0-or-later
"""Sphere deck tests and closed-surface invariant/anti-invariant coupling.

Ported from compute_n7_lower_20260915.py, branch_coupling_exact.py and
anti8_base_coupling_diagnostic_20261010.py. No server paths or classifiers.
"""
from collections import Counter
from itertools import product
from .lefschetz import (partitions, compatible, branch_singular_compatible,
                       equivariant_regular_compatible, traces)
from .hyperelliptic_mod2 import branch_permutation_charpoly, _bits
from .spectral import ordinary_orbits


class ModelBudgetExceeded(RuntimeError):
    pass


class Budget:
    def __init__(self, limit):
        self.limit, self.used = limit, 0

    def tick(self):
        if self.limit is not None and self.used >= self.limit:
            raise ModelBudgetExceeded('orbit-model budget exhausted; no exclusion')
        self.used += 1


def disc_cover(n, interior, outer):
    groups, branch_indices, colors = [], [], [n]
    for p, number in sorted(Counter(interior).items()):
        if p % 2:
            branch_indices.append(len(groups))
            groups.append([2*p, number]); colors.append(number)
        else:
            groups.append([p, 2*number])
    if outer >= 3:
        if outer % 2:
            branch_indices.append(len(groups)); groups.append([2*outer, 1])
        else:
            groups.append([outer, 2])
    branch_count = sum(colors)+outer % 2
    genus = (branch_count-2)//2
    if branch_count % 2 or sum((p-2)*number for p, number in groups) != 4*genus-4:
        raise ValueError('invalid orienting-cover Euler data')
    return dict(genus=genus, singular_groups=groups, branch_singular_groups=branch_indices,
                branch_colors=colors, fixed_branch_points=outer % 2,
                branch_regular=n+int(outer == 1), branch_fixed_regular=int(outer == 1),
                marked_free_regular_pair=outer == 2)


def disc_patterns(p, case, horizon, budget=None):
    cover, outer = case['cover'], case['outer_prongs']
    marked = lambda k: 1 if outer == 1 else 2 if outer == 2 and k % 2 == 0 else 0
    return [compatible(p, cover['singular_groups'], sign, marked, horizon, budget=budget)
            for sign in (False, True)]


def branch_regular_compatible(a, b, lengths):
    counts = Counter(lengths)
    stop = len(a['regular_orbits'])
    left = [0]+[v-counts[k] for k, v in enumerate(a['regular_orbits'], 1)]
    right = [0]+[v-counts[k] for k, v in enumerate(b['regular_orbits'], 1)]
    if any(v < 0 for v in left+right):
        return False
    for k in range(1, stop+1, 2):
        if left[k] % 2 or right[k] % 2:
            return False
        if 2*k <= stop:
            x, y = left[2*k]-right[k]//2, right[2*k]-left[k]//2
            if x < 0 or x != y or x % 2:
                return False
    return all(left[k] == right[k] for k in range(4, stop+1, 4))


def disc_deck_feasible(p, case, positive, negative, budget, couple_mod2=True):
    cv = case['cover']
    for a in positive:
        for b in negative:
            budget.tick()
            if not branch_singular_compatible(a, b, cv['singular_groups'], cv['branch_singular_groups']):
                continue
            if not equivariant_regular_compatible(a, b, cv['branch_regular'],
                                                  cv['branch_fixed_regular'], cv['marked_free_regular_pair']):
                continue
            # Couple the SAME branch permutation to the regular deck orbits
            # and mod-2 action, rather than checking their existence separately.
            if not couple_mod2:
                return {'positive': a, 'negative': b}
            singular = tuple(length for i in cv['branch_singular_groups']
                             for length, shift in a['singular_orbits'][i])
            for part in partitions(case['n']):
                budget.tick()
                regular = part+((1,) if case['outer_prongs'] == 1 else ())
                if not branch_regular_compatible(a, b, regular):
                    continue
                if _bits(branch_permutation_charpoly(regular+singular)) == _bits(p):
                    return {'positive': a, 'negative': b, 'regular_branch_cycles': list(regular)}
    return None


def joint_group(prongs, number, horizon):
    """Correlated (S,Z,E) vectors; preserve correlations between all periods."""
    for lengths in partitions(number):
        for rotations in product(range(prongs), repeat=len(lengths)):
            vector = []
            for k in range(1, horizon+1):
                s = z = e = 0
                for length, rotation in zip(lengths, rotations):
                    if k % length:
                        continue
                    r = rotation*(k//length) % prongs
                    z += length
                    if prongs % 2:
                        s += length*prongs if r == 0 else 0
                    else:
                        s += length*(prongs-1 if r == 0 else -1 if r % 2 == 0 else 1)
                        e += length*(1 if r % 2 == 0 else -1)
                vector.append((s, z, e))
            yield vector


def joint_models(orders, horizon, budget):
    groups = sorted(Counter(p+2 for p in orders).items())
    def visit(index, vector):
        if index == len(groups):
            budget.tick()
            yield vector
            return
        p, n = groups[index]
        for row in joint_group(p, n, horizon):
            budget.tick()
            joined = [tuple(a+b for a, b in zip(x, y)) for x, y in zip(vector, row)]
            yield from visit(index+1, joined)
    yield from visit(0, [(0, 0, 0)]*horizon)


def regular_from_anti(anti, joint):
    fixed = [t-v[0] for t, v in zip(anti, joint)]
    orbits, failure = ordinary_orbits(fixed)
    return (fixed, orbits), failure


def coupled_orbits(anti, base, joint, regular):
    fixed, orbits = regular
    negative_fixed, negative_orbits = [], []
    for k, (ta, tb, (s, z, e)) in enumerate(zip(anti, base, joint), 1):
        numerator = 2-tb-z+ta+e
        if numerator % 2:
            return None, {'period': k, 'reason': 'nonintegral_negative_fixed_count'}
        value = numerator//2
        if not 0 <= value <= fixed[k-1]:
            return None, {'period': k, 'reason': 'negative_fixed_count_out_of_bounds'}
        negative_fixed.append(value)
        # A negative return sign remains negative only at ODD returns.
        net = value-sum(d*negative_orbits[d-1] for d in range(1, k)
                        if k % d == 0 and (k//d) % 2)
        if net % k or not 0 <= net//k <= orbits[k-1]:
            return None, {'period': k, 'reason': 'negative_orbits_out_of_bounds'}
        negative_orbits.append(net//k)
    return {'regular_orbits': orbits, 'negative_orbits': negative_orbits,
            'negative_fixed_counts': negative_fixed}, None
