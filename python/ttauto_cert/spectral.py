# SPDX-License-Identifier: GPL-3.0-or-later
"""Exact arithmetic migrated from the large-disc and genus-three proof jobs.

All polynomial lists are descending. Catalogues are necessary supersets;
neither a spectrum nor a compatible orbit model is a surface realization.
"""
from fractions import Fraction
from functools import lru_cache
from math import gcd
from .lefschetz import traces


def multiply(a, b):
    result = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return result


def opposite(p):
    return [c*(-1)**i for i, c in enumerate(p)]


def schur_inside(p, radius):
    """Strict Schur recursion, with integer scaling and content removal."""
    r = Fraction(radius)
    if r <= 0:
        raise ValueError('positive disk radius required')
    degree = len(p)-1
    a = [c*r.numerator**(degree-i)*r.denominator**i for i, c in enumerate(p)]
    while len(a) > 1:
        first, last = a[0], a[-1]
        if abs(last) >= abs(first):
            return False
        b = [first*a[i]-last*a[-1-i] for i in range(len(a)-1)]
        divisor = 0
        for c in b:
            divisor = gcd(divisor, abs(c))
        if not divisor:
            return False
        a = [c//divisor for c in b]
    return bool(a[0])


def trace_bounds(rank, upper, horizon):
    if type(rank) is not int or rank < 2 or rank % 2:
        raise ValueError('positive even reciprocal degree required')
    r = Fraction(upper)
    if r <= 1:
        raise ValueError('upper must exceed one')
    return [int((rank//2)*(r**k+r**(-k))) for k in range(1, horizon+1)]


def coefficient_range(coeff, powers, bounds):
    k = len(coeff)
    known = -sum(coeff[j]*powers[k-j] for j in range(1, k))
    lo = -(-(known-bounds[k-1])//k)
    hi = (known+bounds[k-1])//k
    return lo, hi, known


def reciprocal_polynomials(rank, upper):
    """Full finite Newton containing catalogue, including nonexpanding spectra."""
    half = rank//2
    bounds = trace_bounds(rank, upper, half)
    def visit(coeff, powers):
        if len(coeff) == half+1:
            yield coeff+coeff[-2::-1]
            return
        lo, hi, known = coefficient_range(coeff, powers, bounds)
        k = len(coeff)
        for c in range(lo, hi+1):
            yield from visit(coeff+[c], powers+[known-k*c])
    yield from visit([1], [rank])


@lru_cache(maxsize=16)
def stable_catalogue(rank, upper):
    return tuple(tuple(p) for p in reciprocal_polynomials(rank, upper)
                 if schur_inside(p, upper))


def ordinary_orbits(fixed):
    orbits = []
    for k, count in enumerate(fixed, 1):
        if type(count) is not int or count < 0:
            return None, {'period': k, 'reason': 'negative_or_nonintegral_fixed_count'}
        remaining = count-sum(d*orbits[d-1] for d in range(1, k) if k % d == 0)
        if remaining < 0 or remaining % k:
            return None, {'period': k, 'reason': 'negative_or_nonintegral_orbits'}
        orbits.append(remaining//k)
    return orbits, None


def peripheral_obstruction(p, dimension, horizon):
    """Necessary trace/orbit bounds for chi_M=chi_H*chi_B, dim B=dimension."""
    t = traces(p, horizon)[1:]
    formal = []
    for k, value in enumerate(t, 1):
        net = value-sum(d*formal[d-1] for d in range(1, k) if k % d == 0)
        if net % k:
            raise ArithmeticError('integral characteristic polynomial has fractional formal orbits')
        formal.append(net//k)
        bound = (dimension if k % 2 else 2*dimension)//k
        if value < -dimension or formal[-1] < -bound:
            return {'period': k, 'trace': value, 'formal_orbits': formal[-1],
                    'peripheral_dimension': dimension, 'orbit_lower_bound': -bound}
    return None
