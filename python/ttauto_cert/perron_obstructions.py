"""Exact obstructions to a simple, strictly dominant positive real eigenvalue.

These are necessary conditions for the homology action of the positive lift of
an orientable pseudo-Anosov map. They must not classify arbitrary irreducible
nonnegative matrices, which can have equal-modulus eigenvalues. A negative
answer means only that these particular obstructions were not found.
"""
from math import gcd
from .exact_math import (
    _derivative, _division, _gcd, _trim, compare_largest_real_root,
    compare_largest_real_roots, count_real_roots)


def _integral(coefficients):
    if any(coefficient.denominator != 1 for coefficient in coefficients):
        raise ArithmeticError('a monic integer polynomial acquired a nonintegral factor')
    return [int(coefficient) for coefficient in coefficients]


def _exact_quotient(dividend, divisor):
    quotient, remainder = _division(dividend, divisor)
    if any(remainder):
        raise ArithmeticError('claimed polynomial factor does not divide exactly')
    return _integral(quotient)


def strict_perron_obstruction(polynomial):
    """Return a reproducible record with ``excluded`` and an exact reason.

    Input: a monic integer polynomial whose largest real root is greater than
    one. First test whether that root occurs in gcd(P,P'). Then remove all
    factors x-1 and x+1. If the remaining factor is Q(x^d), d>1, and has the
    same largest real root, multiplication by d-th roots of unity gives other
    roots of the same modulus. Neither test uses numerical root estimates.
    """
    polynomial = list(polynomial)
    if (len(polynomial) < 2 or polynomial[0] != 1 or
            any(type(coefficient) is not int for coefficient in polynomial)):
        raise ValueError('a monic nonconstant integer polynomial is required')
    if (not count_real_roots(polynomial) or
            compare_largest_real_root(polynomial, 1) <= 0):
        raise ValueError('a largest real root strictly greater than one is required')
    exact = _trim(polynomial)  # The internal gcd expects rational coefficients.
    derivative = _derivative(exact)
    common = _gcd(exact, derivative)
    if (len(common) > 1 and count_real_roots(common) and
            compare_largest_real_roots(common, polynomial) == 0):
        return {'excluded': True, 'reason': 'repeated_dominant_real_root',
                'polynomial': polynomial, 'gcd_with_derivative': _integral(common),
                'polynomial_divided_by_gcd': _exact_quotient(polynomial, common),
                'derivative_divided_by_gcd': _exact_quotient(derivative, common),
                'largest_real_root_equality_verified': True,
                'strict_dominance_or_simplicity_violated': True}

    remaining = polynomial
    stripped = []
    for root in (1, -1):
        factor = [1, -root]
        multiplicity = 0
        while len(remaining) > 1:
            quotient, remainder = _division(remaining, factor)
            if any(remainder):
                break
            remaining = _integral(quotient)
            multiplicity += 1
        if multiplicity:
            stripped.append({'factor': factor, 'multiplicity': multiplicity})
    rotation_order = 0
    for exponent, coefficient in enumerate(reversed(remaining)):
        if coefficient:
            rotation_order = gcd(rotation_order, exponent)
    if (rotation_order > 1 and count_real_roots(remaining) and
            compare_largest_real_roots(remaining, polynomial) == 0):
        compressed = remaining[::rotation_order]
        reconstructed = [0] * len(remaining)
        reconstructed[::rotation_order] = compressed
        if reconstructed != remaining:
            raise ArithmeticError('rotation-factor substitution identity failed')
        return {'excluded': True, 'reason': 'dominant_root_rotation_companions',
                'polynomial': polynomial, 'stripped_linear_factors': stripped,
                'remaining_factor': remaining, 'rotation_order': rotation_order,
                'compressed_polynomial': compressed,
                'substitution_identity_verified': True,
                'largest_real_root_equality_verified': True,
                'strict_dominance_or_simplicity_violated': True}
    # A rotating dominant factor may be hidden by other cyclotomic factors.
    # Write P(x)=sum_r x^r Q_r(x^d). A common divisor G of all Q_r gives
    # G(x^d)|P(x), exactly, without numerical roots or a factorization package.
    degree = len(polynomial)-1
    ascending = polynomial[::-1]
    for order in range(2, degree+1):
        residue_factors = []
        for residue in range(order):
            coefficients = ascending[residue::order][::-1]
            if any(coefficients):
                residue_factors.append(_trim(coefficients))
        common = residue_factors[0]
        for q in residue_factors[1:]:
            common = _gcd(common, q)
            if len(common) == 1:
                break
        if len(common) == 1:
            continue
        # Normalize a possibly uncombined first residue before integrality.
        common = [c/common[0] for c in common]
        compressed = _integral(common)
        inflated = [0]*((len(compressed)-1)*order+1)
        inflated[::order] = compressed
        if count_real_roots(inflated) and compare_largest_real_roots(inflated, polynomial) == 0:
            quotient = _exact_quotient(polynomial, inflated)
            return {'excluded': True, 'reason': 'dominant_rotating_factor',
                    'polynomial': polynomial, 'rotation_order': order,
                    'compressed_polynomial': compressed, 'rotating_factor': inflated,
                    'exact_quotient': quotient, 'largest_real_root_equality_verified': True}
    return {'excluded': False, 'reason': 'no_obstruction_found',
            'polynomial': polynomial, 'stripped_linear_factors': stripped,
            'remaining_factor': remaining,
            'strict_perron_property_certified': False,
            'realizability_claimed': False}
