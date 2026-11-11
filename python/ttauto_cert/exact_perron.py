"""Exact Perron roots and right-vector comparisons for irreducible matrices.

All polynomials have descending coefficients.  A certificate represents a
strictly positive right Perron vector by a column of ``adj(x I - A)``, evaluated
at the largest real root of the characteristic polynomial.  No factorization,
floating-point approximation, primitivity, or topological conclusion is used.

Equality of an algebraic expression is decided by polynomial gcd and Sturm
root counting *before* interval refinement.  Consequently equal coordinates
cannot cause an infinite precision loop.  See the accompanying research note
for the termination and positivity arguments.
"""
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache

from .exact_math import (
    _division, _gcd, _trim, charpoly, compare_largest_real_root,
    compare_largest_real_roots, count_real_roots, identity,
    isolate_largest_real_root, matrix_multiply, polynomial_value,
    sturm_sequence,
)


def _validated_matrix(matrix):
    if not isinstance(matrix, Sequence) or isinstance(matrix, (str, bytes)):
        raise TypeError('matrix must be a finite sequence of finite rows')
    if not matrix:
        raise ValueError('matrix must be nonempty')
    rows = []
    for row in matrix:
        if not isinstance(row, Sequence) or isinstance(row, (str, bytes)):
            raise TypeError('matrix rows must be finite sequences')
        if len(row) != len(matrix):
            raise ValueError('matrix must be square, without ragged rows')
        if any(type(entry) is not int for entry in row):
            raise TypeError('matrix entries must be exact Python integers, not bool or float')
        if any(entry < 0 for entry in row):
            raise ValueError('matrix entries must be nonnegative')
        rows.append(tuple(row))
    result = tuple(rows)
    # Strong connectivity is tested in both directions, without recursion.
    # The 1-by-1 zero matrix is irreducible under the usual PF convention.
    for transpose in (False, True):
        seen, pending = {0}, [0]
        while pending:
            vertex = pending.pop()
            for other in range(len(result)):
                entry = result[other][vertex] if transpose else result[vertex][other]
                if entry and other not in seen:
                    seen.add(other)
                    pending.append(other)
        if len(seen) != len(result):
            raise ValueError('matrix must be irreducible (strongly connected)')
    return result


def _rational(value):
    if type(value) is not int and not isinstance(value, Fraction):
        raise TypeError('use an exact Python integer or Fraction, not bool or float')
    return Fraction(value)


def _polynomial(coefficients):
    if not isinstance(coefficients, Sequence) or isinstance(coefficients, (str, bytes)):
        raise TypeError('polynomial coefficients must be a finite sequence')
    if not coefficients:
        raise ValueError('polynomial coefficient list must be nonempty')
    return _trim(tuple(_rational(value) for value in coefficients))


def _interval_value(coefficients, left, right):
    """An outward exact rational Horner enclosure; no rounding occurs."""
    lower = upper = Fraction(0)
    for coefficient in coefficients:
        products = (lower * left, lower * right, upper * left, upper * right)
        lower, upper = min(products) + coefficient, max(products) + coefficient
    return lower, upper


@lru_cache(maxsize=128)
def _certificate_data(matrix):
    coefficients = tuple(charpoly(matrix))
    size = len(matrix)
    auxiliary = identity(size)
    # adj(x I-A) = B_0 x^(n-1) + ... + B_(n-1),
    # B_0=I and B_k=A B_(k-1)+c_k I.
    column_coefficients = [[] for _ in matrix]
    for degree in range(size):
        for row in range(size):
            column_coefficients[row].append(auxiliary[row][0])
        auxiliary = matrix_multiply(matrix, auxiliary)
        for row in range(size):
            auxiliary[row][row] += coefficients[degree + 1]
    if any(entry for row in auxiliary for entry in row):
        raise ArithmeticError('adjugate recurrence failed its Cayley-Hamilton check')
    return coefficients, tuple(tuple(row) for row in column_coefficients)


@dataclass(frozen=True, init=False)
class PerronCertificate:
    """Immutable exact data for a nonnegative irreducible integer matrix.

    ``right_vector_polynomials[i]`` evaluated at ``rho(A)`` is its i-th
    positive right-vector coordinate; the common normalization is arbitrary.
    Public indices are zero-based. Comparisons return -1, 0, or +1.
    The constructor copies its input, so later caller mutation has no effect.
    """

    matrix: tuple
    characteristic_polynomial: tuple
    right_vector_polynomials: tuple

    def __init__(self, matrix):
        exact = _validated_matrix(matrix)
        coefficients, vector = _certificate_data(exact)
        object.__setattr__(self, 'matrix', exact)
        object.__setattr__(self, 'characteristic_polynomial', coefficients)
        object.__setattr__(self, 'right_vector_polynomials', vector)

    def root_interval(self, bits=80):
        """Closed rational isolating interval, of width at most 2**(-bits)."""
        if type(bits) is not int or bits < 0:
            raise ValueError('bits must be a nonnegative exact integer')
        return isolate_largest_real_root(self.characteristic_polynomial, bits)

    def sign_polynomial(self, coefficients):
        """Exact sign of a rational polynomial evaluated at this Perron root.

        This also handles reducible or nonsquare-free characteristic
        polynomials. A gcd root in the Perron isolating interval detects
        equality without finding the root's minimal polynomial.
        """
        expression = _polynomial(coefficients)
        defining = sturm_sequence(self.characteristic_polynomial)[0]
        expression = _division(expression, defining)[1]
        if expression == (0,):
            return 0
        if len(expression) == 1:
            return 1 if expression[0] > 0 else -1
        common = _gcd(defining, expression)
        left, right = self.root_interval(0)
        if (len(common) > 1 and
                count_real_roots(common, left, right, include_left=True)):
            return 0
        bits = 0
        while True:
            if left == right:
                value = polynomial_value(expression, left)
                if not value:
                    raise ArithmeticError('gcd equality check missed a rational root')
                return 1 if value > 0 else -1
            lower, upper = _interval_value(expression, left, right)
            if lower > 0:
                return 1
            if upper < 0:
                return -1
            bits = 4 if bits == 0 else 2 * bits
            left, right = self.root_interval(bits)

    def compare_entries(self, left, right):
        """Compare two coordinates of the positive right Perron vector."""
        for index in (left, right):
            if type(index) is not int or not 0 <= index < len(self.matrix):
                raise ValueError('vector index must be an exact integer in range')
        if left == right:
            return 0
        difference = tuple(a - b for a, b in zip(
            self.right_vector_polynomials[left], self.right_vector_polynomials[right]))
        return self.sign_polynomial(difference)

    def compare_root_rational(self, bound):
        """Compare rho(A) with an exact integer or Fraction."""
        return compare_largest_real_root(self.characteristic_polynomial, _rational(bound))

    def compare_root(self, other):
        """Compare rho(A) with another certificate's root or matrix's rho."""
        other = other if isinstance(other, PerronCertificate) else PerronCertificate(other)
        return compare_largest_real_roots(
            self.characteristic_polynomial, other.characteristic_polynomial)


def perron_certificate(matrix):
    return PerronCertificate(matrix)


def compare_perron_entries(matrix, left, right):
    return PerronCertificate(matrix).compare_entries(left, right)


def compare_perron_roots(left, right):
    return PerronCertificate(left).compare_root(right)


def compare_perron_to_rational(matrix, bound):
    """Irreducible-matrix counterpart of exact_math's broader PF comparison."""
    return PerronCertificate(matrix).compare_root_rational(bound)
