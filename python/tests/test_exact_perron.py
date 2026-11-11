from dataclasses import FrozenInstanceError
from fractions import Fraction as Q
import itertools
import unittest

from ttauto_cert.exact_math import polynomial_value
from ttauto_cert.exact_perron import (
    PerronCertificate, compare_perron_entries, compare_perron_roots,
    compare_perron_to_rational, perron_certificate,
)


def _determinant(matrix):
    result = 0
    for permutation in itertools.permutations(range(len(matrix))):
        inversions = sum(a > b for i, a in enumerate(permutation)
                         for b in permutation[i + 1:])
        term = (-1) ** inversions
        for row, column in enumerate(permutation):
            term *= matrix[row][column]
        result += term
    return result


class ExactPerronTests(unittest.TestCase):
    def test_adjugate_polynomials_against_independent_cofactor_determinants(self):
        for size in range(1, 6):
            matrix = [[1 + (i * j + i + 2 * j) % 5 for j in range(size)]
                      for i in range(size)]
            certificate = PerronCertificate(matrix)
            for value in (-3, 0, 1, 4):
                shifted = [[int(i == j) * value - matrix[i][j]
                            for j in range(size)] for i in range(size)]
                self.assertEqual(polynomial_value(certificate.characteristic_polynomial, value),
                                 _determinant(shifted))
                for row, polynomial in enumerate(certificate.right_vector_polynomials):
                    # adj(B)[row,0] = cofactor(B)[0,row].
                    minor = [[shifted[i][j] for j in range(size) if j != row]
                             for i in range(1, size)]
                    self.assertEqual(polynomial_value(polynomial, value),
                                     (-1) ** row * _determinant(minor))
                    self.assertEqual(certificate.sign_polynomial(polynomial), 1)

    def test_imprimitive_and_singleton_matrices(self):
        certificate = perron_certificate([[0, 2], [1, 0]])
        self.assertEqual(certificate.compare_entries(0, 1), 1)
        self.assertEqual(certificate.compare_entries(1, 0), -1)
        self.assertEqual(certificate.compare_entries(1, 1), 0)
        self.assertEqual(certificate.compare_root_rational(Q(7, 5)), 1)
        self.assertEqual(certificate.compare_root_rational(Q(3, 2)), -1)
        self.assertEqual(compare_perron_to_rational([[0, 3], [3, 0]], 3), 0)
        for entry in (0, 1, 17):
            one = PerronCertificate([[entry]])
            self.assertEqual(one.right_vector_polynomials, ((1,),))
            self.assertEqual(one.compare_root_rational(entry), 0)
            self.assertEqual(one.compare_entries(0, 0), 0)
        cycle = PerronCertificate([[0, 2, 0], [0, 0, 3], [4, 0, 0]])
        self.assertEqual(cycle.compare_entries(0, 1), -1)
        self.assertEqual(cycle.compare_entries(1, 2), 1)

    def test_exact_equal_coordinates_with_nonzero_difference_polynomial(self):
        # rho=1+sqrt(3), and v_0=v_1. The characteristic polynomial is
        # x (x^2-2x-2), so equality requires the root's proper gcd factor.
        certificate = PerronCertificate([[1, 1, 1], [1, 1, 1], [1, 1, 0]])
        self.assertNotEqual(certificate.right_vector_polynomials[0],
                            certificate.right_vector_polynomials[1])
        self.assertEqual(certificate.compare_entries(0, 1), 0)
        self.assertEqual(certificate.compare_entries(0, 2), 1)
        self.assertEqual(certificate.sign_polynomial([1, -2, -2]), 0)
        self.assertEqual(certificate.sign_polynomial([1, 0]), 1)
        self.assertEqual(certificate.sign_polynomial([1, -2, -2, 0]), 0)
        self.assertEqual(certificate.sign_polynomial([0, 0]), 0)
        # Repeated non-Perron roots must be square-free reduced correctly.
        repeated = PerronCertificate([[2, 1, 1], [1, 2, 1], [1, 1, 2]])
        for left in range(3):
            for right in range(3):
                self.assertEqual(repeated.compare_entries(left, right), 0)

    def test_shared_factor_is_not_mistaken_for_equal_perron_roots(self):
        three = [[2, 1], [1, 2]]
        five = [[3, 2], [2, 3]]
        self.assertEqual(compare_perron_roots(three, five), -1)
        self.assertEqual(compare_perron_roots(five, three), 1)
        self.assertEqual(compare_perron_roots(three, [[1, 1, 1]] * 3), 0)
        self.assertEqual(PerronCertificate(three).compare_root(PerronCertificate(three)), 0)

    def test_huge_integers_and_distinct_roots_closer_than_double_precision(self):
        huge = 2 ** 80
        smaller = PerronCertificate([[huge, 1], [1, 0]])
        larger = PerronCertificate([[huge, 1], [2, 0]])
        # Roots differ by approximately 1/huge; doubles identify both with huge.
        self.assertEqual(smaller.compare_root(larger), -1)
        self.assertEqual(smaller.compare_root_rational(huge), 1)
        self.assertEqual(smaller.compare_root_rational(Q(2 * huge * huge + 3, 2 * huge)), -1)
        self.assertEqual(smaller.sign_polynomial([1, -huge]), 1)
        self.assertEqual(smaller.sign_polynomial([-1, huge]), -1)
        # Coordinates differ by a relative amount below 2^-80.
        self.assertEqual(compare_perron_entries([[huge, huge], [huge, huge + 1]], 0, 1), -1)
        self.assertEqual(compare_perron_entries([[huge + 1, huge], [huge, huge]], 0, 1), 1)

    def test_positive_integer_eigenvectors_from_exact_diagonal_conjugation(self):
        # Independently prescribed PF data, including repeated coordinates:
        # if C has constant row sum L, A=D C D^-1 has right vector diag(D).
        for diagonal in ((1, 2, 3), (3, 3, 1), (7, 1, 7, 2)):
            common = 1
            for entry in diagonal:
                common *= entry
            size = len(diagonal)
            constant_rows = [[common * (1 + (j - i) % size)
                              for j in range(size)] for i in range(size)]
            matrix = [[constant_rows[i][j] * diagonal[i] // diagonal[j]
                       for j in range(size)] for i in range(size)]
            certificate = PerronCertificate(matrix)
            self.assertEqual(certificate.compare_root_rational(sum(constant_rows[0])), 0)
            for left, right in itertools.product(range(size), repeat=2):
                expected = (diagonal[left] > diagonal[right]) - (diagonal[left] < diagonal[right])
                self.assertEqual(certificate.compare_entries(left, right), expected)

    def test_input_is_copied_and_certificate_is_immutable(self):
        matrix = [[2, 1], [1, 1]]
        certificate = PerronCertificate(matrix)
        matrix[0][0] = -100
        self.assertEqual(certificate.matrix, ((2, 1), (1, 1)))
        with self.assertRaises(FrozenInstanceError):
            certificate.matrix = ((0,),)
        left, right = certificate.root_interval(100)
        self.assertLessEqual(right - left, Q(1, 2 ** 100))

    def test_strict_matrix_and_argument_validation(self):
        for matrix in ([], [[1, 2]], [[1, 2], [3]], [[-1]],
                       [[1, 0], [0, 1]], [[1, 1], [0, 1]]):
            with self.subTest(matrix=matrix), self.assertRaises(ValueError):
                PerronCertificate(matrix)
        for matrix in ([[True]], [[1.0]], [[Q(1)]], '1', [1], [[False, 1], [1, 0]]):
            with self.subTest(matrix=matrix), self.assertRaises(TypeError):
                PerronCertificate(matrix)
        certificate = PerronCertificate([[2, 1], [1, 1]])
        for index in (True, 0.0, -1, 2):
            with self.assertRaises(ValueError):
                certificate.compare_entries(index, 0)
        for value in (True, 1.0, '1'):
            with self.assertRaises(TypeError):
                certificate.compare_root_rational(value)
            with self.assertRaises(TypeError):
                certificate.sign_polynomial([value])
        for bits in (True, -1, 1.0):
            with self.assertRaises(ValueError):
                certificate.root_interval(bits)
        with self.assertRaises(ValueError):
            certificate.sign_polynomial([])


if __name__ == '__main__':
    unittest.main()
