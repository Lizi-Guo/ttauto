import unittest
from ttauto_cert.perron_obstructions import strict_perron_obstruction


def multiply(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i+j] += a*b
    return result


class PerronObstructionTests(unittest.TestCase):
    def test_repeated_largest_real_root(self):
        quartic = [1, -1, -2, -1, 1]
        polynomial = [1, -2, -3, 2, 8, 2, -3, -2, 1]
        self.assertEqual(multiply(quartic, quartic), polynomial)
        result = strict_perron_obstruction(polynomial)
        self.assertTrue(result['excluded'])
        self.assertEqual(result['reason'], 'repeated_dominant_real_root')
        self.assertEqual(result['gcd_with_derivative'], quartic)
        self.assertEqual(result['polynomial_divided_by_gcd'], quartic)
        self.assertTrue(result['largest_real_root_equality_verified'])
        simple_square = strict_perron_obstruction([1, -4, 4])
        self.assertTrue(simple_square['excluded'])
        self.assertEqual(simple_square['gcd_with_derivative'], [1, -2])

    def test_four_masked_rotation_obstructions(self):
        for k in (4, 6, 8, 10):
            factor = [1, 0, 0, -k, 0, 0, 1]
            polynomial = [1, -2, 1, -k, 2*k, -k, 1, -2, 1]
            self.assertEqual(multiply([1, -2, 1], factor), polynomial)
            result = strict_perron_obstruction(polynomial)
            self.assertTrue(result['excluded'])
            self.assertEqual(result['reason'], 'dominant_root_rotation_companions')
            self.assertEqual(result['rotation_order'], 3)
            self.assertEqual(result['remaining_factor'], factor)
            self.assertEqual(result['compressed_polynomial'], [1, -k, 1])
            self.assertEqual(result['stripped_linear_factors'],
                             [{'factor': [1, -1], 'multiplicity': 2}])

    def test_both_linear_factors_and_unmasked_rotation(self):
        factor = [1, 0, -3, 0, 1]
        for extra in ([1], [1, -1], [1, 1], [1, 0, -1]):
            result = strict_perron_obstruction(multiply(extra, factor))
            self.assertTrue(result['excluded'])
            self.assertEqual(result['rotation_order'], 2)
            self.assertEqual(result['remaining_factor'], factor)

    def test_repeated_nondominant_real_root_is_not_rejected(self):
        polynomial = multiply([1, -3, 1], [1, -2, 1])
        result = strict_perron_obstruction(polynomial)
        self.assertFalse(result['excluded'])
        self.assertEqual(result['remaining_factor'], [1, -3, 1])
        self.assertFalse(result['strict_perron_property_certified'])

    def test_repeated_nonreal_roots_are_not_rejected(self):
        polynomial = multiply([1, -3, 1], multiply([1, 0, 1], [1, 0, 1]))
        self.assertFalse(strict_perron_obstruction(polynomial)['excluded'])

    def test_known_strict_examples_survive(self):
        for polynomial in ([1, -3, 1], [1, -1, -2, -1, 1],
                           [1, -5, 8, -3, -2, -3, 8, -5, 1]):
            result = strict_perron_obstruction(polynomial)
            self.assertFalse(result['excluded'])
            self.assertFalse(result['realizability_claimed'])

    def test_invalid_inputs_fail_closed(self):
        for polynomial in ([], [1], [2, -3, 1], [1, -3.0, 1],
                           [True, -3, 1], [1, 0, 1], [1, -2, 1]):
            with self.assertRaises(ValueError):
                strict_perron_obstruction(polynomial)


if __name__ == '__main__':
    unittest.main()
