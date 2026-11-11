import unittest
from ttauto_cert.exact_math import charpoly
from ttauto_cert.hyperelliptic_mod2 import (
    branch_cycle_types, branch_permutation_charpoly,
    compatible_branch_cycle_types, polynomial_divides_mod2)


def explicit_quotient_matrix(cycles):
    # Independent linear construction; no polynomial multiplication or division.
    permutation = []
    for length in cycles:
        start = len(permutation)
        permutation.extend(start + (i + 1) % length for i in range(length))
    n = len(permutation)
    d = n - 2
    matrix = [[0] * d for _ in range(d)]
    for j in range(d):
        ambient = [0] * n
        ambient[permutation[j]] ^= 1
        ambient[permutation[n - 1]] ^= 1
        # Basis e_i+e_(n-1), i=0..n-3. The omitted basis vector is the
        # sum of those retained, since we quotient by the all-ones vector.
        for i in range(d):
            matrix[i][j] = ambient[i] ^ ambient[n - 2]
    return matrix


class HyperellipticMod2Tests(unittest.TestCase):
    def test_exact_divisibility_in_characteristic_two(self):
        self.assertTrue(polynomial_divides_mod2([1, 1], [1, 0, 0, 0, 1]))
        self.assertTrue(polynomial_divides_mod2([1, 0, 1], [1, 0, 0, 0, 1]))
        self.assertTrue(polynomial_divides_mod2([1, 1, 1], [1, 0, 0, 1]))
        self.assertFalse(polynomial_divides_mod2([1, 0, 1], [1, 0, 0, 1]))
        self.assertFalse(polynomial_divides_mod2([1, 1, 1], [1, 0, 1]))
        self.assertTrue(polynomial_divides_mod2([1, 1], [0]))
        with self.assertRaises(ZeroDivisionError):
            polynomial_divides_mod2([2], [1, 1])

    def test_cycle_formula_against_explicit_module(self):
        checked = 0
        # Every uncolored cycle type on N=4,6,8,10 branch points.
        for n in (4, 6, 8, 10):
            for row in branch_cycle_types((n,)):
                expected = [coefficient % 2 for coefficient in
                            charpoly(explicit_quotient_matrix(row['branch_cycles']))]
                self.assertEqual(row['mod2_polynomial'], expected)
                self.assertTrue(row['exact_division_verified'])
                checked += 1
        self.assertEqual(checked, 80)  # p(4)+p(6)+p(8)+p(10)=5+11+22+42.

    def test_colored_groups_and_fixed_point_are_preserved(self):
        rows = list(branch_cycle_types((6, 3), fixed_points=1))
        self.assertEqual(len(rows), 33)
        self.assertEqual(len({tuple(row['mod2_polynomial']) for row in rows}), 6)
        for row in rows:
            self.assertEqual([sum(group) for group in row['group_cycles']], [6, 3])
            self.assertEqual(row['fixed_points'], 1)
            self.assertEqual(row['branch_cycles'][-1], 1)
        self.assertEqual(len(list(branch_cycle_types((), fixed_points=4))), 1)

    def test_n6_rejected_residue_and_known_retained_pattern(self):
        rejected = [1, -4, 6, -8, 11, -8, 6, -4, 1]
        self.assertEqual(compatible_branch_cycle_types(rejected, (6, 3), 1), [])
        # The same obstruction holds with unrestricted branch permutations.
        self.assertEqual(compatible_branch_cycle_types(rejected, (10,)), [])
        retained = [1, -5, 8, -3, -2, -3, 8, -5, 1]
        self.assertTrue(compatible_branch_cycle_types(retained, (6, 3), 1))

    def test_p0_for_two_three_prong_singularities(self):
        p0 = [1, -3, 2, 0, 2, -3, 1]
        rows = compatible_branch_cycle_types(p0, (6, 2))
        self.assertEqual(len(rows), 2)
        self.assertEqual({tuple(row['group_cycles'][0]) for row in rows}, {(1, 5)})
        self.assertEqual({tuple(row['group_cycles'][1]) for row in rows}, {(1, 1), (2,)})

    def test_actual_torus_characteristic_polynomials(self):
        for a in range(1, 5):
            for b in range(1, 5):
                # A=[[1+ab,a],[b,1]], tau=-I; the origin is a fixed branch point.
                self.assertTrue(compatible_branch_cycle_types([1, -(2 + a*b), 1], (3,), 1))

    def test_invalid_inputs(self):
        for cycles in ((1, 2), (0, 4), (True, 3)):
            with self.assertRaises(ValueError):
                branch_permutation_charpoly(cycles)
        with self.assertRaises(ValueError):
            list(branch_cycle_types((3,), 0))
        with self.assertRaises(ValueError):
            list(branch_cycle_types((4,), -1))
        with self.assertRaises(ValueError):
            compatible_branch_cycle_types([1, -3.0, 1], (4,))
        with self.assertRaises(ValueError):
            compatible_branch_cycle_types([1, -3, 1], (6,))


if __name__ == '__main__':
    unittest.main()
