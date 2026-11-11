import unittest
from fractions import Fraction
from itertools import product
from ttauto_cert import lefschetz as lf
from ttauto_cert.exact_math import count_real_roots


class LefschetzTests(unittest.TestCase):
    def test_self_checks_remain_active_under_optimization(self):
        lf.self_test()

    def test_newton_enumeration_matches_independent_coefficient_box(self):
        # Genus two, U=2: the independent Vieta box is |a|<=8, |b|<=24.
        U=Fraction(2)
        expected=[]
        for a,b in product(range(-8,9),range(-24,25)):
            p=[1,a,b,a,1]
            trace=lf.traces(p,12)
            if any(abs(trace[k])>2*(U**k+U**(-k)) for k in range(1,13)):continue
            if not count_real_roots(p,1,U):continue
            if count_real_roots(p,U,None) or count_real_roots(p,None,-U,include_right=False):continue
            expected.append(p)
        self.assertEqual(sorted(lf.reciprocal_superset(2,U,12)),sorted(expected))

    def test_known_lower_bounds_and_unresolved_status(self):
        result=lf.run(include_star=True)
        for name in ['s1','s2','s3','s6']:
            self.assertTrue(result[name]['lower_bound_established'])
            self.assertFalse(result[name]['realizability_claimed'])
            self.assertEqual(result[name]['survivors'],[])
        bad=[1,-3,2,0,2,-3,1]
        for name in ['s4','s5']:
            self.assertFalse(result[name]['lower_bound_established'])
            self.assertEqual([s['polynomial'] for s in result[name]['survivors']],[bad])

    def test_invalid_exact_inputs_fail_explicitly(self):
        with self.assertRaises(TypeError):list(lf.reciprocal_superset(2,1.9))
        with self.assertRaises(ValueError):list(lf.reciprocal_superset(0,Fraction(2)))
        with self.assertRaises(ValueError):lf.traces([1,1.0,1],8)
        with self.assertRaises(ValueError):lf.compatible([1,-3,1],[(4,1)])
        with self.assertRaises(ValueError):list(lf.singular_patterns([(5,1)],False,8))
        with self.assertRaises(ValueError):list(lf.partitions(2,0))

    def test_marked_unbranched_fixed_point_is_not_a_fixed_puncture(self):
        # One free phi 2-cycle fixed setwise by tau gives two psi fixed points.
        positive={'regular_orbits':[0,1,0,0]}
        negative={'regular_orbits':[2,0,0,0]}
        self.assertTrue(lf.equivariant_regular_compatible(positive,negative,0,marked_free_pair=True))
        self.assertFalse(lf.equivariant_regular_compatible({'regular_orbits':[0,0,0,1]},
                         {'regular_orbits':[0,0,0,1]},0,marked_free_pair=True))


if __name__=='__main__':unittest.main()
