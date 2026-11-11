import itertools
import unittest
from fractions import Fraction as Q

from ttauto_cert.exact_math import (
    charpoly, compare_largest_real_root, compare_largest_real_roots, compare_perron_to_rational,
    count_real_roots, identity, isolate_largest_real_root,
    matrix_multiply, polynomial_value,
)


def multiply_poly(a,b):
    result=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            result[i+j]+=x*y
    return result


def permutation_determinant(matrix):
    result=0
    for perm in itertools.permutations(range(len(matrix))):
        inversions=sum(perm[i]>perm[j] for i in range(len(perm)) for j in range(i+1,len(perm)))
        term=(-1)**inversions
        for i,j in enumerate(perm):
            term*=matrix[i][j]
        result+=term
    return result


class ExactMathTests(unittest.TestCase):
    def test_matrix_products_and_charpoly_against_determinants(self):
        for size in range(1,6):
            matrix=[[2*i-3*j+i*j for j in range(size)] for i in range(size)]
            self.assertEqual(matrix_multiply(identity(size),matrix),matrix)
            coefficients=charpoly(matrix)
            for x in range(-3,5):
                shifted=[[int(i==j)*x-matrix[i][j] for j in range(size)] for i in range(size)]
                self.assertEqual(polynomial_value(coefficients,x),permutation_determinant(shifted))
        huge=10**100
        self.assertEqual(matrix_multiply([[huge]],[[huge]]),[[10**200]])
        self.assertEqual(charpoly([[huge,1],[1,huge]]),[1,-2*huge,huge*huge-1])
        self.assertEqual(charpoly([]),[1])

    def test_repeated_roots_and_every_endpoint_convention(self):
        # (-2)^multiplicity 3, 0^multiplicity 2, 1^multiplicity 4.
        poly=[1]
        for root in [-2,-2,-2,0,0,1,1,1,1]:
            poly=multiply_poly(poly,[1,-root])
        distinct={-2,0,1}
        self.assertEqual(count_real_roots(poly),3)
        for left in range(-3,3):
            for right in range(left,3):
                for closed_left,closed_right in itertools.product([False,True],repeat=2):
                    expected=sum((root>left or closed_left and root==left) and
                                 (root<right or closed_right and root==right) for root in distinct)
                    self.assertEqual(count_real_roots(poly,left,right,include_left=closed_left,
                                                     include_right=closed_right),expected)
        self.assertEqual(compare_largest_real_root(poly,1),0)
        self.assertEqual(isolate_largest_real_root(poly),(Q(1),Q(1)))

    def test_known_dilatations_and_interval_isolation(self):
        torus=[1,-3,1]
        left,right=isolate_largest_real_root(torus,bits=100)
        self.assertLessEqual(right-left,Q(1,2**100))
        self.assertEqual(count_real_roots(torus,left,right),1)
        self.assertLess(polynomial_value(torus,left),0)
        self.assertGreater(polynomial_value(torus,right),0)
        self.assertEqual(compare_perron_to_rational([[2,1],[1,1]],Q(2618,1000)),1)
        self.assertEqual(compare_perron_to_rational([[2,1],[1,1]],Q(2619,1000)),-1)
        self.assertEqual(compare_perron_to_rational([[0,2],[2,0]],2),0)
        self.assertEqual(compare_perron_to_rational([[0,1],[0,0]],0),0)
        self.assertEqual(compare_perron_to_rational([[2,0],[0,3]],Q(5,2)),1)
        # Archived primitive six-puncture counterexample polynomial.
        bad=multiply_poly([1,-1],[1,-1,0,-4,0,-1,1])
        self.assertEqual(compare_largest_real_root(bad,Q(2015357181280,10**12)),1)
        self.assertEqual(compare_largest_real_root(bad,Q(2015357181281,10**12)),-1)

    def test_negative_roots_close_roots_and_no_real_roots(self):
        negative_left,negative_right=isolate_largest_real_root([1,3,2])
        self.assertLessEqual(negative_left,-1)
        self.assertGreaterEqual(negative_right,-1)
        self.assertEqual(count_real_roots([1,3,2],negative_left,negative_right,include_left=True),1)
        eps=Q(1,2**120)
        close=multiply_poly([1,-1],[1,-1-eps])
        left,right=isolate_largest_real_root(close,bits=20)
        self.assertEqual(count_real_roots(close,left,right,include_left=True),1)
        self.assertLessEqual(left,1+eps)
        self.assertGreaterEqual(right,1+eps)
        self.assertEqual(count_real_roots([1,0,1]),0)
        with self.assertRaises(ValueError):
            isolate_largest_real_root([1,0,1])
        with self.assertRaises(ValueError):
            compare_largest_real_root([5],1)
        with self.assertRaises(ValueError):
            count_real_roots([0])
        with self.assertRaises(TypeError):
            compare_largest_real_root([1,-3,1],2.618)
        with self.assertRaises(ValueError):
            compare_perron_to_rational([[1,-1],[0,1]],1)

    def test_exact_algebraic_root_comparison(self):
        polynomial=[1,-3,1]
        self.assertEqual(compare_largest_real_roots(polynomial,multiply_poly(polynomial,[1,-1])),0)
        self.assertEqual(compare_largest_real_roots(polynomial,multiply_poly(polynomial,[1,-5])),-1)
        self.assertEqual(compare_largest_real_roots([1,-1-Q(1,2**120)],[1,-1]),1)
        self.assertEqual(compare_largest_real_roots([-1,3,-1],polynomial),0)


if __name__=="__main__":
    unittest.main()
