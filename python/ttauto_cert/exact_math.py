"""Exact matrix and real-algebraic comparisons for folding certificates.

Polynomials use DESCENDING coefficients. Matrix entries are Python integers;
polynomial and root endpoints are integers or fractions.Fraction. Floats are
rejected at the exact-arithmetic boundary. A real-root interval may collapse
to a rational point. Repeated roots are counted once, via square-free reduction.
No numerical eigenvalue computation is used in an acceptance decision.
"""
from fractions import Fraction
from functools import lru_cache


def _q(value):
    if isinstance(value,float):
        raise TypeError("use Fraction or an integer, not an approximate float")
    return Fraction(value)


def identity(size):
    if not isinstance(size,int) or size<0:
        raise ValueError("matrix dimension must be a nonnegative integer")
    return [[int(i==j) for j in range(size)] for i in range(size)]


def _shape(matrix):
    rows=len(matrix)
    cols=len(matrix[0]) if rows else 0
    if any(len(row)!=cols for row in matrix):
        raise ValueError("ragged matrix")
    if any(not isinstance(a,int) for row in matrix for a in row):
        raise TypeError("matrix entries must be exact Python integers")
    return rows,cols


def matrix_multiply(left,right):
    m,k=_shape(left)
    k2,n=_shape(right)
    if k!=k2:
        raise ValueError("incompatible matrix dimensions")
    return [[sum(left[i][r]*right[r][j] for r in range(k))
             for j in range(n)] for i in range(m)]


def charpoly(matrix):
    """Faddeev--LeVerrier, with checked exact division at every coefficient."""
    n,m=_shape(matrix)
    if n!=m:
        raise ValueError("characteristic polynomial requires a square matrix")
    result=[1]
    auxiliary=identity(n)
    for k in range(1,n+1):
        auxiliary=matrix_multiply(matrix,auxiliary)
        trace=sum(auxiliary[i][i] for i in range(n))
        coefficient,remainder=divmod(-trace,k)
        if remainder:
            raise ArithmeticError("characteristic-polynomial division is not exact")
        result.append(coefficient)
        for i in range(n):
            auxiliary[i][i]+=coefficient
    return result


def _trim(poly):
    p=tuple(_q(a) for a in poly)
    if not p:
        raise ValueError("empty polynomial coefficient list")
    start=0
    while start<len(p)-1 and not p[start]:
        start+=1
    return p[start:]


def polynomial_value(coefficients,x):
    x=_q(x)
    value=Fraction(0)
    for coefficient in coefficients:
        value=value*x+_q(coefficient)
    return value


def _derivative(poly):
    degree=len(poly)-1
    return tuple(a*(degree-i) for i,a in enumerate(poly[:-1])) or (Fraction(0),)


def _division(dividend,divisor):
    divisor=_trim(divisor)
    if divisor==(0,):
        raise ZeroDivisionError("zero polynomial divisor")
    remainder=list(_trim(dividend))
    quotient=[Fraction(0)]*max(1,len(remainder)-len(divisor)+1)
    while len(remainder)>=len(divisor) and any(remainder):
        offset=len(remainder)-len(divisor)
        factor=remainder[0]/divisor[0]
        quotient[len(quotient)-1-offset]=factor
        for i,coefficient in enumerate(divisor):
            remainder[i]-=factor*coefficient
        remainder=list(_trim(remainder))
    return _trim(quotient),_trim(remainder)


def _monic(poly):
    return tuple(a/poly[0] for a in poly)


def _gcd(left,right):
    while right!=(0,):
        left,right=right,_division(left,right)[1]
    return _monic(left)


@lru_cache(maxsize=256)
def _sturm_cached(coefficients):
    p=_trim(coefficients)
    if p==(0,):
        raise ValueError("zero polynomial has no finite root set")
    if len(p)==1:
        return ((Fraction(1),),)
    divisor=_gcd(p,_derivative(p))
    p,remainder=_division(p,divisor)
    if any(remainder):
        raise ArithmeticError("square-free division failed")
    p=_monic(p)
    sequence=[p,_derivative(p)]
    while True:
        remainder=_division(sequence[-2],sequence[-1])[1]
        if not any(remainder):
            break
        # Positive rescaling controls coefficient growth and preserves signs.
        divisor=abs(remainder[0])
        sequence.append(tuple(-a/divisor for a in remainder))
    return tuple(sequence)


def sturm_sequence(coefficients):
    """Sturm chain of the monic square-free part (immutable exact tuples)."""
    return _sturm_cached(_trim(coefficients))


def _variations(sequence,x,positive_infinity=False):
    signs=[]
    for p in sequence:
        if x is None:
            value=p[0]
            if not positive_infinity and (len(p)-1)%2:
                value=-value
        else:
            value=polynomial_value(p,x)
        if value:
            signs.append(1 if value>0 else -1)
    return sum(a!=b for a,b in zip(signs,signs[1:]))


def _count(sequence,left,right,include_left=False,include_right=True):
    if left is not None and right is not None:
        if left>right:
            raise ValueError("left endpoint exceeds right endpoint")
        if left==right:
            return int(include_left and include_right and
                       polynomial_value(sequence[0],left)==0)
    count=_variations(sequence,left)-_variations(sequence,right,True)
    # With zero entries omitted, V(left)-V(right) counts (left,right].
    if include_left and left is not None and not polynomial_value(sequence[0],left):
        count+=1
    if not include_right and right is not None and not polynomial_value(sequence[0],right):
        count-=1
    return count


def count_real_roots(coefficients,left=None,right=None,*,include_left=False,include_right=True):
    """Count DISTINCT real roots in the selected interval; None denotes infinity.

    Default convention is (left,right]. Infinite endpoints are never roots.
    For a singleton interval, the root counts only if both ends are included.
    """
    left=None if left is None else _q(left)
    right=None if right is None else _q(right)
    return _count(sturm_sequence(coefficients),left,right,include_left,include_right)


def compare_largest_real_root(coefficients,bound):
    """Return -1, 0 or +1 according as the largest real root is <, = or > bound."""
    bound=_q(bound)
    sequence=sturm_sequence(coefficients)
    if not _count(sequence,None,None):
        raise ValueError("polynomial has no real roots")
    if _count(sequence,bound,None):
        return 1
    if not polynomial_value(sequence[0],bound):
        return 0
    return -1


def isolate_largest_real_root(coefficients,bits=80):
    """Return rational endpoints of width <= 2**(-bits) containing only the
    largest distinct real root. If the root is hit exactly, return (root,root).
    """
    if not isinstance(bits,int) or bits<0:
        raise ValueError("bits must be a nonnegative integer")
    sequence=sturm_sequence(coefficients)
    p=sequence[0]
    if not _count(sequence,None,None):
        raise ValueError("polynomial has no real roots")
    # Strict Cauchy bound, leaving both endpoints outside the full root set.
    radius=2+max((abs(a/p[0]) for a in p[1:]),default=Fraction(0))
    left,right=-radius,radius
    tolerance=Fraction(1,1<<bits)
    while right-left>tolerance or _count(sequence,left,right,include_left=True)!=1:
        middle=(left+right)/2
        if _count(sequence,middle,None):
            left=middle
        elif not polynomial_value(p,middle):
            return middle,middle
        else:
            right=middle
    return left,right


def compare_perron_to_rational(matrix,bound):
    """Exact spectral-radius comparison for a nonempty NONNEGATIVE matrix.

    Perron--Frobenius guarantees that its spectral radius is a real root and
    the largest real root of its characteristic polynomial. No irreducibility
    or pseudo-Anosov conclusion is inferred here.
    """
    n,m=_shape(matrix)
    if not n or n!=m or any(a<0 for row in matrix for a in row):
        raise ValueError("Perron comparison requires a nonempty nonnegative square matrix")
    return compare_largest_real_root(charpoly(matrix),bound)


def compare_largest_real_roots(left_coefficients,right_coefficients):
    """Compare two largest real algebraic roots exactly, including equality.

    Isolating intervals eventually separate unequal algebraic numbers. For
    equality, the polynomial gcd has a root in both isolating intervals;
    each interval contains exactly one root of its original polynomial.
    """
    left_poly=sturm_sequence(left_coefficients)[0]
    right_poly=sturm_sequence(right_coefficients)[0]
    common=_gcd(left_poly,right_poly)
    bits=0
    while True:
        a,b=isolate_largest_real_root(left_poly,bits)
        c,d=isolate_largest_real_root(right_poly,bits)
        if b<c:
            return -1
        if d<a:
            return 1
        if len(common)>1 and count_real_roots(common,max(a,c),min(b,d),include_left=True):
            return 0
        bits=4 if bits==0 else 2*bits
