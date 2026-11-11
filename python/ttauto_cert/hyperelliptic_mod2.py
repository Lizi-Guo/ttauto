"""Exact necessary homology conditions from hyperelliptic branch permutations.

For a double cover of the sphere with N=2g+2 branch points, H1(F2) is the
even-weight permutation module modulo its all-ones vector. Hence the branch
permutation characteristic polynomial is (x+1)^2 times the homology polynomial.
This is a necessary condition, never a realizability or pseudo-Anosov test.
All public polynomial coefficient lists are in descending order.
"""
from itertools import product

SOURCE = 'https://londmathsoc.onlinelibrary.wiley.com/doi/full/10.1112/plms.12565'


def _bits(coefficients):
    if not coefficients or any(type(value) is not int for value in coefficients):
        raise ValueError('a nonempty exact integer polynomial is required')
    result = 0
    for value in coefficients:
        result = (result << 1) | (value % 2)
    return result


def _coefficients(bits):
    return [(bits >> i) & 1 for i in reversed(range(bits.bit_length()))] or [0]


def _multiply(left, right):
    result = 0
    while right:
        if right & 1:
            result ^= left
        left <<= 1
        right >>= 1
    return result


def _divide(dividend, divisor):
    if divisor == 0:
        raise ZeroDivisionError('zero polynomial modulo 2')
    quotient = 0
    while dividend and dividend.bit_length() >= divisor.bit_length():
        shift = dividend.bit_length() - divisor.bit_length()
        quotient ^= 1 << shift
        dividend ^= divisor << shift
    return quotient, dividend


def polynomial_divides_mod2(divisor, dividend):
    """Whether divisor divides dividend in F2[x]; a zero divisor is invalid."""
    return _divide(_bits(dividend), _bits(divisor))[1] == 0


def branch_permutation_charpoly(cycle_lengths):
    """Return product(x^length+1)/(x+1)^2 in F2[x], checking exact division.

    The cycle lengths describe all branch points, including fixed branch points.
    An even total N>=4 is required; the returned polynomial has degree N-2.
    """
    cycle_lengths = tuple(cycle_lengths)
    if any(type(length) is not int or length < 1 for length in cycle_lengths):
        raise ValueError('branch cycle lengths must be positive integers')
    size = sum(cycle_lengths)
    if size < 4 or size % 2:
        raise ValueError('the branch count must be even and at least four')
    full = 1
    for length in cycle_lengths:
        full = _multiply(full, (1 << length) | 1)
    quotient, remainder = _divide(full, 5)  # (x+1)^2=x^2+1 in F2.
    if remainder or quotient.bit_length() != size - 1:
        raise ArithmeticError('branch characteristic polynomial quotient failed')
    return _coefficients(quotient)


def _partitions(total, first=1):
    if total == 0:
        yield ()
    for size in range(first, total + 1):
        for tail in _partitions(total - size, size):
            yield (size,) + tail


def branch_cycle_types(groups, fixed_points=0):
    """Enumerate independent permutation groups and distinguished fixed points.

    groups=(6,3), fixed_points=1 describes six points permuted among themselves,
    three further points permuted among themselves, and one fixed point. Colors
    are specified by the group positions; different groups never interchange.
    """
    groups = tuple(groups)
    if any(type(size) is not int or size < 1 for size in groups):
        raise ValueError('branch group sizes must be positive integers')
    if type(fixed_points) is not int or fixed_points < 0:
        raise ValueError('the fixed branch-point count must be nonnegative')
    total = sum(groups) + fixed_points
    if total < 4 or total % 2:
        raise ValueError('the branch count must be even and at least four')
    for selected in product(*[tuple(_partitions(size)) for size in groups]):
        cycles = tuple(length for group in selected for length in group) + (1,) * fixed_points
        yield {'group_cycles': [list(group) for group in selected],
               'fixed_points': fixed_points, 'branch_cycles': list(cycles),
               'mod2_polynomial': branch_permutation_charpoly(cycles),
               'exact_division_verified': True}


def compatible_branch_cycle_types(polynomial, groups, fixed_points=0):
    """Return every permitted branch cycle type matching polynomial modulo 2.

    An empty list is an exact obstruction for this branch configuration. A
    nonempty list is only a necessary-condition witness, not a surface map.
    """
    groups = tuple(groups)
    allowed = list(branch_cycle_types(groups, fixed_points))
    residue = _bits(polynomial)
    degree = sum(groups) + fixed_points - 2
    if len(polynomial) != degree + 1 or polynomial[0] != 1:
        raise ValueError('a monic polynomial of degree branch_count-2 is required')
    result = []
    for cycle_type in allowed:
        if residue != _bits(cycle_type['mod2_polynomial']):
            continue
        # Check the equivalent full permutation identity as well as quotient
        # equality, including multiplicities in characteristic two.
        full = 1
        for length in cycle_type['branch_cycles']:
            full = _multiply(full, (1 << length) | 1)
        if _multiply(residue, 5) != full:
            raise ArithmeticError('inconsistent exact branch polynomial identity')
        result.append(cycle_type)
    return result
