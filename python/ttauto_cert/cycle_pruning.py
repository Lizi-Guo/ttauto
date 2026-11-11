"""Exact cycle-growth bounds for entrywise nonnegative matrix completions.

This transfers the draft's lower curve-graph method to partially assigned
edge words. No floral geometry, known minimum, reciprocity or determinant
assumption is needed. If B is the partial matrix and A >= B is a completion,
rho(A) >= rho(B) >= growth(H) for any selected actual cycle graph H of B.
Only a STRICT negative sign at 1/U rejects the inclusive search rho(A) <= U.
Cycle-discovery limits weaken the test; they never establish exhaustion.
"""
from fractions import Fraction
from functools import lru_cache
from math import prod


def _matrix(matrix):
    matrix = tuple(tuple(row) for row in matrix)
    n = len(matrix)
    if not n or any(len(row) != n for row in matrix):
        raise ValueError('a nonempty square matrix is required')
    if any(type(x) is not int or x < 0 for row in matrix for x in row):
        raise ValueError('matrix entries must be nonnegative integers')
    return matrix


def clique_coefficients(cycles):
    """Ascending Q(t); one route represents all its parallel-edge choices."""
    cliques = [(0, 0, 1)]
    polynomial = {0: 1}
    for cycle in cycles:
        vertices = cycle['vertices']
        mask = sum(1 << v for v in vertices)
        added = [(used | mask, degree+len(vertices), -coefficient*cycle['multiplicity'])
                 for used, degree, coefficient in cliques if not used & mask]
        for _, degree, coefficient in added:
            polynomial[degree] = polynomial.get(degree, 0)+coefficient
        cliques.extend(added)
    return [polynomial.get(k, 0) for k in range(max(polynomial)+1)]


def _integer_value(polynomial, upper):
    a, b = upper.numerator, upper.denominator
    d = len(polynomial)-1
    return sum(c*b**k*a**(d-k) for k, c in enumerate(polynomial))


class CycleGrowthPruner:
    """Bounded cycle discovery, exact integer signs, and a bounded LRU cache."""
    def __init__(self, upper, *, max_length=4, max_cycles=64, max_steps=4096,
                 clique_limit=12):
        if isinstance(upper, float):
            raise TypeError('use an exact rational upper bound')
        self.upper = Fraction(upper)
        if self.upper <= 1:
            raise ValueError('upper must exceed one')
        for value in (max_length, max_cycles, max_steps, clique_limit):
            if type(value) is not int or value < 1:
                raise ValueError('positive integer discovery limits required')
        self.max_length, self.max_cycles = max_length, max_cycles
        self.max_steps, self.clique_limit = max_steps, clique_limit
        self._cached = lru_cache(maxsize=8192)(self._find)

    def witness(self, matrix):
        return self._cached(_matrix(matrix))

    def cache_info(self):
        return self._cached.cache_info()._asdict()

    def excluded_leading_columns(self, columns, assigned):
        """Trusted count-column fast path used by the edge-word enumerator.

        Return a principal-block index tuple with a growth obstruction, or
        None. The caller constructed nonnegative integer counts. Up to size
        three, negative principal det(U I-B) gives the same exact clique sign
        without building cycles or allocating a cache key at every word.
        A missed sign is harmless; equality is always retained.
        """
        a, b = self.upper.numerator, self.upper.denominator
        diagonal = [a-b*columns[i][i] for i in range(assigned)]
        for i, value in enumerate(diagonal):
            if value < 0:
                return (i,)
        if assigned <= 3:
            for i in range(assigned):
                for j in range(i):
                    if diagonal[i]*diagonal[j]-b*b*columns[i][j]*columns[j][i] < 0:
                        return (j, i)
            if assigned == 3:
                x, y, z = diagonal
                c = columns
                determinant = (x*y*z-b*b*(x*c[1][2]*c[2][1]+y*c[0][2]*c[2][0]+z*c[0][1]*c[1][0])
                               -b**3*(c[0][1]*c[1][2]*c[2][0]+c[0][2]*c[2][1]*c[1][0]))
                if determinant < 0:
                    return (0, 1, 2)
            return None
        matrix = tuple(tuple(columns[j][i] for j in range(assigned)) for i in range(assigned))
        return tuple(range(assigned)) if self._cached(matrix) is not None else None

    def _find(self, matrix):
        n = len(matrix)
        limit = min(n, self.max_length)
        a, b = self.upper.numerator, self.upper.denominator
        weights = [b**k*a**(limit-k) for k in range(limit+1)]
        anchor_sums = [0]*n
        adjacency = [[v for v in range(n) if matrix[v][u]] for u in range(n)]
        cycles, steps = [], 0
        for start in range(n):
            # A route is generated only at its minimum vertex, so rotations
            # are not counted twice. Opposite directions remain distinct.
            stack = [(start, (start,), 1, 1 << start)]
            while stack and steps < self.max_steps and len(cycles) < self.max_cycles:
                current, path, multiplicity, used = stack.pop()
                for following in adjacency[current]:
                    steps += 1
                    if steps > self.max_steps:
                        break
                    mult = multiplicity*matrix[following][current]
                    if following == start:
                        cycle = {'vertices': list(path), 'multiplicity': mult}
                        cycles.append(cycle)
                        for vertex in path:
                            anchor_sums[vertex] += mult*weights[len(path)]
                            if anchor_sums[vertex] > weights[0]:
                                chosen = [c for c in cycles if vertex in c['vertices']]
                                return self._record(matrix, chosen, 'common_vertex_cycles', vertex)
                        if len(cycles) >= self.max_cycles:
                            break
                    elif following > start and not used & (1 << following) and len(path) < limit:
                        stack.append((following, path+(following,), mult, used | (1 << following)))
            if steps >= self.max_steps or len(cycles) >= self.max_cycles:
                break
        chosen = sorted(cycles, key=lambda c: (len(c['vertices']), c['vertices']))[:self.clique_limit]
        if chosen and _integer_value(clique_coefficients(chosen), self.upper) < 0:
            return self._record(matrix, chosen, 'induced_cycle_clique', None)
        return None

    def _record(self, matrix, cycles, rule, anchor):
        polynomial = clique_coefficients(cycles)
        value = _integer_value(polynomial, self.upper)
        if value >= 0:
            raise ArithmeticError('a non-strict sign cannot prune an inclusive upper bound')
        return {'rule': rule, 'anchor': anchor, 'partial_matrix': [list(row) for row in matrix],
                'cycles': cycles, 'clique_coefficients_ascending': polynomial,
                'upper': str(self.upper), 'negative_integer_value': value,
                'completion_requirement': 'entrywise_nonnegative_extension',
                'floral_assumption_used': False}


def verify_cycle_witness(witness):
    """Independent subset/Fraction replay, including actual directed routes."""
    matrix = _matrix(witness['partial_matrix'])
    upper = Fraction(witness['upper'])
    if upper <= 1 or witness['completion_requirement'] != 'entrywise_nonnegative_extension':
        raise ValueError('invalid witness scope')
    cycles = witness['cycles']
    seen = set()
    for cycle in cycles:
        path = tuple(cycle['vertices'])
        if (not path or len(set(path)) != len(path) or
                any(type(v) is not int or not 0 <= v < len(matrix) for v in path)):
            raise ValueError('cycle must be embedded')
        canonical = min(path[k:]+path[:k] for k in range(len(path)))
        if canonical in seen:
            raise ValueError('duplicate cycle rotation')
        seen.add(canonical)
        actual = prod(matrix[path[(j+1) % len(path)]][v] for j, v in enumerate(path))
        if type(cycle['multiplicity']) is not int or not 0 < cycle['multiplicity'] <= actual:
            raise ValueError('cycle multiplicity is not supplied by matrix edges')
    # A common-anchor certificate has no cliques of size >= 2. This is
    # checked directly to avoid enumerating 2**64 subsets of an exit family.
    if witness['rule'] == 'common_vertex_cycles':
        if not cycles or not all(witness['anchor'] in c['vertices'] for c in cycles):
            raise ValueError('cycles do not share the asserted anchor')
        coefficients = [1]+[0]*max(len(c['vertices']) for c in cycles)
        for cycle in cycles:
            coefficients[len(cycle['vertices'])] -= cycle['multiplicity']
    elif witness['rule'] == 'induced_cycle_clique':
        coefficients = [0]*(len(matrix)+1)
        for mask in range(1 << len(cycles)):
            used, degree, coefficient = set(), 0, 1
            for j, cycle in enumerate(cycles):
                if mask >> j & 1:
                    support = set(cycle['vertices'])
                    if used & support:
                        break
                    used |= support
                    degree += len(support)
                    coefficient *= -cycle['multiplicity']
            else:
                coefficients[degree] += coefficient
        while len(coefficients) > 1 and coefficients[-1] == 0:
            coefficients.pop()
    else:
        raise ValueError('unknown cycle rule')
    expected = list(witness['clique_coefficients_ascending'])
    while len(expected) > 1 and expected[-1] == 0:
        expected.pop()
    if coefficients != expected:
        raise ValueError('incorrect clique polynomial')
    if sum(Fraction(c)*upper**(-k) for k, c in enumerate(coefficients)) >= 0:
        raise ValueError('no strict negative rational witness')
    if _integer_value(witness['clique_coefficients_ascending'], upper) != witness['negative_integer_value']:
        raise ValueError('incorrect integer sign witness')
    return True
