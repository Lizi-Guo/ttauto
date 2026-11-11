"""Finite standard-map enumeration, independent of fold reachability.

For a supplied genus-zero standard multigon tree with m real edges, every
real edge image is R(IR)*. The Ham--Song bound sum(M) <= floor(U**m)+m-1
therefore bounds complete edge words as well as real counts. Hironaka--Tsang
provides such a map for f itself after puncturing its interior singularities
and choosing its fixed outer puncture as the sole outer puncture.

This module exhausts ONE supplied state. Callers must independently enumerate
all state types and classify the emitted surface maps; a surface automorphism
with a primitive real matrix is not automatically pseudo-Anosov. Explicit
budgets and early iterator closure never report completion.
"""
from collections import Counter, deque
from fractions import Fraction
from functools import lru_cache
from itertools import permutations
from math import gcd

from .exact_math import charpoly, compare_largest_real_root
from .cycle_pruning import CycleGrowthPruner
from .ribbon_verify import (cyclic_word, inverse, reduce_word, substitute,
                            validate_graph)


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def lazy_product(factories):
    """Cartesian product without pooling factorial-size input iterators.

    Each zero-argument factory must return a fresh iterator. The explicit
    stack also avoids a call-depth limit on the number of factors.
    """
    factories = tuple(factories)
    if not factories:
        yield ()
        return
    values = []
    stack = [iter(factories[0]())]
    while stack:
        depth = len(stack)-1
        try:
            value = next(stack[-1])
        except StopIteration:
            stack.pop()
            del values[depth:]
            continue
        del values[depth:]
        values.append(value)
        if len(stack) == len(factories):
            yield tuple(values)
        else:
            stack.append(iter(factories[len(stack)]()))


def exact_sum_bound(upper, dimension):
    if isinstance(upper, float):
        raise TypeError('upper bound must be an exact rational, not float')
    upper = Fraction(upper)
    _require(upper > 1, 'spectral upper bound must exceed one')
    _require(type(dimension) is int and dimension > 0, 'positive matrix dimension required')
    power = upper**dimension
    return power.numerator//power.denominator+dimension-1


def primitive(matrix):
    """Strong connectivity and directed cycle gcd, using integer arithmetic."""
    n = len(matrix)
    _require(n > 0 and all(len(row) == n for row in matrix), 'nonempty square matrix required')
    _require(all(type(x) is int and x >= 0 for row in matrix for x in row),
             'nonnegative integer matrix required')
    adjacency = [[i for i in range(n) if matrix[i][j]] for j in range(n)]
    reverse = [[j for j in range(n) if matrix[i][j]] for i in range(n)]

    def distances(graph):
        result = {0: 0}
        queue = deque([0])
        while queue:
            v = queue.popleft()
            for w in graph[v]:
                if w not in result:
                    result[w] = result[v]+1
                    queue.append(w)
        return result
    depth = distances(adjacency)
    if len(depth) != n or len(distances(reverse)) != n:
        return False
    period = 0
    for v, row in enumerate(adjacency):
        for w in row:
            period = gcd(period, depth[v]+1-depth[w])
    return period == 1


def surjective_free_group(words, rank):
    """Stallings folds of the image subgroup, not just its abelianization."""
    _require(type(rank) is int and rank > 0 and len(words) == rank, 'invalid free-group rank')
    _require(all(type(d) is int and 0 < abs(d) <= rank for w in words for d in w),
             'free-group letter outside the alphabet')
    parent, edges = [0], []

    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v

    for word in words:
        reduced = reduce_word(word)
        current = 0
        for i, letter in enumerate(reduced):
            following = 0
            if i+1 != len(reduced):
                following = len(parent)
                parent.append(following)
            edges.extend([(current, letter, following), (following, -letter, current)])
            current = following
    changed = True
    while changed:
        changed = False
        outgoing = {}
        for u, label, v in edges:
            u, v = find(u), find(v)
            key = (u, label)
            if key in outgoing:
                other = find(outgoing[key])
                if other != v:
                    parent[v] = other
                    changed = True
            else:
                outgoing[key] = v
    base = find(0)
    return ({find(i) for i in range(len(parent))} == {base} and
            {(find(u), label, find(v)) for u, label, v in edges} ==
            {(base, sign*j, base) for j in range(1, rank+1) for sign in (-1, 1)})


class StandardState:
    """Validated colored standard tree and an explicit graph/free-group basis."""
    def __init__(self, graph):
        self.graph = graph
        self.origins, boundaries = validate_graph(graph)
        self.m = graph['main_edges']
        self.q = len(graph['polygons'])
        self.e = len(graph['endpoints'])
        _require(self.m == self.q-1 and self.m > 0, 'require one-outer-region sphere tree')
        polygon_darts = [d for cycle in graph['polygons'] for d in cycle]
        _require(sorted(map(abs, polygon_darts)) == list(range(self.m+1, self.e+1)),
                 'infinitesimal polygons must partition the infinitesimal edges')
        for star in graph['cyclic_directions']:
            kinds = [abs(d) > self.m for d in star]
            _require(sum(kinds) == 2 and len(kinds) > 2,
                     'each standard switch needs two infinitesimal and positive real incidence')
            _require(sum(kinds[i] != kinds[i-1] for i in range(len(kinds))) == 2,
                     'real and infinitesimal half-edges must form the two smoothing sides')
        self.omitted = {abs(p[-1]): (i+1)*(1 if p[-1] > 0 else -1)
                        for i, p in enumerate(graph['polygons'])}
        # The real tree together with all but one side of each polygon is a
        # spanning tree, including for a monogon. Verify this combinatorially.
        paths = {0: []}
        queue = deque([0])
        tree_edges = set()
        while queue:
            v = queue.popleft()
            for d in graph['cyclic_directions'][v]:
                if abs(d) in self.omitted:
                    continue
                w = self.origins[-d]
                tree_edges.add(abs(d))
                if w not in paths:
                    paths[w] = paths[v]+[d]
                    queue.append(w)
        _require(len(paths) == len(graph['cyclic_directions']) and
                 len(tree_edges) == len(paths)-1, 'complement of omitted polygon edges is not a tree')
        self.basis = [paths[self.origins[p[-1]]]+[p[-1]]+
                      inverse(paths[self.origins[-p[-1]]]) for p in graph['polygons']]
        for i, word in enumerate(self.basis, 1):
            _require(self.project(word) == [i], 'incorrect graph/free-group basis')
        self.outer = next((w for w in boundaries if any(abs(d) <= self.m for d in w)), None)
        _require(self.outer is not None, 'missing outer boundary')
        self.outer_class = cyclic_word(self.project(self.outer))

    def project(self, word):
        return reduce_word([(1 if d > 0 else -1)*self.omitted[abs(d)]
                            for d in word if abs(d) in self.omitted])

    def peripheral_transports(self):
        """Every color-preserving polygon permutation and oriented rotation."""
        polygons = self.graph['polygons']
        groups = {}
        for i, polygon in enumerate(polygons):
            key = (len(polygon), self.graph['punctured'][i], self.graph['labels'][i])
            groups.setdefault(key, []).append(i)
        blocks = list(groups.values())
        for block_images in lazy_product(lambda block=block: permutations(block) for block in blocks):
            polygon_images = [None]*self.q
            for block, targets in zip(blocks, block_images):
                for a, b in zip(block, targets):
                    polygon_images[a] = b
            for rotations in lazy_product(lambda p=p: iter(range(len(p))) for p in polygons):
                vertices = [None]*len(self.graph['cyclic_directions'])
                words = [None]*self.e
                for p, target in enumerate(polygon_images):
                    for j, d in enumerate(polygons[p]):
                        image = polygons[target][(j+rotations[p]) % len(polygons[p])]
                        vertices[self.origins[d]] = self.origins[image]
                        words[abs(d)-1] = [image if d > 0 else -image]
                _require(sorted(vertices) == list(range(len(vertices))), 'polygon transport is not bijective')
                yield vertices, words, polygon_images

    def surface_automorphism(self, words, polygon_images):
        images = [self.project(substitute(words, word)) for word in self.basis]
        if not surjective_free_group(images, self.q):
            return None, 'not_automorphism'
        for i, word in enumerate(images):
            if cyclic_word(word) != (polygon_images[i]+1,):
                return None, 'peripheral_mismatch'
        if cyclic_word(self.project(substitute(words, self.outer))) != self.outer_class:
            return None, 'outer_mismatch'
        # Artin's characterization applied after conjugating the image of the
        # outer word to the outer word. This is an outer surface automorphism.
        return images, None


class _BudgetExhausted(Exception):
    pass


class StandardMapSearch:
    """Stream bounded surface candidates from one state; completion is explicit.

    max_maps bounds endpoint-compatible full maps; max_prefixes bounds real
    word extensions, including extensions which never complete a map. Neither
    cutoff is a mathematical exhaustion certificate. This object is single-use.
    """
    def __init__(self, graph, upper, *, max_maps=None, max_prefixes=None, cycle_pruning=False):
        self.state = StandardState(graph)
        self.bound = exact_sum_bound(upper, self.state.m)
        self.upper = Fraction(upper)
        _require(type(cycle_pruning) is bool, 'cycle_pruning must be boolean')
        self.cycle_pruner = CycleGrowthPruner(self.upper) if cycle_pruning else None
        self.cycle_witnesses = []
        for limit in (max_maps, max_prefixes):
            _require(limit is None or type(limit) is int and limit >= 0, 'nonnegative budget required')
        self.max_maps, self.max_prefixes = max_maps, max_prefixes
        self.counts = Counter()
        self.status = 'not_started'

    def _increment(self, key, limit=None):
        if limit is not None and self.counts[key] >= limit:
            raise _BudgetExhausted(key)
        self.counts[key] += 1

    def _words(self, start, finish, length, column_filter=None):
        s = self.state
        _require(type(length) is int and length >= 1, 'positive real word length required')

        def extensions(v, remaining):
            for real in s.graph['cyclic_directions'][v]:
                if abs(real) > s.m:
                    continue
                self._increment('word_prefixes', self.max_prefixes)
                reached = s.origins[-real]
                if remaining == 1:
                    yield reached, (real,)
                    continue
                for inf in s.graph['cyclic_directions'][reached]:
                    if abs(inf) > s.m:
                        yield s.origins[-inf], (real, inf)

        # Each frame remembers the prefix length on entry. Reusing one word
        # buffer avoids quadratic copying while descending a long path.
        prefix = []
        column = [0]*s.m if column_filter is not None else None

        def truncate(base):
            if column is not None:
                for dart in prefix[base:]:
                    if abs(dart) <= s.m:
                        column[abs(dart)-1] -= 1
            del prefix[base:]

        stack = [(iter(extensions(start, length)), length, 0)]
        while stack:
            choices, remaining, base = stack[-1]
            try:
                reached, segment = next(choices)
            except StopIteration:
                stack.pop()
                truncate(base)
                continue
            truncate(base)
            prefix.extend(segment)
            if column is not None:
                column[abs(segment[0])-1] += 1
                if not column_filter(column):
                    continue
            if remaining == 1:
                if reached == finish:
                    yield list(prefix)
            else:
                stack.append((iter(extensions(reached, remaining-1)),
                              remaining-1, len(prefix)))

    def _assign(self, vertices, words, total):
        """All endpoint-compatible full maps with this real-letter total."""
        s = self.state
        columns = [None]*s.m
        pruner = getattr(self, 'cycle_pruner', None)

        def allowed_column(column, edge, lookup):
            columns[edge] = column
            self.counts['cycle_growth_checks'] += 1
            support = lookup(tuple(column))
            if support is None:
                return True
            self.counts['cycle_growth_prunes'] += 1
            self.counts['cycle_prunes_at_assigned_columns_'+str(edge+1)] += 1
            if len(self.cycle_witnesses) < 8:
                partial = tuple(tuple(columns[j][i] for j in support) for i in support)
                witness = pruner.witness(partial)
                if witness is None:
                    raise ArithmeticError('fast principal sign lacks a cycle witness')
                if witness not in self.cycle_witnesses:
                    self.cycle_witnesses.append(witness)
            return False

        def choices(edge, remaining):
            a, b = s.graph['endpoints'][edge]
            # Earlier columns stay fixed throughout this generator's life.
            # Many signed words have the same current count column.
            @lru_cache(maxsize=512)
            def lookup(column):
                self.counts['cycle_column_cache_misses'] += 1
                return pruner.excluded_leading_columns(columns, edge+1)
            check = (lambda column: allowed_column(column, edge, lookup)) if pruner is not None else None
            for length in range(1, remaining-(s.m-edge-1)+1):
                iterator = (self._words(vertices[a], vertices[b], length, column_filter=check)
                            if check is not None else self._words(vertices[a], vertices[b], length))
                for word in iterator:
                    yield word, length

        stack = [(iter(choices(0, total)), total)]
        while stack:
            edge = len(stack)-1
            options, remaining = stack[-1]
            try:
                word, length = next(options)
            except StopIteration:
                stack.pop()
                continue
            words[edge] = word
            if edge+1 == s.m:
                if remaining == length:
                    yield [list(w) for w in words]
            else:
                stack.append((iter(choices(edge+1, remaining-length)), remaining-length))

    def run(self):
        _require(self.status == 'not_started', 'search object is single-use')
        self.status = 'running'
        s = self.state

        @lru_cache(maxsize=8192)
        def allowed(flat):
            polynomial = charpoly([list(flat[i*s.m:(i+1)*s.m]) for i in range(s.m)])
            valid = (compare_largest_real_root(polynomial, Fraction(1)) > 0 and
                     compare_largest_real_root(polynomial, self.upper) <= 0)
            return valid, polynomial

        try:
            for vertices, infinitesimal_words, polygon_images in s.peripheral_transports():
                for total in range(s.m, self.bound+1):
                    for words in self._assign(vertices, infinitesimal_words, total):
                        self._increment('endpoint_compatible_maps', self.max_maps)
                        matrix = [[sum(abs(d) == i for d in word) for word in words[:s.m]]
                                  for i in range(1, s.m+1)]
                        if not primitive(matrix):
                            continue
                        self.counts['primitive_maps'] += 1
                        valid, polynomial = allowed(tuple(x for row in matrix for x in row))
                        if not valid:
                            continue
                        self.counts['primitive_below_upper'] += 1
                        images, rejection = s.surface_automorphism(words, polygon_images)
                        if rejection == 'not_automorphism':
                            continue
                        self.counts['free_group_automorphisms'] += 1
                        if rejection:
                            self.counts[rejection] += 1
                            continue
                        self.counts['oriented_surface_automorphisms'] += 1
                        yield {'vertex_images': list(vertices), 'edge_images': words,
                               'polygon_images': list(polygon_images), 'matrix': matrix,
                               'charpoly': polynomial, 'free_group_images': images,
                               'pseudo_anosov_certified': False,
                               'surface_map_recognition': 'Stallings+oriented peripheral words'}
            self.status = 'completed'
        except _BudgetExhausted as error:
            self.status = 'truncated'
            self.limit_reason = str(error)
        except GeneratorExit:
            self.status = 'aborted'
            raise
        except BaseException:
            self.status = 'failed'
            raise

    def summary(self):
        return {'schema': 'bounded-standard-state-search-v1', 'status': self.status,
                'completed': self.status == 'completed', 'upper_bound': str(self.upper),
                'real_dimension': self.state.m, 'real_sum_bound': self.bound,
                'full_real_image_length_bound': 2*self.bound-self.state.m,
                'max_maps': self.max_maps, 'max_prefixes': self.max_prefixes,
                'cycle_growth_pruning': self.cycle_pruner is not None,
                'cycle_pruning_cache': self.cycle_pruner.cache_info() if self.cycle_pruner else None,
                'cycle_pruning_witness_samples': self.cycle_witnesses,
                'counts': dict(self.counts), 'uses_fold_reachability': False,
                'all_stratum_states_certified': False, 'pseudo_anosov_certified': False}
