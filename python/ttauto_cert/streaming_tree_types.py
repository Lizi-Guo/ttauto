"""Lazy version of the planted multigon-tree grammar.

Only the recursion stack is kept. Canonical codes can repeat for different
univalent roots; callers must quotient these (on disk for large searches).
The grammar and rooted count are the same as tree_types.py, independent of
any growth cutoff. No truncated stream is an exhaustive state catalogue.
"""
from itertools import product
from math import comb, factorial, prod

from .tree_types import _validate_stratum, canonical_code, positive_compositions


def rooted_count(n, interior_prongs):
    prongs = _validate_stratum(n, interior_prongs)
    colors = (1,) + tuple(sorted(set(prongs)))
    counts = (n-1,) + tuple(prongs.count(p) for p in colors[1:])
    total = sum(counts)
    defect = sum(c*(2-p) for c, p in zip(counts, colors))
    numerator = factorial(total-1)*comb(3*total-2, defect-1)
    answer, remainder = divmod(numerator, prod(factorial(c) for c in counts))
    assert remainder == 0
    return answer


def iter_rooted_codes(n, interior_prongs):
    """Yield the exact canonical code of each rooted decorated tree lazily."""
    prongs = _validate_stratum(n, interior_prongs)
    colors = (1,) + tuple(sorted(set(prongs)))
    inventory = (n-1,) + tuple(prongs.count(p) for p in colors[1:])
    # This explicit implementation uses nested Python generators. Reject
    # unsupported depth instead of silently returning a partial catalogue.
    if sum(inventory) > 100:
        raise ValueError('streaming grammar supports at most 101 multigons')

    def defect(counts):
        return sum(c*(2-p) for c, p in zip(counts, colors))

    def forests(counts, minimum):
        if not any(counts):
            if minimum == 0:
                yield ()
            return
        if defect(counts) < max(1, minimum):
            return
        # Larger first components make a first tree available immediately;
        # no smaller inventory is materialized as a table first.
        for first in product(*(range(c, -1, -1) for c in counts)):
            if not any(first) or defect(first) < 1:
                continue
            rest = tuple(a-b for a, b in zip(counts, first))
            if any(rest) and defect(rest) < max(1, minimum-1):
                continue
            if not any(rest) and minimum > 1:
                continue
            for tree in trees(first):
                for tail in forests(rest, max(0, minimum-1)):
                    yield (tree,) + tail

    def trees(counts):
        if not any(counts) or defect(counts) < 1:
            return
        # High-prong roots first: useful for low-dilatation strata, while
        # preserving every root color and every ordered forest.
        for color in reversed(range(len(colors))):
            p = colors[color]
            if counts[color] == 0:
                continue
            rest = list(counts)
            rest[color] -= 1
            if sum(rest) < p-1:
                continue
            for children in forests(tuple(rest), p-1):
                length = len(children)
                for before in range(length-p+2):
                    if p == 1:
                        yield (p, children[:before], children[before:], ())
                        continue
                    for after in range(length-before-p+2):
                        for parts in positive_compositions(length-before-after, p-1):
                            offset = before+after
                            other = []
                            for size in parts:
                                other.append(children[offset:offset+size])
                                offset += size
                            yield (p, children[:before], children[before:before+after], tuple(other))

    def to_orders(tree):
        orders = [((1,),), None]
        stack = [(tree, 1, 0)]
        while stack:
            node, here, parent = stack.pop()
            p, before, after, other = node
            groups = []
            for group in (before+(None,)+after,) + other:
                neighbors = []
                for child in group:
                    if child is None:
                        neighbors.append(parent)
                    else:
                        index = len(orders)
                        orders.append(None)
                        neighbors.append(index)
                        stack.append((child, index, here))
                groups.append(tuple(neighbors))
            assert len(groups) == p
            orders[here] = tuple(groups)
        return tuple(orders)

    for tree in trees(inventory):
        yield canonical_code(to_orders(tree))


def unique_codes_on_disk(connection, n, prongs, statistics):
    """Quotient root duplicates in SQLite, rebuilding the stream on resume."""
    connection.execute('CREATE TEMP TABLE stream_seen (code TEXT PRIMARY KEY)')
    statistics.update(algorithm='lazy_planted_tree', rooted_generated=0,
                      types=0, rooted_count_formula=rooted_count(n, prongs),
                      completed=False)
    try:
        for code in iter_rooted_codes(n, prongs):
            statistics['rooted_generated'] += 1
            key = ' '.join(map(str, code))
            cursor = connection.execute('INSERT OR IGNORE INTO stream_seen VALUES (?)', (key,))
            if cursor.rowcount:
                statistics['types'] += 1
                yield code
        if statistics['rooted_generated'] != statistics['rooted_count_formula']:
            raise ArithmeticError('streamed tree count disagrees with Lagrange inversion')
        statistics['completed'] = True
    finally:
        connection.execute('DROP TABLE stream_seen')
