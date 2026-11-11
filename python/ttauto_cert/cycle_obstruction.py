# SPDX-License-Identifier: GPL-3.0-or-later
"""Abstract cycle-graph exclusion, migrated from the large-disc proof jobs.

Enumerates a containing family of cycle disjointness graphs of a putative
primitive matrix with a prescribed characteristic polynomial. A survivor is
NOT a matrix, a ribbon graph, or a surface map. No logarithmic monotonicity
assumption is used. Budgets and connected survivors never certify exclusion.
"""
from hashlib import sha256
import json
from .lefschetz import traces, partitions
from .spectral import ordinary_orbits, multiply
from .cover_orbits import ModelBudgetExceeded


def signed_peripheral_polynomials(case, budget):
    from collections import Counter
    from itertools import product
    groups = list(Counter(p for p in case['interior_prongs'] if p % 2 == 0).values())
    if case['outer_prongs'] % 2 == 0:
        groups.append(1)  # distinguished outer face cannot be exchanged
    def group_factors(number):
        for lengths in partitions(number):
            for signs in product((-1, 1), repeat=len(lengths)):
                budget.tick()
                result = [1]
                for length, sign in zip(lengths, signs):
                    result = multiply(result, [1]+[0]*(length-1)+[-sign])
                yield result
    def visit(i, polynomial):
        if i == len(groups):
            budget.tick()
            yield polynomial
            return
        for factor in group_factors(groups[i]):
            yield from visit(i+1, multiply(polynomial, factor))
    yield from visit(0, [1])


def cycle_search(target, max_nodes=200000):
    n = len(target)-1
    caps, failure = ordinary_orbits(traces(target, n)[1:])
    if failure:
        return {'completed': True, 'excluded': True, 'reason': 'negative_primitive_cycle_count',
                'failure': failure, 'nodes': 0, 'target': target}
    state = {'nodes': 0, 'capacity_prunes': 0, 'matching_graphs': 0}
    transcript = sha256()
    found = None
    def connected(adj):
        seen, pending = {0}, [0]
        while pending:
            i = pending.pop()
            for j in range(len(adj)):
                if j not in seen and not (adj[i] >> j) & 1:
                    seen.add(j); pending.append(j)
        return len(seen) == len(adj)
    def visit(lengths, adj, cliques, p):
        nonlocal found
        if max_nodes is not None and state['nodes'] >= max_nodes:
            raise ModelBudgetExceeded('cycle search incomplete')
        state['nodes'] += 1
        transcript.update(json.dumps([lengths, adj, p], separators=(',', ':')).encode()+b'\n')
        k = next((i for i in range(1, n+1) if p[i] != target[i]), None)
        if k is None:
            state['matching_graphs'] += 1
            if lengths and connected(adj):
                found = {'lengths': lengths, 'disjointness': adj, 'polynomial': p}
                return True
            return False
        if p[k] < target[k] or lengths.count(k) >= caps[k-1]:
            return False
        if lengths and k < lengths[-1]:
            raise ArithmeticError('cycle coefficient order invariant')
        allowed = [i for i, d in enumerate(lengths) if d+k <= n]
        newbit = 1 << len(lengths)
        for mask in range(1 << len(allowed)):
            neighbors = sum(1 << j for i, j in enumerate(allowed) if mask >> i & 1)
            added, q = [], p[:]
            for bits, weight, sign in cliques:
                if bits & ~neighbors:
                    continue
                if weight+k > n:
                    state['capacity_prunes'] += 1
                    break
                added.append((bits | newbit, weight+k, -sign))
                q[weight+k] -= sign
            else:
                newadj = [a | (newbit if neighbors >> i & 1 else 0)
                          for i, a in enumerate(adj)]+[neighbors]
                if visit(lengths+[k], newadj, cliques+added, q):
                    return True
        return False
    try:
        visit([], [], [(0, 0, 1)], [1]+[0]*n)
    except ModelBudgetExceeded:
        return {**state, 'completed': False, 'excluded': False, 'reason': 'cycle_budget',
                'target': target, 'transcript_sha256': transcript.hexdigest()}
    return {**state, 'completed': found is None, 'excluded': found is None,
            'reason': 'no_connected_cycle_graph' if found is None else 'connected_abstract_survivor',
            'target': target, 'survivor': found, 'transcript_sha256': transcript.hexdigest()}
