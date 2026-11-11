# SPDX-License-Identifier: GPL-3.0-or-later
"""Experimental closed-surface automata with multiple inner/outer polygons.

General connected ribbon-state backend, all face transports, and full
homology quotient after filling every puncture. The general macro coverage
proof is in closed_surface_general_prefix_proof_20260925_zh.md; this program
builds one prescribed topological type, not every type on a closed surface.
"""
from collections import deque
from fractions import Fraction
from itertools import product
from pathlib import Path
import argparse
import json
import random
import time

from ttauto_cert.exact_math import charpoly, identity, matrix_multiply
from ttauto_cert.standard_maps import primitive
from ttauto_cert.ribbon_verify import cyclic_word, inverse, reduce_word, substitute
from .ribbon_states import enumerate_states, graph_from_fans, symmetry_orders, rooted_transports
from .combinatorics import require, rref


def real_matrix(words, m):
    return [[sum(abs(d) == i+1 for d in word) for word in words[:m]] for i in range(m)]



class SurfaceState:
    def __init__(self, graph):
        self.graph = graph
        self.m = graph['main_edges']
        self.k = len(graph['real_fans'])
        self.e = self.m+self.k
        self.origins = {d: v for v, star in enumerate(graph['cyclic_directions']) for d in star}
        self.faces = graph['boundary_cycles']
        self.outer = [w for w in self.faces if any(abs(d) <= self.m for d in w)]
        paths, tree, queue = {0: []}, set(), deque([0])
        while queue:
            v = queue.popleft()
            for d in graph['cyclic_directions'][v]:
                other = self.origins[-d]
                if other not in paths:
                    paths[other] = paths[v]+[d]
                    tree.add(abs(d))
                    queue.append(other)
        require(len(paths) == self.k, 'connected source graph')
        kept = [e for e in range(1, self.e+1) if e not in tree]
        self.generators = {e: i+1 for i, e in enumerate(kept)}
        self.basis = [paths[self.origins[e]]+[e]+inverse(paths[self.origins[-e]]) for e in kept]
        self.rank = len(kept)
        require(self.rank == self.m+1, 'graph rank')
        relations = [self.abelian(self.project(w)) for w in self.faces]
        reduced, self.pivots = rref(relations, self.rank)
        self.relations = reduced[:len(self.pivots)]
        self.free = [i for i in range(self.rank) if i not in self.pivots]
        require(len(self.free) == 2*graph['genus'], 'closed homology rank')
        require(all(x.denominator == 1 for row in self.relations for x in row), 'integral peripheral quotient')
        self.relations = [[int(x) for x in row] for row in self.relations]

    def project(self, word):
        return reduce_word([self.generators[abs(d)]*(1 if d > 0 else -1)
                            for d in word if abs(d) in self.generators])

    def abelian(self, word):
        vector = [0]*self.rank
        for d in word:
            vector[abs(d)-1] += 1 if d > 0 else -1
        return vector

    def quotient(self, vector):
        return [vector[j]-sum(vector[p]*row[j] for p, row in zip(self.pivots, self.relations))
                for j in self.free]

    def homology(self, words):
        columns = [self.abelian(self.project(substitute(words, w))) for w in self.basis]
        for relation in self.relations:
            image = [sum(columns[j][i]*relation[j] for j in range(self.rank)) for i in range(self.rank)]
            require(not any(self.quotient(image)), 'filled relations not invariant')
        result = [[self.quotient(columns[j])[i] for j in self.free] for i in range(len(self.free))]
        return result, {'graph_basis': self.basis, 'free_graph_edges': list(self.generators),
                        'filled_relations_rref': self.relations, 'pivot_indices': self.pivots,
                        'kept_basis_indices': self.free}


def canonical_transport(fans, m, orders):
    records = []
    for order in orders:
        names, code = {}, []
        for v in order:
            row = []
            for d in fans[v]:
                if abs(d) not in names:
                    names[abs(d)] = (len(names)+1)*(1 if d > 0 else -1)
                row.append(names[abs(d)]*(1 if d > 0 else -1))
            code.append(tuple(row))
        vertices = [order.index(v) for v in range(len(fans))]
        words = [[names[e]] for e in range(1, m+1)]+[[m+v+1] for v in vertices]
        records.append((tuple(code), vertices, words))
    return min(records)


def slides(state, orders=None):
    for boundary in state.outer:
        for outer in (boundary, inverse(boundary)):
            for i, incoming in enumerate(outer):
                b, j = outer[(i+1) % len(outer)], outer[(i+2) % len(outer)]
                a = -incoming
                if not (abs(a) <= state.m and abs(b) <= state.m < abs(j) and abs(a) != abs(b)):
                    continue
                fans = [fan[:] for fan in state.graph['real_fans']]
                fans[state.origins[a]].remove(a)
                require(all(fans), 'empty real fan')
                destination = state.origins[-j]
                fans[destination].insert(len(fans[destination]) if j < 0 else 0, a)
                middle = graph_from_fans(state.graph['inner_prongs'], fans)
                require(middle and middle['genus'] == state.graph['genus'] and
                        middle['outer_prongs'] == state.graph['outer_prongs'], 'slide changes surface')
                q = [[e] for e in range(1, state.e+1)]
                q[abs(a)-1] = [b, j, a] if a > 0 else inverse([b, j, a])
                back = [[e] for e in range(1, state.e+1)]
                back[abs(a)-1] = [-j, -b, a] if a > 0 else inverse([-j, -b, a])
                require(all(reduce_word(substitute(q, back[e-1])) == [e] and reduce_word(substitute(back, q[e-1])) == [e]
                            for e in range(1, state.e+1)), 'two-sided slide inverse')
                require(sorted(cyclic_word(substitute(q, w)) for w in state.faces) ==
                        sorted(cyclic_word(w) for w in middle['boundary_cycles']), 'all oriented peripherals')
                code, vertices, naming = min(rooted_transports(state.graph['inner_prongs'], state.graph['inner_colors'], fans))
                words = [substitute(naming, w) for w in q]
                yield code, {'long_edge': a, 'short_edge': b, 'infinitesimal_edge': j,
                             'vertex_images': vertices, 'edge_images': words,
                             'matrix': real_matrix(words, state.m)}


def build(inner, outer, genus=2, colors=None, max_decorations=None):
    graphs, counts = enumerate_states(genus, inner, outer, colors, max_decorations)
    states = [SurfaceState(g) for g in graphs]
    indices = {tuple(map(tuple, s.graph['real_fans'])): i for i, s in enumerate(states)}
    arcs = []
    for source in states:
        row = []
        for code, arc in slides(source):
            require(code in indices, 'state closure')
            arc['target'] = indices[code]
            target = states[arc['target']]
            for e, word in enumerate(arc['edge_images'], 1):
                require(target.origins[word[0]] == arc['vertex_images'][source.origins[e]] and
                        target.origins[-word[-1]] == arc['vertex_images'][source.origins[-e]], 'map endpoints')
                require(all(target.origins[-x] == target.origins[y] for x, y in zip(word, word[1:])), 'path continuity')
            target_faces = [cyclic_word(w) for w in target.faces]
            arc['face_transport'] = [target_faces.index(cyclic_word(substitute(arc['edge_images'], w))) for w in source.faces]
            require(sorted(arc['face_transport']) == list(range(len(source.faces))), 'peripheral bijection')
            row.append(arc)
        arcs.append(row)
    return states, arcs, counts

