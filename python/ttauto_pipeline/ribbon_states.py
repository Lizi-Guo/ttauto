# SPDX-License-Identifier: GPL-3.0-or-later
"""General standard ribbon states with several inner infinitesimal polygons.

Finite exhaustive labelled pairing enumeration, modulo polygon rotations,
permutations of equal-type polygons, and signed real-edge renaming.
Allows loops and parallel edges; disconnected pairings are discarded.
"""
from itertools import permutations, product
from pathlib import Path
import argparse
import json
import time

from .combinatorics import boundary_cycles, compositions, matchings, require


def polygon_vertices(prongs):
    blocks, at = [], 0
    for p in prongs:
        blocks.append(list(range(at, at+p)))
        at += p
    return blocks


def graph_from_fans(prongs, fans):
    m, k = sum(map(len, fans))//2, len(fans)
    blocks = polygon_vertices(prongs)
    stars = [None]*k
    inners = []
    for block in blocks:
        inners.append([m+v+1 for v in block])
        for i, v in enumerate(block):
            stars[v] = [m+v+1, -(m+block[i-1]+1), *fans[v]]
    origins = {d: v for v, star in enumerate(stars) for d in star}
    seen, todo = {0}, [0]
    while todo:
        v = todo.pop()
        for d in stars[v]:
            w = origins[-d]
            if w not in seen:
                seen.add(w)
                todo.append(w)
    if len(seen) != k:
        return None
    faces = boundary_cycles(stars)
    require(all(inner in faces for inner in inners), 'inner face cycles')
    outer = [w for w in faces if w not in inners]
    require(all(not (abs(w[i]) > m and abs(w[(i+1) % len(w)]) > m)
                for w in outer for i in range(len(w))), 'outer II corner')
    out_prongs = sorted(sum(abs(w[i]) <= m and abs(w[(i+1) % len(w)]) <= m
                            for i in range(len(w))) for w in outer)
    chi = -m+len(faces)
    require((2-chi) % 2 == 0, 'integer genus')
    return {'inner_prongs': list(prongs), 'outer_prongs': out_prongs,
            'main_edges': m, 'genus': (2-chi)//2, 'real_fans': [list(f) for f in fans],
            'cyclic_directions': stars, 'boundary_cycles': faces, 'connected': True}


def symmetry_orders(prongs, colors=None):
    colors = tuple(colors or ["filled"]*len(prongs))
    blocks = polygon_vertices(prongs)
    for perm in permutations(range(len(prongs))):
        if any((prongs[i], colors[i]) != (prongs[j], colors[j]) for i, j in enumerate(perm)):
            continue
        for rotations in product(*(range(p) for p in prongs)):
            order = []
            for i, shift in zip(perm, rotations):
                block = blocks[i]
                order.extend(block[shift:]+block[:shift])
            yield order


def canonical(fans, orders):
    best = None
    for order in orders:
        names, rows = {}, []
        for v in order:
            row = []
            for d in fans[v]:
                if abs(d) not in names:
                    names[abs(d)] = (len(names)+1)*(1 if d > 0 else -1)
                row.append(names[abs(d)]*(1 if d > 0 else -1))
            rows.append(tuple(row))
        code = tuple(rows)
        if best is None or code < best:
            best = code
    return best


def rooted_transports(prongs, colors, fans):
    """Canonicalize a connected ribbon graph with O(vertices) rooted walks.

    A root fixes the first polygon rotation. At first crossing to a new
    polygon, its entry prong fixes its rotation. This avoids materializing
    factorially many permutations of equal polygons (notably punctures).
    The color and prong size remain in the code and polygon blocks are then
    stably grouped by type, so only allowed isomorphisms are identified.
    """
    blocks = polygon_vertices(prongs)
    owner = {v: p for p, block in enumerate(blocks) for v in block}
    origin = {d: v for v, fan in enumerate(fans) for d in fan}
    for root in range(len(fans)):
        seen, ordered_blocks, pending = set(), [], [root]
        for first in pending:
            p = owner[first]
            if p in seen:
                continue
            seen.add(p)
            block = blocks[p]
            at = block.index(first)
            rotated = block[at:]+block[:at]
            ordered_blocks.append((prongs[p], colors[p], rotated))
            pending.extend(origin[-d] for v in rotated for d in fans[v] if owner[origin[-d]] not in seen)
        if len(seen) != len(blocks):
            continue
        ordered_blocks.sort(key=lambda row: row[:2])
        order = [v for _, _, block in ordered_blocks for v in block]
        names, code = {}, []
        for v in order:
            row = []
            for d in fans[v]:
                if abs(d) not in names:
                    names[abs(d)] = (len(names)+1)*(1 if d > 0 else -1)
                row.append(names[abs(d)]*(1 if d > 0 else -1))
            code.append(tuple(row))
        vertices = [order.index(v) for v in range(len(fans))]
        m = len(names)
        words = [[names[e]] for e in range(1, m+1)]+[[m+v+1] for v in vertices]
        yield tuple(code), vertices, words


class StateBudgetExceeded(RuntimeError):
    """State construction stopped; no completeness claim is possible."""


def enumerate_states(genus, inner, outer, colors=None, max_decorations=None):
    colors = tuple(colors or ["filled"]*len(inner))
    require(list(zip(inner, colors)) == sorted(zip(inner, colors)), 'inner data must be sorted')
    inner, outer = tuple(inner), tuple(sorted(outer))
    require(inner and outer and all(p > 0 for p in inner+outer), 'nonempty positive puncture partitions')
    require(sum(p-2 for p in inner+outer) == 4*genus-4, 'Euler prong relation')
    m, k = 2*genus-2+len(inner)+len(outer), sum(inner)
    require(k <= 2*m, 'positive fans impossible')
    accepted = disconnected = examined = 0
    found = {}
    for degrees in compositions(2*m, k):
        for matching in matchings(tuple(range(2*m))):
            if max_decorations is not None and examined >= max_decorations:
                raise StateBudgetExceeded(f"state decoration budget exhausted at {examined}")
            examined += 1
            slots = [None]*(2*m)
            for e, (a, b) in enumerate(matching, 1):
                slots[a], slots[b] = e, -e
            fans, at = [], 0
            for d in degrees:
                fans.append(slots[at:at+d])
                at += d
            graph = graph_from_fans(inner, fans)
            if graph is None:
                disconnected += 1
                continue
            if graph['genus'] != genus or graph['outer_prongs'] != list(outer):
                continue
            accepted += 1
            code = min(rooted_transports(inner, colors, fans))[0]
            if code not in found:
                found[code] = graph_from_fans(inner, code)
    graphs = [found[c] for c in sorted(found)]
    for graph in graphs:
        graph['inner_colors'] = list(colors)
    return graphs, {'labelled_pairing_decorations': examined,
            'disconnected_decorations': disconnected, 'accepted_labelled_decorations': accepted,
            'oriented_isomorphism_classes': len(found), 'canonicalization': 'rooted_oriented_polygon_walks'}

