# SPDX-License-Identifier: GPL-3.0-or-later
"""Finite standard ribbon states and exact proper-prefix cone experiments.

Enumerates every state with ONE infinitesimal polygon, prescribed genus and
outer prongs, modulo oriented polygon rotations and real-edge renaming.
This is a topological state enumeration, not pA-map enumeration.
"""
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path


def require(value, message):
    if not value:
        raise ValueError(message)


def compositions(total, size):
    if size == 1:
        if total >= 1:
            yield (total,)
        return
    for first in range(1, total-size+2):
        for tail in compositions(total-first, size-1):
            yield (first,)+tail


def matchings(items):
    if not items:
        yield ()
        return
    a = items[0]
    for i in range(1, len(items)):
        for rest in matchings(items[1:i]+items[i+1:]):
            yield ((a, items[i]),)+rest


def boundary_cycles(stars):
    predecessor = {d: star[i-1] for star in stars for i, d in enumerate(star)}
    require(len(predecessor) == sum(map(len, stars)), 'duplicate dart')
    require(set(predecessor) == {-d for d in predecessor}, 'unpaired dart')
    todo, result = set(predecessor), []
    while todo:
        start = min(todo)
        word, d = [], start
        while d in todo:
            todo.remove(d)
            word.append(d)
            d = predecessor[-d]
        require(d == start, 'noncyclic boundary permutation')
        result.append(word)
    return result


def canonical_fans(fans):
    codes = []
    for shift in range(len(fans)):
        naming, code = {}, []
        for fan in fans[shift:]+fans[:shift]:
            row = []
            for d in fan:
                if abs(d) not in naming:
                    naming[abs(d)] = (len(naming)+1)*(1 if d > 0 else -1)
                row.append(naming[abs(d)]*(1 if d > 0 else -1))
            code.append(tuple(row))
        codes.append(tuple(code))
    return min(codes)


def graph_from_fans(fans):
    k = len(fans)
    m = sum(map(len, fans))//2
    stars = [[m+v+1, -(m+(v-1) % k+1), *fan] for v, fan in enumerate(fans)]
    cycles = boundary_cycles(stars)
    inner = list(range(m+1, m+k+1))
    require(inner in cycles, 'positive inner polygon missing')
    outer = [word for word in cycles if word != inner]
    require(all(not (abs(w[i]) > m and abs(w[(i+1) % len(w)]) > m)
                for w in outer for i in range(len(w))), 'outer I-I turn')
    prongs = sorted(sum(abs(w[i]) <= m and abs(w[(i+1) % len(w)]) <= m
                        for i in range(len(w))) for w in outer)
    euler = k-(m+k)+len(cycles)
    require((2-euler) % 2 == 0, 'nonintegral genus')
    return {'real_fans': [list(f) for f in fans], 'cyclic_directions': stars,
            'main_edges': m, 'inner_prongs': [k], 'outer_prongs': prongs,
            'genus': (2-euler)//2, 'boundary_cycles': cycles}


def enumerate_states(genus, inner_prongs, outer_prongs):
    # Single inner polygon makes connectedness automatic; real slots are
    # distinguishable until the explicit oriented isomorphism quotient.
    k = inner_prongs
    m = 2*genus-2+1+len(outer_prongs)
    require(k >= 1 and m >= 1, 'invalid dimensions')
    require(k-2+sum(p-2 for p in outer_prongs) == 4*genus-4, 'prong Euler relation')
    found, raw, matching_count = {}, 0, 0
    all_matchings = list(matchings(tuple(range(2*m))))
    for degrees in compositions(2*m, k):
        for pairing in all_matchings:
            matching_count += 1
            slots = [None]*(2*m)
            for e, (a, b) in enumerate(pairing, 1):
                slots[a], slots[b] = e, -e
            fans, at = [], 0
            for degree in degrees:
                fans.append(tuple(slots[at:at+degree]))
                at += degree
            graph = graph_from_fans(fans)
            if graph['genus'] != genus or graph['outer_prongs'] != sorted(outer_prongs):
                continue
            raw += 1
            code = canonical_fans(fans)
            if code not in found:
                found[code] = graph_from_fans(list(code))
    return [found[c] for c in sorted(found)], {
        'slot_pairings_per_composition': len(all_matchings),
        'labelled_pairing_decorations': matching_count,
        'accepted_labelled_decorations': raw,
        'oriented_isomorphism_classes': len(found)}


def rref(rows, width):
    a = [[Fraction(x) for x in row] for row in rows]
    pivot_cols, at = [], 0
    for col in range(width):
        pivot = next((r for r in range(at, len(a)) if a[r][col]), None)
        if pivot is None:
            continue
        a[at], a[pivot] = a[pivot], a[at]
        scalar = a[at][col]
        a[at] = [v/scalar for v in a[at]]
        for r in range(len(a)):
            if r != at and a[r][col]:
                scalar = a[r][col]
                a[r] = [x-scalar*y for x, y in zip(a[r], a[at])]
        pivot_cols.append(col)
        at += 1
        if at == len(a):
            break
    return a, pivot_cols

