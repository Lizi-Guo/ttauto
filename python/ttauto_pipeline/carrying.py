# SPDX-License-Identifier: GPL-3.0-or-later
"""Ordered ribbon-band carrying certificates allowing repeated route vertices.

Each target band has ordered parallel lanes. At a target switch, transit
strands must admit successive adjacent-pair removal, leaving exactly one
source switch star in its original cyclic order. This is a finite planar
disk-embedding certificate, including loops and return-to-source slides.
"""
from itertools import permutations, product
from .combinatorics import require as need

def cyclic(w):
    w = tuple(w)
    return min((w[i:] + w[:i] for i in range(len(w)))) if w else ()

def check_orders(source, target, arc, orders):
    words = arc['edge_images']
    vertex_map = arc['vertex_images']
    origins = {d: v for v, s in enumerate(target['cyclic_directions']) for d in s}
    records = []
    for v, star in enumerate(target['cyclic_directions']):
        ports = []
        for d in star:
            lanes = orders[abs(d)] if d > 0 else orders[abs(d)][::-1]
            for edge, j in lanes:
                word = words[edge - 1]
                start = word[j] == d
                need(start or word[j] == -d, 'lane does not traverse this band')
                if start:
                    label = ('H', edge) if j == 0 else ('T', edge, j - 1)
                else:
                    label = ('H', -edge) if j == len(word) - 1 else ('T', edge, j)
                ports.append(label)
        remaining = ports[:]
        removed = []
        while True:
            match = next((i for i, p in enumerate(remaining) if p[0] == 'T' and p == remaining[(i + 1) % len(remaining)]), None)
            if match is None:
                break
            removed.append(remaining[match])
            if match == len(remaining) - 1:
                remaining = remaining[1:-1]
            else:
                remaining = remaining[:match] + remaining[match + 2:]
        if any((p[0] != 'H' for p in remaining)):
            return None
        source_vertex = vertex_map.index(v)
        if cyclic((p[1] for p in remaining)) != cyclic(source['cyclic_directions'][source_vertex]):
            return None
        records.append(dict(target_vertex=v, source_vertex=source_vertex, ports=ports, adjacent_transit_pair_removals=removed, remaining_hub_ports=remaining))
    for edge, word in enumerate(words, 1):
        src = next((v for v, s in enumerate(source['cyclic_directions']) if edge in s))
        dst = next((v for v, s in enumerate(source['cyclic_directions']) if -edge in s))
        need(origins[word[0]] == vertex_map[src] and origins[-word[-1]] == vertex_map[dst], 'edge image has the wrong endpoints')
        need(all((origins[-a] == origins[b] for a, b in zip(word, word[1:]))), 'edge image is not a continuous path')
    return records

def certify(source, target, arc):
    m = source['main_edges']
    words = arc['edge_images']
    e = len(words)
    need(target['main_edges'] == m, 'source and target have different real dimensions')
    need(sorted(arc['vertex_images']) == list(range(len(source['cyclic_directions']))), 'vertex images are not a permutation')
    need(all(((abs(word[0]) <= m) == (i < m) and (abs(word[-1]) <= m) == (i < m) for i, word in enumerate(words))), 'real edges must begin and end on real bands, infinitesimal edges on infinitesimal bands')
    need(all(((abs(a) <= m) != (abs(b) <= m) for w in words for a, b in zip(w, w[1:]))), 'edge image does not alternate real and infinitesimal bands')
    uses = {d: [] for d in range(1, e + 1)}
    for edge, word in enumerate(words, 1):
        for j, d in enumerate(word):
            uses[abs(d)].append((edge, j))
    need(sorted(map(len, uses.values())) == [1] * (e - 2) + [2, 2], 'single paired slide only')
    choices = [list(permutations(uses[d])) for d in range(1, e + 1)]
    for assignment in product(*choices):
        orders = {d: list(row) for d, row in enumerate(assignment, 1)}
        record = check_orders(source, target, arc, orders)
        if record is not None:
            return dict(positive_band_lane_orders=orders, switch_disks=record, repeated_vertices_allowed=True, all_transit_chords_removed=True, all_source_cyclic_stars_preserved=True, smoothing_types_preserved=True)
    raise ValueError('no noncrossing ordered-band certificate')
