# SPDX-License-Identifier: GPL-3.0-or-later
"""Canonical gate audit adapted from the independently checked filled-spine audit."""
from collections import Counter
from ttauto_cert.exact_math import matrix_multiply as mm, charpoly as cp
from ttauto_cert.ribbon_verify import cyclic_word as conjugacy, substitute as apply
from .ribbon import SurfaceState
from .topology import verify_map as _verify

def topology(graph):
    return SurfaceState(graph)

def verify_map(source, target, vertices, words):
    return _verify(source, target, vertices, words)

def need(condition, message):
    if not condition:
        raise ValueError(message)

def components(n, edges):
    labels = list(range(n))
    for a, b in edges:
        x, y = (labels[a], labels[b])
        labels = [x if z == y else z for z in labels]
    return sorted([sorted((i for i, z in enumerate(labels) if z == x)) for x in set(labels)])

def audit(graph, words, vertices):
    m = graph['main_edges']
    inner = graph['inner_prongs']
    k = sum(inner)
    state = topology(graph)
    orig, faces = (state.origins, state.faces)
    verify_map(state, state, vertices, words)
    need(len(words) == m + k and len(vertices) == k, 'canonical gate invariant')
    need(sorted((conjugacy(apply(words, w)) for w in faces)) == sorted((conjugacy(w) for w in faces)), 'canonical gate invariant')
    need(all((len(w) == 1 and abs(w[0]) > m for w in words[m:])), 'canonical gate invariant')
    need(sorted((abs(w[0]) for w in words[m:])) == list(range(m + 1, m + k + 1)), 'canonical gate invariant')
    stars = []
    vertex_of = {}
    old_vertex = {}
    at = 0
    for v, p in enumerate(inner):
        star = []
        for i in range(at, at + p):
            old_vertex[i] = v
            star += graph['real_fans'][i]
        stars.append(star)
        at += p
        for d in star:
            need(d not in vertex_of, 'canonical gate invariant')
            vertex_of[d] = v
    need(Counter(vertex_of.keys()) == Counter({d: 1 for e in range(1, m + 1) for d in (e, -e)}), 'canonical gate invariant')
    next_dart = {d: star[(i - 1) % len(star)] for star in stars for i, d in enumerate(star)}
    unvisited = set(vertex_of)
    boundaries = []
    while unvisited:
        start = min(unvisited)
        d = start
        w = []
        while d in unvisited:
            unvisited.remove(d)
            w.append(d)
            d = next_dart[-d]
        need(d == start, 'canonical gate invariant')
        boundaries.append(w)
    need(len(stars) - m + len(boundaries) == 2 - 2 * graph['genus'], 'canonical gate invariant')
    need(sorted((conjugacy([d for d in w if abs(d) <= m]) for w in faces if any((abs(d) <= m for d in w)))) == sorted((conjugacy(w) for w in boundaries)), 'canonical gate invariant')
    collapsed_vm = []
    for v, p in enumerate(inner):
        targets = {old_vertex[vertices[i]] for i in range(sum(inner[:v]), sum(inner[:v]) + p)}
        need(len(targets) == 1, 'canonical gate invariant')
        collapsed_vm.append(targets.pop())
    need(sorted(collapsed_vm) == list(range(len(inner))), 'canonical gate invariant')
    images = {e: [d for d in words[e - 1] if abs(d) <= m] for e in range(1, m + 1)}
    need(all(images.values()), 'canonical gate invariant')
    for e, w in list(images.items()):
        images[-e] = [-d for d in w[::-1]]
    for d, w in images.items():
        need(vertex_of[w[0]] == collapsed_vm[vertex_of[d]] and vertex_of[-w[-1]] == collapsed_vm[vertex_of[-d]], 'canonical gate invariant')
        need(all((vertex_of[-a] == vertex_of[b] for a, b in zip(w, w[1:]))), 'canonical gate invariant')
        need(all((a != -b for a, b in zip(w, w[1:]))), 'canonical gate invariant')
    derivative = {d: w[0] for d, w in images.items()}
    turns = {tuple(sorted((-a, b))) for e in range(1, m + 1) for a, b in zip(images[e], images[e][1:])}
    queue = list(turns)
    while queue:
        a, b = queue.pop()
        need(a != b, 'canonical gate invariant')
        t = tuple(sorted((derivative[a], derivative[b])))
        need(t[0] != t[1], 'canonical gate invariant')
        if t not in turns:
            turns.add(t)
            queue.append(t)

    def equivalent(a, b):
        seen = set()
        while a != b and (a, b) not in seen:
            seen.add((a, b))
            a, b = (derivative[a], derivative[b])
        return a == b
    local = []
    for v, star in enumerate(stars):
        remaining = set(star)
        gates = []
        while remaining:
            a = min(remaining)
            g = sorted((b for b in remaining if equivalent(a, b)))
            need(a in g, 'canonical gate invariant')
            remaining -= set(g)
            gates.append(g)
        lookup = {d: i for i, g in enumerate(gates) for d in g}
        sequence = [lookup[d] for d in star]
        runs = [sequence[i] for i in range(len(sequence)) if sequence[i] != sequence[i - 1]]
        if not runs:
            runs = [sequence[0]]
        need(Counter(runs) == Counter(range(len(gates))), 'noncontiguous derivative gate')
        positions = {g: i for i, g in enumerate(runs)}
        edges = {tuple(sorted((lookup[a], lookup[b]))) for a, b in turns if vertex_of[a] == v}
        need(all((a != b for a, b in edges)), 'canonical gate invariant')
        for a, b in edges:
            x, y = sorted((positions[a], positions[b]))
            for c, d in edges:
                if len({a, b, c, d}) < 4:
                    continue
                need((x < positions[c] < y) == (x < positions[d] < y), 'crossing canonical infinitesimal chords')
        for a, b in ((a, b) for a in star for b in star):
            need(equivalent(a, b) == equivalent(derivative[a], derivative[b]), 'canonical gate invariant')
        local.append(dict(vertex=v, cyclic_directions=star, gates=gates, cyclic_gate_order=runs, canonical_taken_turn_edges=sorted(edges), canonical_gate_components=components(len(gates), edges), contiguous_gate_blocks=True, infinitesimal_chords_noncrossing=True))
    matrix = [[sum((abs(d) == i + 1 for d in images[j + 1])) for j in range(m)] for i in range(m)]
    power = [[int(i == j) for j in range(m)] for i in range(m)]
    positive = None
    for n in range(1, (m - 1) ** 2 + 2):
        power = mm(matrix, power)
        if all((x > 0 for row in power for x in row)):
            positive = n
            break
    need(positive is not None, 'canonical gate invariant')
    return dict(passed=True, genus=graph['genus'], filled_spine_vertices=len(inner), filled_spine_edges=m, filled_spine_faces=boundaries, main_edge_images={e: images[e] for e in range(1, m + 1)}, turn_orbit=sorted(turns), vertices=local, real_polynomial=list(cp(matrix)), positive_matrix_power=positive, all_edge_iterates_immersed=True, all_canonical_gate_graphs_connected=all((len(v['canonical_gate_components']) == 1 for v in local)), criterion_acceptance_pending=True, actual_stratum_certified=False, minimum_proved=False)
