# SPDX-License-Identifier: GPL-3.0-or-later
"""Signed maps, terminal symmetries and the tangent-orientation double cover."""
from collections import Counter, deque
from ttauto_cert.exact_math import charpoly, identity, matrix_multiply
from ttauto_cert.ribbon_verify import automorphisms, cyclic_word, substitute
from .combinatorics import boundary_cycles, require
from .ribbon import SurfaceState, real_matrix


def ribbon_graph(graph):
    m = graph["main_edges"]
    stars = graph["cyclic_directions"]
    origins = {d: v for v, star in enumerate(stars) for d in star}
    polygons, at = [], 0
    for p in graph["inner_prongs"]:
        polygons.append(list(range(m+at+1, m+at+p+1)))
        at += p
    return {"main_edges": m, "cyclic_directions": stars,
            "endpoints": [[origins[e], origins[-e]] for e in range(1, m+len(stars)+1)],
            "polygons": polygons, "punctured": [int(c == "puncture") for c in graph["inner_colors"]],
            "labels": [0]*len(polygons)}


def verify_map(source, target, vertices, words):
    require(len(words) == source.e == target.e, "map edge count")
    require(sorted(vertices) == list(range(target.k)), "map vertex permutation")
    for e, word in enumerate(words, 1):
        require(word and all(type(d) is int and d in target.origins for d in word), "invalid signed word")
        require(target.origins[word[0]] == vertices[source.origins[e]] and
                target.origins[-word[-1]] == vertices[source.origins[-e]], "map endpoints")
        require(all(target.origins[-a] == target.origins[b] for a, b in zip(word, word[1:])), "discontinuous map")
    require(sorted(cyclic_word(substitute(words, w)) for w in source.faces) ==
            sorted(cyclic_word(w) for w in target.faces), "oriented peripheral transport")
    return True


def terminals(state):
    result = []
    for naming in automorphisms(ribbon_graph(state.graph)):
        words = [[e] for e in naming]
        vertices = [state.origins[naming[star[0]-1]] if star[0] > 0
                    else state.origins[-naming[-star[0]-1]]
                    for star in state.graph["cyclic_directions"]]
        verify_map(state, state, vertices, words)
        result.append({"vertex_images": vertices, "edge_images": words,
                       "matrix": real_matrix(words, state.m)})
    require(bool(result), "identity terminal missing")
    return result


def track_orientation(state):
    colors, todo = {0: 0}, deque([0])
    while todo:
        v = todo.popleft()
        for d in state.graph["cyclic_directions"][v]:
            w = state.origins[-d]
            if w in colors and colors[w] != 1-colors[v]:
                return None
            if w not in colors:
                colors[w] = 1-colors[v]
                todo.append(w)
    return [colors[v] for v in range(state.k)]


def orientation_homology(state, words, vertices):
    """Fill all lifted faces; positive lift preserves the tangent orientation.

    At every split switch all real directions lie on one side and all
    infinitesimal directions on the other. Crossing any band flips the local
    orientation label. Thus the double cover is the bipartite double cover of
    the full ribbon graph, not of its unsigned real transition matrix.
    """
    colors = track_orientation(state)
    if colors is not None:
        matrix, _ = state.homology(words)
        sign = 1 if colors[vertices[0]] == colors[0] else -1
        require(all((colors[vertices[v]] == colors[v]) == (sign == 1) for v in range(state.k)),
                "tangent orientation transport is inconsistent")
        positive = [[sign*x for x in row] for row in matrix]
        return {"genus": state.graph["genus"], "connected_double_cover": False,
                "base_track_orientable": True, "positive_lift_matrix": positive,
                "actual_base_matrix": matrix,
                "polynomial": charpoly(positive) if positive else [1],
                "base_action_sign": sign, "lifted_prongs": state.graph["inner_prongs"]+state.graph["outer_prongs"]}
    m, k, e = state.m, state.k, state.e
    def edge_label(edge, sheet):
        return 2*(edge-1)+sheet+1
    def dart(d, sheet):
        # Positive lifted edge starts on `sheet`; reversing it starts on the other.
        return edge_label(d, sheet) if d > 0 else -edge_label(-d, 1-sheet)
    stars = [[dart(d, sheet) for d in star]
             for star in state.graph["cyclic_directions"] for sheet in (0, 1)]
    faces = boundary_cycles(stars)
    genus_numerator = 2-2*k+2*e-len(faces)
    require(genus_numerator % 2 == 0, "cover genus")
    genus = genus_numerator//2
    inner = [w for w in faces if all(abs(d) > 2*m for d in w)]
    outer = [w for w in faces if w not in inner]
    prongs = [len(w) for w in inner]+[
        sum(abs(w[i]) <= 2*m and abs(w[(i+1) % len(w)]) <= 2*m for i in range(len(w)))
        for w in outer]
    graph = {"main_edges": 2*m, "genus": genus, "cyclic_directions": stars,
             "real_fans": [[d for d in star if abs(d) <= 2*m] for star in stars],
             "boundary_cycles": faces, "inner_prongs": [len(w) for w in inner],
             "outer_prongs": prongs[len(inner):], "inner_colors": ["filled"]*len(inner)}
    cover = SurfaceState(graph)
    lifted_words = []
    for word in words:
        for sheet in (0, 1):
            lifted_words.append([dart(d, (sheet+i) % 2) for i, d in enumerate(word)])
    lifted_vertices = [2*v+sheet for v in vertices for sheet in (0, 1)]
    verify_map(cover, cover, lifted_vertices, lifted_words)
    matrix, _ = cover.homology(lifted_words)
    deck_words = [[j+1 if j % 2 else j-1] for j in range(1, 2*e+1)]
    deck_vertices = [v ^ 1 for v in range(2*k)]
    verify_map(cover, cover, deck_vertices, deck_words)
    deck, _ = cover.homology(deck_words)
    require(matrix_multiply(deck, deck) == identity(len(deck)), "deck involution on H1")
    other = matrix_multiply(deck, matrix)
    require(other == matrix_multiply(matrix, deck), "lift must commute with deck action")
    require(len(deck)+sum(deck[i][i] for i in range(len(deck))) == 4*state.graph['genus'],
            "deck-invariant homology rank must equal base H1 rank")
    require(Counter(prongs) == Counter(q for p in state.graph["inner_prongs"]+state.graph["outer_prongs"]
                                      for q in ([2*p] if p % 2 else [p, p])), "cover prongs")
    return {"genus": genus, "connected_double_cover": True, "base_track_orientable": False,
            "positive_lift_matrix": matrix, "polynomial": charpoly(matrix),
            "deck_matrix": deck, "other_lift_matrix": other,
            "other_lift_polynomial": charpoly(other),
            "lifted_prongs": prongs, "base_action_sign": 1}
