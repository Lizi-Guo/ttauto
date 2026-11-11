# SPDX-License-Identifier: GPL-3.0-or-later
"""Geometric validation is mandatory for attainment, independent of pruning."""
from collections import deque
from ttauto_cert.exact_math import matrix_multiply, compare_largest_real_roots
from ttauto_cert.ribbon_verify import cyclic_word, substitute
from .combinatorics import boundary_cycles, require
from .gates import audit as gate_audit


def canonical_carrier(context):
    """Reconstruct the canonical BH infinitesimal graph, not saved face labels.

    Supported local patterns: a cycle or a cycle with one edge omitted,
    including the two-gate single-edge junction. Other patterns remain
    unresolved; they are never discarded as a claimed non-pA map.
    """
    state = context.state
    words, vertices = context.mapping
    try:
        audit = gate_audit(state.graph, words, vertices)
    except ValueError as error:
        return {"status": "unresolved", "reason": "filled_spine_not_certified", "detail": str(error)}
    if not audit["all_canonical_gate_graphs_connected"]:
        return {"status": "unresolved", "reason": "canonical_gate_graph_disconnected", "audit": audit}
    m = state.m
    stars, dart_gate, side_lookup, inf_edges = [], {}, {}, []
    for record in audit["vertices"]:
        gates, order = record["gates"], record["cyclic_gate_order"]
        n, offset = len(gates), len(stars)
        if n < 2:
            return {"status": "unresolved", "reason": "fewer_than_two_gates"}
        expected = {tuple(sorted((order[i], order[(i+1) % n]))) for i in range(n)}
        actual = set(map(tuple, record["canonical_taken_turn_edges"]))
        if not actual <= expected or len(actual) not in ({1} if n == 2 else {n-1, n}):
            return {"status": "unresolved", "reason": "unsupported_canonical_gate_pattern", "audit": audit}
        for gi in order:
            members = set(gates[gi])
            star = record["cyclic_directions"]
            start = next(i for i, d in enumerate(star) if d in members and star[i-1] not in members)
            fan = [d for d in star[start:]+star[:start] if d in members]
            vertex = len(stars)
            stars.append(fan)
            dart_gate.update({d: vertex for d in fan})
        added = set()
        for i, gi in enumerate(order):
            gj = order[(i+1) % n]
            pair = tuple(sorted((gi, gj)))
            if pair not in actual or pair in added:
                continue
            added.add(pair)
            a, b = offset+i, offset+(i+1) % n
            label = m+len(inf_edges)+1
            inf_edges.append((a, b))
            # Outgoing side, reversed incoming side, real fan: the original
            # cyclic smoothing convention, with absent sides simply removed.
            stars[a].insert(0, label)
            insertion = 1 if stars[b] and stars[b][0] > m else 0
            stars[b].insert(insertion, -label)
            side_lookup[a, b], side_lookup[b, a] = label, -label
    real_words = audit["main_edge_images"]
    derivative = {d: (real_words[d][0] if d > 0 else -real_words[-d][-1]) for d in dart_gate}
    vertex_images = []
    for star in stars:
        images = {dart_gate[derivative[d]] for d in star if abs(d) <= m}
        require(len(images) == 1, "canonical gate image")
        vertex_images.append(images.pop())
    require(sorted(vertex_images) == list(range(len(stars))), "canonical gate permutation")
    canonical_words = []
    for e in range(1, m+1):
        old = real_words[e]
        word = [old[0]]
        for a, b in zip(old, old[1:]):
            word += [side_lookup[dart_gate[-a], dart_gate[b]], b]
        canonical_words.append(word)
    for a, b in inf_edges:
        canonical_words.append([side_lookup[vertex_images[a], vertex_images[b]]])
    origins = {d: v for v, star in enumerate(stars) for d in star}
    for e, word in enumerate(canonical_words, 1):
        require(origins[word[0]] == vertex_images[origins[e]] and
                origins[-word[-1]] == vertex_images[origins[-e]], "canonical endpoints")
        require(all(origins[-a] == origins[b] for a, b in zip(word, word[1:])), "canonical continuity")
    faces = boundary_cycles(stars)
    require(len(stars)-len(canonical_words)+len(faces) == 2-2*state.graph["genus"], "canonical genus")
    keys = [cyclic_word(w) for w in faces]
    require(sorted(cyclic_word(substitute(canonical_words, w)) for w in faces) == sorted(keys), "canonical peripherals")
    full = [[sum(abs(d) == i+1 for d in word) for word in canonical_words] for i in range(len(canonical_words))]
    switch = [[sum((1 if abs(d) <= m else -1) for d in star if abs(d) == j+1)
               for j in range(len(full))] for star in stars]
    vm = [[int(i == vertex_images[j]) for j in range(len(stars))] for i in range(len(stars))]
    require(matrix_multiply(switch, full) == matrix_multiply(vm, switch), "canonical switch intertwining")
    reached = set(range(m))
    while True:
        following = reached | {i for i in range(len(full)) for j in reached if full[i][j]}
        if following == reached:
            break
        reached = following
    require(len(reached) == len(full), "zero-width infinitesimal band")
    prongs = [sum((abs(w[i]) <= m) == (abs(w[(i+1) % len(w)]) <= m) for i in range(len(w))) for w in faces]
    if any(p < 2 for p in prongs):
        return {"status": "unresolved", "reason": "canonical_face_cannot_be_filled", "actual_prongs": prongs}
    # Tangent orientability is computed on the reconstructed carrier.
    colors, pending, orientable = {0: 0}, deque([0]), True
    while pending:
        v = pending.popleft()
        for d in stars[v]:
            w = origins[-d]
            if w in colors:
                orientable &= colors[w] != colors[v]
            else:
                colors[w] = 1-colors[v]
                pending.append(w)
    target = [p for p in state.graph["inner_prongs"]+state.graph["outer_prongs"] if p > 2]
    actual = sorted(p for p in prongs if p > 2)
    desired = context.config["surface"]["foliation"]
    matched = actual == sorted(target) and (desired == "any" or orientable == (desired == "orientable"))
    return {"status": "attained" if matched else "rejected", "reason": "canonical_BH_carrier" if matched else "actual_stratum_mismatch",
            "actual_singular_prongs": actual, "actual_foliation_orientable": orientable,
            "canonical_stars": stars, "canonical_edge_images": canonical_words,
            "canonical_vertex_images": vertex_images, "canonical_faces": faces,
            "canonical_face_prongs": prongs, "all_iterated_turns_legal": True,
            "positive_full_PF_measure": True,
            "theorem_dependency": "Bestvina--Handel (1995), sections 3.3--3.4: efficient surface carrier and canonical filling train track",
            "scope": "this signed map only; no normal-form coverage or global minimum follows"}


def pillowcase(context):
    """Certify only the sphere with exactly four marked one-prong points.

    Its branched double cover is a torus. Hyperbolicity of the induced integral
    H1 action classifies this mapping class; this is not a criterion for higher
    genus or additional punctures. The cover is built from the actual signed map.
    """
    surface = context.config["surface"]
    if not (surface["kind"] == "disc" and surface["punctures"] == 3 and
            surface["inner_prongs"] == [1, 1, 1] and surface["outer_prongs"] == [1]
            and surface["inner_colors"] == ["puncture"]*3):
        return {"status": "unresolved", "reason": "unsupported_disc_geometry",
                "signed_surface_map_verified": True,
                "scope": "no general punctured-disc classification is implemented"}
    cover = context.cover
    require(cover["connected_double_cover"] and cover["genus"] == 1 and
            sorted(cover["lifted_prongs"]) == [2]*4, "pillowcase torus cover")
    matrix = cover["positive_lift_matrix"]
    require(len(matrix) == 2 and all(len(row) == 2 for row in matrix) and
            all(type(x) is int for row in matrix for x in row), "integral rank-two cover action")
    a, b = matrix[0]
    c, d = matrix[1]
    require(a*d-b*c == 1, "torus action determinant")
    trace = a+d
    evidence = {"cover": cover, "trace": trace, "determinant": 1,
                "theorem_dependency": "Farb--Margalit, A Primer on Mapping Class Groups, sections 2.2.5 and 13.1: four-punctured-sphere branched cover and torus classification",
                "scope": "this signed D3 map, with the boundary treated as the fixed fourth mark; no global minimum follows"}
    if abs(trace) <= 2:
        return {**evidence, "status": "rejected", "reason": "nonhyperbolic_torus_action"}
    polynomial = [1, -abs(trace), 1]
    if compare_largest_real_roots(polynomial, context.candidate["polynomial"]) != 0:
        return {**evidence, "status": "rejected", "reason": "actual_dilatation_mismatch",
                "actual_dilatation_polynomial": polynomial}
    if surface["foliation"] == "orientable":
        return {**evidence, "status": "rejected", "reason": "actual_stratum_mismatch"}
    return {**evidence, "status": "attained", "reason": "pillowcase_hyperbolic_torus_action",
            "actual_puncture_prongs": [1]*3, "actual_outer_prongs": 1,
            "actual_interior_prongs": [], "actual_foliation_orientable": False,
            "actual_dilatation_polynomial": polynomial}


def classify(context):
    if context.config["geometry"]["kind"] == "audit":
        return {"status": "unresolved", "reason": "audit_only",
                "signed_surface_map_verified": True}
    if context.config["surface"]["kind"] == "closed":
        return canonical_carrier(context)
    return pillowcase(context)
