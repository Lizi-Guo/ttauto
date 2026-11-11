"""Exact gate blow-up and complementary-region certificates.

The gate routines do NOT classify a mapping class as pseudo-Anosov. They reconstruct
the cusped ribbon graph determined by a graph map and checks a reported stratum
against that graph.  No Hall type, growth estimate, gate list, or prong count is
used to construct it.  The restricted local model is a cyclic interval gate
partition whose used-turn graph is the boundary of a polygon, possibly with
one edge missing. Unsupported local models raise UnsupportedTopology. A separate
routine certifies explicitly given hyperbolic torus/pillowcase affine maps.

The proof and the remaining realization/provenance hypotheses are documented in
review/topology_kernel_20260915_zh.md.  In particular, connected local Whitehead
graphs are NOT asserted to be a general pseudo-Anosov criterion.

SPDX-License-Identifier: GPL-3.0-or-later
"""
from itertools import combinations
from fractions import Fraction


class TopologyCertificateError(ValueError):
    """Malformed graph, inconsistent report, or failed necessary check."""


class UnsupportedTopology(TopologyCertificateError):
    """A local model outside the explicitly proved reconstruction class."""


def _require(condition, message):
    if not condition:
        raise TopologyCertificateError(message)


def _supported(condition, message):
    if not condition:
        raise UnsupportedTopology(message)


def _integer(value):
    return type(value) is int


def _inverse(word):
    return [-d for d in reversed(word)]


def _cyclic_key(word):
    """Oriented cyclic word; reversal is deliberately NOT identified."""
    word = tuple(word)
    return min(word[i:] + word[:i] for i in range(len(word))) if word else ()


def _cyclic_reduce(word):
    reduced = []
    for d in word:
        if reduced and reduced[-1] == -d:
            reduced.pop()
        else:
            reduced.append(d)
    left, right = 0, len(reduced)
    while right - left >= 2 and reduced[left] == -reduced[right - 1]:
        left += 1
        right -= 1
    return reduced[left:right]


def _boundary_cycles(rotations):
    successor = {}
    for cyclic in rotations.values():
        _require(bool(cyclic), 'empty rotation at a vertex')
        for i, d in enumerate(cyclic):
            _require(d not in successor, 'direction occurs in two rotations')
            successor[d] = cyclic[(i + 1) % len(cyclic)]
    _require(set(successor) == {-d for d in successor}, 'a direction is missing')
    unused = set(successor)
    boundaries = []
    while unused:
        first = min(unused)
        d = first
        boundary = []
        while d in unused:
            unused.remove(d)
            boundary.append(d)
            d = successor[-d]
        _require(d == first, 'boundary permutation failed to close')
        boundaries.append(list(_cyclic_key(boundary)))
    return sorted(boundaries), successor


def _permutation_orbits(mapping):
    _require(set(mapping.values()) == set(mapping), 'region action is not a permutation')
    remaining = set(mapping)
    orbits = []
    while remaining:
        first = min(remaining)
        orbit = []
        point = first
        while point in remaining:
            remaining.remove(point)
            orbit.append(point)
            point = mapping[point]
        _require(point == first, 'region orbit did not close')
        orbits.append(orbit)
    return orbits


def reconstruct_prong_certificate(graph, punctures):
    """Reconstruct regions from a final graph, not from a classification report.

    Input marking: type-1 edges on one specified peripheral boundary carry its
    puncture number.  This function verifies that those boundaries are distinct
    and exhaust exactly ``punctures`` boundaries.  The remaining boundary is
    distinguished as outer.  It cannot verify transport from an earlier input.

    The return value ALWAYS has ``pseudo_anosov_certified=False``.  Its prongs
    are cusp counts of the reconstructed track; they are foliation prongs only
    under a separate invariant-foliation realization theorem/certificate.
    """
    _require(_integer(punctures) and punctures >= 1, 'invalid puncture count')
    _require(isinstance(graph, dict), 'graph must be an object')
    raw_vertices, raw_edges = graph.get('vertices'), graph.get('edges')
    _require(isinstance(raw_vertices, list) and isinstance(raw_edges, list),
             'vertices and edges must be lists')
    vertices, edges = {}, {}
    for v in raw_vertices:
        _require(isinstance(v, dict) and _integer(v.get('label')), 'invalid vertex label')
        _require(v['label'] not in vertices, 'duplicate vertex label')
        vertices[v['label']] = v
    for e in raw_edges:
        _require(isinstance(e, dict) and _integer(e.get('label')) and e['label'] > 0,
                 'invalid edge label')
        _require(e['label'] not in edges, 'duplicate edge label')
        _require(_integer(e.get('start')) and e['start'] in vertices
                 and _integer(e.get('end')) and e['end'] in vertices, 'edge endpoint absent')
        _require(_integer(e.get('type')) and e['type'] in (0, 1, 2), 'invalid edge type')
        if e['type'] == 1:
            _require(_integer(e.get('puncture')) and 1 <= e['puncture'] <= punctures,
                     'invalid peripheral marking')
        edges[e['label']] = e
    _require(bool(vertices) and bool(edges), 'empty graph')

    def start(d):
        return edges[abs(d)]['start' if d > 0 else 'end']

    def word(d):
        w = edges[abs(d)]['image']
        return w if d > 0 else _inverse(w)

    rotations, origins = {}, {}
    for label, v in vertices.items():
        _require(_integer(v.get('image')) and v['image'] in vertices, 'vertex image absent')
        ds = v.get('cyclic_directions')
        _require(isinstance(ds, list) and bool(ds), 'invalid cyclic directions')
        for d in ds:
            _require(_integer(d) and d and abs(d) in edges and d not in origins
                     and start(d) == label, 'invalid rotation system')
            origins[d] = label
        rotations[label] = ds[:]
    _require(len(origins) == 2 * len(edges), 'rotation system misses a direction')
    adjacency = {v: set() for v in vertices}
    for e in edges.values():
        adjacency[e['start']].add(e['end'])
        adjacency[e['end']].add(e['start'])
    reached = {next(iter(vertices))}
    pending = list(reached)
    while pending:
        for v in adjacency[pending.pop()] - reached:
            reached.add(v)
            pending.append(v)
    _require(len(reached) == len(vertices), 'disconnected graph')
    original_boundaries, _ = _boundary_cycles(rotations)
    _require(len(original_boundaries) == punctures + 1
             and len(vertices) - len(edges) + len(original_boundaries) == 2,
             'ribbon graph is not a spine of the prescribed punctured sphere')

    derivative, turns = {}, set()
    for e in edges:
        w = edges[e].get('image')
        _require(isinstance(w, list) and bool(w)
                 and all(_integer(d) and d in origins for d in w), 'invalid or collapsed edge image')
        _require(start(w[0]) == vertices[start(e)]['image']
                 and start(-w[-1]) == vertices[start(-e)]['image'], 'edge image endpoints mismatch')
        for a, b in zip(w, w[1:]):
            _require(start(-a) == start(b), 'discontinuous edge image')
            turns.add(tuple(sorted((-a, b))))
        derivative[e], derivative[-e] = w[0], -w[-1]
    pending = list(turns)
    while pending:
        a, b = pending.pop()
        _require(a != b, 'an iterate backtracks')
        _require(origins[a] == origins[b], 'turn has different origins')
        image = tuple(sorted((derivative[a], derivative[b])))
        if image not in turns:
            turns.add(image)
            pending.append(image)

    gates, direction_gate, whitehead = {}, {}, {}
    for v, ds in rotations.items():
        parent = {d: d for d in ds}

        def root(a):
            while parent[a] != a:
                a = parent[a]
            return a

        for a, b in combinations(ds, 2):
            x, y, seen = a, b, set()
            while x != y and (x, y) not in seen:
                seen.add((x, y))
                x, y = derivative[x], derivative[y]
            if x == y:
                parent[root(a)] = root(b)
        roots = [root(d) for d in ds]
        _supported(len(set(roots)) >= 2, 'fewer than two gates at a vertex')
        cut = next(i for i in range(len(ds)) if roots[i] != roots[i - 1])
        ordered = ds[cut:] + ds[:cut]
        blocks = []
        for d in ordered:
            if not blocks or root(blocks[-1][0]) != root(d):
                blocks.append([])
            blocks[-1].append(d)
        _supported(len(blocks) == len(set(roots)), 'a gate is not a cyclic interval')
        gates[v] = blocks
        for i, block in enumerate(blocks):
            for d in block:
                direction_gate[d] = (v, i)
        used = set()
        for a, b in turns:
            if origins[a] == v:
                ga, gb = direction_gate[a][1], direction_gate[b][1]
                _require(ga != gb, 'used turn lies within a gate')
                used.add(tuple(sorted((ga, gb))))
        k = len(blocks)
        circular = {tuple(sorted((i, (i + 1) % k))) for i in range(k)}
        _supported(used <= circular, 'used-turn edge is not between adjacent gates')
        _supported(len(used) in (k - 1, k), 'used-turn graph is not a polygon or a polygon minus one edge')
        # For k=2 the unique supported model is one edge, never a digon.
        whitehead[v] = used

    # Place the old directions outwards; the new arcs lie inside the vertex
    # disk in its cyclic order.  The two tangent sides at a switch are exactly
    # "old" and "new", so a boundary turn is a cusp iff its sides agree.
    new_edges, new_rotations, gate_labels, corner_edge = [], {}, {}, {}
    next_edge = max(edges) + 1
    for v in sorted(vertices):
        for i in range(len(gates[v])):
            gate_labels[(v, i)] = len(gate_labels) + 1
        for i, j in sorted(whitehead[v]):
            label = next_edge
            next_edge += 1
            corner_edge[(v, i, j)], corner_edge[(v, j, i)] = label, -label
            new_edges.append({'label': label, 'start': gate_labels[(v, i)],
                              'end': gate_labels[(v, j)], 'kind': 'infinitesimal',
                              'original_vertex': v})
    for v in sorted(vertices):
        k = len(gates[v])
        for i, block in enumerate(gates[v]):
            ds = block[:]
            # k=2 has a single local arc; append it once.
            neighbors = [(i + 1) % k] if k == 2 else [(i + 1) % k, (i - 1) % k]
            for j in neighbors:
                if (v, i, j) in corner_edge:
                    ds.append(corner_edge[(v, i, j)])
            new_rotations[gate_labels[(v, i)]] = ds
    for e in sorted(edges):
        new_edges.append({'label': e, 'start': gate_labels[direction_gate[e]],
                          'end': gate_labels[direction_gate[-e]], 'kind': 'original'})
    boundaries, successor = _boundary_cycles(new_rotations)
    _require(len(new_rotations) - len(new_edges) + len(boundaries) == 2,
             'gate blow-up has nonzero genus')

    boundary_kind, marking = {}, {}
    for boundary in original_boundaries:
        key = _cyclic_key(boundary)
        labels = {edges[abs(d)].get('puncture') for d in boundary}
        if all(edges[abs(d)]['type'] == 1 for d in boundary) and len(labels) == 1:
            puncture = labels.pop()
            _require(puncture not in marking, 'puncture marks two boundary components')
            marking[puncture] = key
            boundary_kind[key] = ('puncture', puncture)
        else:
            boundary_kind[key] = ('outer', None)
    _require(set(marking) == set(range(1, punctures + 1)), 'peripheral markings do not exhaust the punctures')
    _require(sum(kind == 'outer' for kind, _ in boundary_kind.values()) == 1,
             'there must be exactly one unmarked outer boundary')

    new_edge_by_label = {e['label']: e for e in new_edges}
    regions, exterior_regions, interior_regions = [], {}, {}
    for boundary in boundaries:
        cusps = []
        for d in boundary:
            a, b = -d, successor[-d]
            if (abs(a) in edges) == (abs(b) in edges):
                cusps.append([a, b])
        projection = [d for d in boundary if abs(d) in edges]
        region = {'id': len(regions), 'boundary': boundary,
                  'projected_boundary': list(_cyclic_key(projection)),
                  'cusps': cusps, 'prongs': len(cusps)}
        if projection:
            key = _cyclic_key(projection)
            _require(key in boundary_kind and key not in exterior_regions,
                     'blown-up boundary does not project bijectively to original boundaries')
            kind, label = boundary_kind[key]
            region['kind'], region['puncture'] = kind, label
            exterior_regions[key] = region['id']
            _require(len(cusps) >= 1, 'marked complementary region has no cusp')
        else:
            origins_here = {new_edge_by_label[abs(d)]['original_vertex'] for d in boundary}
            _require(len(origins_here) == 1, 'interior region crosses original vertex disks')
            v = origins_here.pop()
            _require(v not in interior_regions and len(cusps) == len(gates[v]) and len(cusps) >= 3,
                     'invalid interior complementary polygon')
            region['kind'], region['original_vertex'] = 'interior', v
            interior_regions[v] = region['id']
        regions.append(region)
    _require(set(exterior_regions) == set(boundary_kind), 'an original boundary was lost')
    _require(sum(r['prongs'] - 2 for r in regions) == -4, 'sphere cusp Euler formula fails')

    action = {}
    for key, region_id in exterior_regions.items():
        image = _cyclic_key(_cyclic_reduce(d for e in key for d in word(e)))
        _require(image in exterior_regions, 'graph map does not preserve oriented boundary classes')
        action[region_id] = exterior_regions[image]
    for v, region_id in interior_regions.items():
        target = vertices[v]['image']
        _require(target in interior_regions, 'interior polygon maps outside the interior polygons')
        k = len(gates[v])
        image_gates = []
        for block in gates[v]:
            images = {direction_gate[derivative[d]] for d in block}
            _require(len(images) == 1, 'derivative does not descend to gates')
            image_gates.append(images.pop()[1])
        _require(k == len(gates[target]) and len(set(image_gates)) == k,
                 'interior polygon gate map is not bijective')
        _require(all(image_gates[(i + 1) % k] == (image_gates[i] + 1) % k for i in range(k)),
                 'interior polygon gate map reverses or scrambles cyclic order')
        action[region_id] = interior_regions[target]
    orbits = _permutation_orbits(action)
    for region in regions:
        image = regions[action[region['id']]]
        _require(region['kind'] == image['kind'] and region['prongs'] == image['prongs'],
                 'region action changes its marking type or cusp count')
        if region['kind'] == 'outer':
            _require(action[region['id']] == region['id'], 'outer marked point is not fixed')
        region['image'] = action[region['id']]
    puncture_prongs = {r['puncture']: r['prongs'] for r in regions if r['kind'] == 'puncture'}
    return {
        'schema': 1,
        'certificate_kind': 'gate_blowup_complementary_regions',
        'pseudo_anosov_certified': False,
        'input_to_final_provenance_certified': False,
        'invariant_foliation_realization_certified': False,
        'graph_endpoints_and_all_iterated_turns_checked': True,
        'marked_genus_zero_gate_blowup_checked': True,
        'uses_backend_type_or_prong_claims': False,
        'punctures': punctures,
        'gate_data': [{'vertex': v, 'gates': gates[v],
                       'used_gate_turns': [list(p) for p in sorted(whitehead[v])]}
                      for v in sorted(vertices)],
        'all_used_turns': [list(p) for p in sorted(turns)],
        'blown_up_graph': {'vertices': [{'label': label, 'cyclic_directions': ds}
                                        for label, ds in sorted(new_rotations.items())],
                          'edges': sorted(new_edges, key=lambda e: e['label'])},
        'regions': regions,
        'region_orbits': orbits,
        'reconstructed_interior_prongs': sorted(r['prongs'] for r in regions if r['kind'] == 'interior'),
        'reconstructed_outer_prongs': next(r['prongs'] for r in regions if r['kind'] == 'outer'),
        'reconstructed_puncture_prongs': [puncture_prongs[i] for i in range(1, punctures + 1)],
    }


def validate_reported_prongs(result):
    """Check every reported location, prong count, and ordered orbit.

    Successful return authenticates the report relative to the reconstructed
    cusped track only.  ``result['type']`` is neither inspected nor trusted.
    An inconsistent report raises TopologyCertificateError, NOT a mathematical
    target-mismatch exclusion.  The historic input is never mutated.
    """
    _require(isinstance(result, dict) and 'final_graph' in result, 'missing final graph')
    certificate = reconstruct_prong_certificate(result['final_graph'], result.get('punctures'))
    regions = certificate['regions']
    lookup = {}
    for r in regions:
        key = ('interior', r['original_vertex']) if r['kind'] == 'interior' else (
            'boundary', tuple(r['projected_boundary']))
        lookup[key] = r['id']
    seen = set()
    reports = result.get('singularity_orbits')
    _require(isinstance(reports, list), 'missing singularity report')
    for orbit in reports:
        _require(isinstance(orbit, dict) and isinstance(orbit.get('singularities'), list)
                 and bool(orbit['singularities']), 'empty or malformed singularity orbit')
        ids = []
        for s in orbit['singularities']:
            _require(isinstance(s, dict) and type(s.get('interior')) is bool,
                     'invalid singularity kind')
            location, prongs = s.get('location'), s.get('prongs')
            _require(isinstance(location, list) and bool(location)
                     and all(_integer(x) for x in location), 'invalid singularity location')
            _require(_integer(prongs) and prongs >= 1, 'invalid reported prong count')
            if s['interior']:
                _require(len(location) == 1, 'interior location must be one original vertex')
                key = ('interior', location[0])
            else:
                key = ('boundary', _cyclic_key(location))
            _require(key in lookup, 'reported singularity does not match a reconstructed region')
            region_id = lookup[key]
            _require(region_id not in seen, 'a complementary region is reported twice')
            _require(prongs == regions[region_id]['prongs'], 'reported prongs disagree with reconstructed cusp count')
            seen.add(region_id)
            ids.append(region_id)
        _require(all(regions[region_id]['image'] == ids[(i + 1) % len(ids)]
                     for i, region_id in enumerate(ids)), 'reported singularity orbit is incorrect')
    _require(seen == set(range(len(regions))), 'report omits reconstructed regions')
    certificate['reported_prongs_locations_and_orbits_match'] = True
    return certificate


def certify_pillowcase_affine(matrix, interior_marks, outer_mark):
    """A COMPLETE pseudo-Anosov certificate for a restricted explicit family.

    ``matrix`` is a hyperbolic integral 2x2 matrix of determinant +1.  It defines
    the map on R^2/Z^2 and commutes with x -> -x.  Its pillowcase quotient is a
    sphere with four one-prong branch points.  All four must be marked, the
    distinguished outer mark must be fixed, and the other marks must be
    permuted.  Extra marks are regular (two-prong) points.

    A rational point is [[x_numerator, x_denominator], [y_numerator,
    y_denominator]].  Coordinates are taken modulo Z^2 and sign.  Floating-point
    input is forbidden.  Marking a point by the same quotient point twice is an
    error.  This certifies the explicitly given quotient map, not its equality
    to a supplied braid or graph map (no such input is accepted).
    """
    _require(isinstance(matrix, list) and len(matrix) == 2
             and all(isinstance(row, list) and len(row) == 2 for row in matrix)
             and all(_integer(x) for row in matrix for x in row),
             'affine matrix must be a 2 by 2 integer matrix')
    a, b = matrix[0]
    c, d = matrix[1]
    _require(a * d - b * c == 1, 'affine matrix must have determinant +1')
    trace = a + d
    _require(abs(trace) > 2, 'affine matrix is not hyperbolic')

    def quotient_key(point):
        reduced = tuple(t % 1 for t in point)
        opposite = tuple((-t) % 1 for t in point)
        return min(reduced, opposite)

    def parse_point(point):
        _require(isinstance(point, list) and len(point) == 2, 'point must have two rational coordinates')
        coordinates = []
        for fraction in point:
            _require(isinstance(fraction, list) and len(fraction) == 2
                     and all(_integer(x) for x in fraction) and fraction[1] != 0,
                     'rational coordinates must be exact numerator-denominator pairs')
            coordinates.append(Fraction(*fraction))
        return quotient_key(coordinates)

    def serialize(point):
        return [[x.numerator, x.denominator] for x in point]

    def image(point):
        x, y = point
        return quotient_key((a * x + b * y, c * x + d * y))

    def prongs(point):
        return 1 if all((2 * t).denominator == 1 for t in point) else 2

    _require(isinstance(interior_marks, list), 'interior marks must be a list')
    interior = [parse_point(p) for p in interior_marks]
    outer = parse_point(outer_mark)
    points = interior + [outer]
    _require(len(set(points)) == len(points), 'duplicate marked quotient point')
    branch_points = {(Fraction(i, 2), Fraction(j, 2)) for i in (0, 1) for j in (0, 1)}
    _require(branch_points <= set(points), 'an unmarked one-prong branch point is forbidden')
    _require(image(outer) == outer, 'outer quotient point is not fixed')
    point_index = {p: i for i, p in enumerate(interior)}
    _require(all(image(p) in point_index for p in interior), 'interior marks are not invariant')
    permutation = {i: point_index[image(p)] for i, p in enumerate(interior)}
    orbits = _permutation_orbits(permutation)
    _require(sum(prongs(p) - 2 for p in points) == -4, 'affine sphere prong formula failed')
    return {
        'schema': 1,
        'certificate_kind': 'explicit_hyperbolic_torus_pillowcase_quotient',
        'pseudo_anosov_certified': True,
        'invariant_foliation_realization_certified': True,
        'uses_backend_type_or_prong_claims': False,
        'claims_equality_to_any_braid_or_graph_map': False,
        'matrix': [row[:] for row in matrix],
        'trace': trace,
        'determinant': 1,
        'dilatation_minimal_polynomial': [1, -abs(trace), 1],
        'dilatation_specification': 'the larger positive root, strictly greater than 1',
        'dilatation_strict_integer_bounds': [abs(trace) - 1, abs(trace)],
        'punctures': len(interior),
        'interior_marks': [serialize(p) for p in interior],
        'outer_mark': serialize(outer),
        'puncture_permutation': [permutation[i] for i in range(len(interior))],
        'puncture_orbits': orbits,
        'interior_unpunctured_prongs': [],
        'puncture_prongs': [prongs(p) for p in interior],
        'outer_prongs': prongs(outer),
        'surface_genus': 0,
    }
