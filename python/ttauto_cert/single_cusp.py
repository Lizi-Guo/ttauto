"""Bind the single-outer-cusp factorization theorem to exported fold words.

The mathematical theorem is in review/n6_s3_move_coverage_zh.md.  This audit
checks its two directed slides at EVERY supplied state; it does not establish
that the supplied states exhaust a stratum, or classify a map as pA.
"""
from .ribbon_verify import cyclic_word, inverse, substitute, verify_graph
from .standard_maps import StandardState


def _require(condition, message):
    if not condition:
        raise ValueError('Single-cusp move audit: '+message)


def single_cusp_directions(graph):
    """Return both ordered (long a, short b, following inf j0) choices.

    Signs are directed GRAPH edges, not letters after collapsing a spanning
    tree.  The reversed outer cycle supplies the opposite slide direction.
    """
    state = StandardState(graph)
    real_degrees = [sum(abs(d) <= state.m for d in star)
                    for star in graph['cyclic_directions']]
    _require(sum(d-1 for d in real_degrees) == 1,
             'requires exactly one outer cusp')
    choices = []
    for direction, cycle in ((1, state.outer), (-1, inverse(state.outer))):
        rr = [i for i in range(len(cycle))
              if abs(cycle[i]) <= state.m and
              abs(cycle[(i+1) % len(cycle)]) <= state.m]
        _require(len(rr) == 1, 'outer boundary must have exactly one R-R turn')
        _require(not any(abs(cycle[i]) > state.m and
                         abs(cycle[(i+1) % len(cycle)]) > state.m
                         for i in range(len(cycle))),
                 'outer boundary cannot have an I-I turn')
        i = rr[0]
        a, b, j = -cycle[i], cycle[(i+1) % len(cycle)], cycle[(i+2) % len(cycle)]
        _require(state.origins[a] == state.origins[b], 'cusp endpoints')
        _require(abs(j) > state.m, 'edge after short edge must be infinitesimal')
        _require(state.origins[-b] == state.origins[j], 'short/inf endpoints')
        _require(state.origins[-a] != state.origins[-b], 'parallel real edges')
        _require(real_degrees[state.origins[-b]] == 1,
                 'short edge must end at a single-real switch')
        choices.append({'outer_orientation': direction, 'long_edge': a,
                        'short_edge': b, 'infinitesimal_edge': j})
    _require(choices[0]['long_edge'] == choices[1]['short_edge'] and
             choices[0]['short_edge'] == choices[1]['long_edge'],
             'reversed outer cycle must exchange the two directions')
    return choices


def audit_single_cusp_moves(data):
    """Validate complete paired-slide coverage on a supplied exported graph.

    Includes the existing independent inverse, boundary, matrix and complete
    terminal-group audit.  Raises ValueError on missing, duplicate or wrong
    moves, including under python -O.  No floating-point checks are used.
    """
    vertices = data['vertices']
    _require(bool(vertices), 'empty state set')
    ribbon_result = verify_graph(data)
    states = [StandardState(v['ribbon_graph']) for v in vertices]
    records = []
    arc_count = 0
    for vertex, state in zip(vertices, states):
        directions = single_cusp_directions(state.graph)
        matches = [[] for _ in directions]
        for arc_index, arc in enumerate(vertex['arcs']):
            words = arc['edge_images']
            target = states[arc['target']]
            _require(target.m == state.m, 'real dimension changed')
            # Keep the intended zero/true-puncture colors as well as all
            # positive peripheral cycles, even though every polygon is a
            # puncture for the free-group part of the theory.
            target_cycles = {cyclic_word(p): i for i, p in enumerate(target.graph['polygons'])}
            for p, cycle in enumerate(state.graph['polygons']):
                image = cyclic_word(substitute(words, cycle))
                _require(image in target_cycles, 'positive infinitesimal polygon transport')
                q = target_cycles[image]
                _require((state.graph['punctured'][p], state.graph['labels'][p]) ==
                         (target.graph['punctured'][q], target.graph['labels'][q]),
                         'polygon color changed')
            found = []
            for k, choice in enumerate(directions):
                a, b, j = (choice[key] for key in
                           ('long_edge', 'short_edge', 'infinitesimal_edge'))
                qa, qb, qj = (substitute(words, [d]) for d in (a, b, j))
                if (len(qa) == 3 and len(qb) == len(qj) == 1 and
                        qa[:2] == qb+qj and abs(qa[2]) <= target.m and
                        all(len(words[e-1]) == 1 for e in range(1, state.m+1)
                            if e != abs(a))):
                    record = dict(choice, arc_index=arc_index,
                                  branch=arc.get('branch'), raw_fold=arc.get('raw_fold'),
                                  target=arc['target'], target_real_suffix=qa[2],
                                  long_image=qa, short_image=qb, inf_image=qj)
                    matches[k].append(record)
                    found.append(k)
            _require(len(found) == 1,
                     'arc must realize exactly one prescribed paired slide at state '+str(vertex['id']))
            arc_count += 1
        _require(all(len(group) == 1 for group in matches),
                 'both directed slides must occur exactly once at state '+str(vertex['id']))
        records.append({'state': vertex['id'], 'outer_cycle': state.outer,
                        'directions': [group[0] for group in matches]})
    return {'schema': 'ttauto-single-cusp-move-audit-v1', 'validated': True,
            'state_count': len(vertices), 'arc_count': arc_count,
            'state_audits': records,
            'ribbon_verification': ribbon_result,
            'all_stratum_states_certified': False,
            'pseudo_anosov_certified': False}
