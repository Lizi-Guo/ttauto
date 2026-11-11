"""Bind all guarded paired slides to full signed words, for any outer cusp count.

The combination of the boundary-defect and real-forest arguments is proved
in review/multicusp_general_prefix_20260915_zh.md. This module audits the
local move inventory and terminal groups of a supplied automaton. It does
not certify that all stratum states were supplied, or classify any map.
"""
from .ribbon_verify import cyclic_word, inverse, substitute, verify_graph
from .standard_maps import StandardState


def _require(condition, message):
    if not condition:
        raise ValueError('Paired-fold audit: '+message)


def paired_slide_directions(graph):
    """Find the exact short-edge guard on each side of every exterior cusp.

    Returning fewer than twice the cusp count is expected: a prospective
    short edge may end at a non-extreme position in its other real fan.
    """
    state = StandardState(graph)
    p = sum(sum(abs(d) <= state.m for d in star)-1
            for star in graph['cyclic_directions'])
    _require(p > 0, 'positive outer cusp count required')
    choices = []
    for orientation, cycle in ((1, state.outer), (-1, inverse(state.outer))):
        positions = [i for i in range(len(cycle))
                     if abs(cycle[i]) <= state.m and
                     abs(cycle[(i+1) % len(cycle)]) <= state.m]
        _require(len(positions) == p, 'outer cusp count mismatch')
        _require(not any(abs(cycle[i]) > state.m and
                         abs(cycle[(i+1) % len(cycle)]) > state.m
                         for i in range(len(cycle))), 'exterior I-I turn')
        for i in positions:
            a, b, j = (-cycle[i], cycle[(i+1) % len(cycle)],
                       cycle[(i+2) % len(cycle)])
            _require(state.origins[a] == state.origins[b], 'cusp endpoints')
            if abs(j) <= state.m:
                continue
            _require(state.origins[-b] == state.origins[j], 'short/inf endpoints')
            _require(state.origins[-a] != state.origins[-b], 'parallel real edges')
            choices.append({'outer_orientation': orientation,
                            'long_edge': a, 'short_edge': b,
                            'infinitesimal_edge': j})
    keys = [(c['long_edge'], c['short_edge'], c['infinitesimal_edge']) for c in choices]
    _require(len(set(keys)) == len(keys), 'duplicate geometric slide')
    return choices


def audit_paired_folds(data):
    """Reject missing, duplicate or incorrect slides, including under -O.

    Labels and raw fold numbers are not trusted as descriptions of a move:
    matching is by its complete signed edge words after normalization.
    """
    vertices = data['vertices']
    _require(bool(vertices), 'empty state set')
    for vertex in vertices:
        for arc in vertex['arcs']:
            _require(type(arc['target']) is int and 0 <= arc['target'] < len(vertices),
                     'invalid target state')
    ribbon = verify_graph(data)
    states = [StandardState(v['ribbon_graph']) for v in vertices]
    m = states[0].m
    _require(all(s.m == m for s in states), 'real dimension changed')
    records = []
    for vertex, state in zip(vertices, states):
        directions = paired_slide_directions(state.graph)
        matches = [[] for _ in directions]
        for arc_index, arc in enumerate(vertex['arcs']):
            words = arc['edge_images']
            target = states[arc['target']]
            target_cycles = {cyclic_word(p): i
                             for i, p in enumerate(target.graph['polygons'])}
            for p, cycle in enumerate(state.graph['polygons']):
                image = cyclic_word(substitute(words, cycle))
                _require(image in target_cycles, 'positive polygon transport')
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
                        qa[:2] == qb+qj and abs(qa[2]) <= m and
                        all(len(words[e-1]) == 1 for e in range(1, m+1)
                            if e != abs(a))):
                    matches[k].append(dict(choice, arc_index=arc_index,
                                           branch=arc.get('branch'), target=arc['target'],
                                           long_image=qa, short_image=qb, inf_image=qj))
                    found.append(k)
            _require(len(found) == 1, 'arc must realize exactly one guarded slide')
        _require(all(len(group) == 1 for group in matches),
                 'each guarded slide must occur exactly once at state '+str(vertex['id']))
        records.append({'state': vertex['id'], 'outer_cycle': state.outer,
                        'outer_cusps': sum(sum(abs(d) <= m for d in star)-1
                                          for star in state.graph['cyclic_directions']),
                        'directions': [group[0] for group in matches]})
    return {'schema': 'ttauto-paired-fold-move-audit-v1', 'validated': True,
            'state_count': len(states), 'arc_count': sum(len(v['arcs']) for v in vertices),
            'state_audits': records, 'ribbon_verification': ribbon,
            'all_stratum_states_certified': False, 'pseudo_anosov_certified': False}
