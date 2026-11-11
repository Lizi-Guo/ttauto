"""Exact greedy paired-fold factorization of a complete standard graph map.

The general prefix theorem proves termination under the four conditions
checked here. This is a map factorizer, never a pseudo-Anosov classifier.
An audited automaton need not contain every state in a stratum: its complete
outgoing guarded moves suffice to follow this particular finite factorization.
"""
from copy import deepcopy

from .paired_folds import audit_paired_folds
from .ribbon_verify import cyclic_word, inverse_elementary, reduce_word, substitute
from .standard_maps import StandardState, surjective_free_group


def _require(condition, message):
    if not condition:
        raise ValueError('Paired factorization: '+message)


def _real_length(words, m):
    return sum(abs(d) <= m for word in words[:m] for d in word)


def validate_standard_map(source, target, vertex_images, edge_images):
    """Check all four hypotheses, including integral free-group surjectivity.

    ``source`` and ``target`` are validated StandardState instances. The free
    groups have equal finite rank, so surjectivity implies isomorphism.
    """
    _require((source.m, source.e, source.q) == (target.m, target.e, target.q),
             'source and target standard dimensions differ')
    v = len(source.graph['cyclic_directions'])
    _require(len(vertex_images) == v and
             all(type(i) is int for i in vertex_images) and
             sorted(vertex_images) == list(range(len(target.graph['cyclic_directions']))),
             'vertex images must be a bijection')
    _require(len(edge_images) == source.e, 'wrong number of edge images')
    for e, word in enumerate(edge_images, 1):
        _require(bool(word) and all(type(d) is int and d in target.origins for d in word),
                 'empty edge image or unknown directed edge')
        _require(target.origins[word[0]] == vertex_images[source.origins[e]] and
                 target.origins[-word[-1]] == vertex_images[source.origins[-e]],
                 'edge image endpoints disagree with vertex map')
        _require(all(target.origins[-a] == target.origins[b] for a, b in zip(word, word[1:])),
                 'discontinuous edge image')
        if e <= source.m:
            _require(len(word) % 2 == 1 and
                     all((abs(d) <= target.m) == (i % 2 == 0) for i, d in enumerate(word)),
                     'real edge image must be R(IR)*')
        else:
            _require(len(word) == 1 and abs(word[0]) > target.m,
                     'infinitesimal image must be one infinitesimal edge')
    target_cycles = {cyclic_word(p): i for i, p in enumerate(target.graph['polygons'])}
    polygon_images = []
    for p, cycle in enumerate(source.graph['polygons']):
        image = cyclic_word(substitute(edge_images, cycle))
        _require(image in target_cycles, 'oriented polygon transport failed')
        q = target_cycles[image]
        _require(len(cycle) == len(target.graph['polygons'][q]) and
                 (source.graph['punctured'][p], source.graph['labels'][p]) ==
                 (target.graph['punctured'][q], target.graph['labels'][q]),
                 'polygon type or color changed')
        polygon_images.append(q)
    _require(sorted(polygon_images) == list(range(target.q)), 'polygon transport is not bijective')
    images = [target.project(substitute(edge_images, word)) for word in source.basis]
    _require(surjective_free_group(images, target.q), 'induced free-group map is not an isomorphism')
    _require(cyclic_word(substitute(edge_images, source.outer)) == cyclic_word(target.outer),
             'oriented outer boundary is not preserved')
    return {'four_conditions_verified': True, 'polygon_images': polygon_images,
            'free_group_images': images, 'pseudo_anosov_certified': False}


def _terminal_isomorphism(source, target, vertices, words):
    _require(all(len(word) == 1 for word in words), 'terminal map has a nonsingleton image')
    _require(sorted(abs(word[0]) for word in words) == list(range(1, target.e+1)),
             'terminal edge map is not bijective')
    for v, star in enumerate(source.graph['cyclic_directions']):
        image = substitute(words, star)
        other = target.graph['cyclic_directions'][vertices[v]]
        _require(len(image) == len(other) and
                 any(image == other[k:]+other[:k] for k in range(len(other))),
                 'terminal map does not preserve ribbon rotation')


class PairedFactorizer:
    """Audit and snapshot an automaton once; factor multiple complete maps."""
    def __init__(self, data):
        self._data = deepcopy(data)
        self.audit = audit_paired_folds(self._data)
        self._states = [StandardState(v['ribbon_graph']) for v in self._data['vertices']]
        self._choices = deepcopy([row['directions'] for row in self.audit['state_audits']])

    def factor(self, source, target, vertex_images, edge_images):
        _require(all(type(i) is int and 0 <= i < len(self._states) for i in (source, target)),
                 'source/target state index is invalid')
        original_vertices = list(vertex_images)
        original_words = [list(w) for w in edge_images]
        vertices, words = list(original_vertices), deepcopy(original_words)
        current, m = source, self._states[source].m
        validation = validate_standard_map(self._states[source], self._states[target], vertices, words)
        initial = _real_length(words, m)
        steps = []
        length = initial
        while length > m:
            choice = None
            for candidate in self._choices[current]:
                a, b, j = (candidate[k] for k in ('long_edge', 'short_edge', 'infinitesimal_edge'))
                fa, fb, fj = (substitute(words, [d]) for d in (a, b, j))
                if len(fa) > len(fb) and fa[:len(fb)] == fb:
                    _require(fa[len(fb):len(fb)+1] == fj,
                             'whole-prefix direction contradicts the boundary theorem')
                    choice = candidate
                    break
            if choice is None:
                raise RuntimeError('Paired factorization: nonterminal four-condition map has no '
                                   'audited whole-prefix factor; mathematical/implementation inconsistency')
            arc_index = choice['arc_index']
            arc = self._data['vertices'][current]['arcs'][arc_index]
            middle = arc['target']
            back = inverse_elementary(arc['edge_images'])
            residual = [reduce_word(substitute(words, w)) for w in back]
            transported = [None]*len(vertices)
            for old, new in enumerate(arc['vertex_images']):
                transported[new] = vertices[old]
            _require([substitute(residual, w) for w in arc['edge_images']] == words,
                     'F=GQ failed as a literal edge-word identity')
            _require([transported[v] for v in arc['vertex_images']] == vertices,
                     'F=GQ failed on vertices')
            validate_standard_map(self._states[middle], self._states[target], transported, residual)
            following_length = _real_length(residual, m)
            short_length = sum(abs(d) <= m for d in substitute(words, [choice['short_edge']]))
            _require(following_length == length-short_length and m <= following_length < length,
                     'paired factor did not decrease the proved length potential')
            steps.append({'source': current, 'arc_index': arc_index, 'target': middle,
                          'long_edge': choice['long_edge'], 'short_edge': choice['short_edge'],
                          'infinitesimal_edge': choice['infinitesimal_edge'],
                          'real_length_before': length, 'real_length_after': following_length})
            current, words, vertices, length = middle, residual, transported, following_length
        _require(len(steps) <= initial-m, 'proved step bound was exceeded')
        _terminal_isomorphism(self._states[current], self._states[target], vertices, words)
        terminal_index = None
        if current == target:
            matches = [i for i, a in enumerate(self._data['vertices'][current]['automorphisms'])
                       if a['vertex_images'] == vertices and a['edge_images'] == words]
            _require(len(matches) == 1, 'terminal isomorphism missing from audited complete group')
            terminal_index = matches[0]
        terminal = {'source': current, 'target': target, 'vertex_images': vertices,
                    'edge_images': words, 'automorphism_index': terminal_index,
                    'ribbon_isomorphism_verified': True}
        # Independently compose the returned factors forward, without any
        # free reduction, and include the actual terminal isomorphism.
        replay_words = [[i+1] for i in range(self._states[source].e)]
        replay_vertices = list(range(len(original_vertices)))
        for step in steps:
            arc = self._data['vertices'][step['source']]['arcs'][step['arc_index']]
            replay_words = [substitute(arc['edge_images'], w) for w in replay_words]
            replay_vertices = [arc['vertex_images'][v] for v in replay_vertices]
        replay_words = [substitute(words, w) for w in replay_words]
        replay_vertices = [vertices[v] for v in replay_vertices]
        _require(replay_words == original_words and replay_vertices == original_vertices,
                 'returned full factorization failed forward literal replay')
        return {'schema': 'ttauto-greedy-paired-factorization-v1', 'completed': True,
                'source': source, 'target': target, 'dimension': m,
                'initial_real_length': initial, 'proved_step_bound': initial-m,
                'step_count': len(steps), 'steps': steps, 'terminal': terminal,
                'input_validation': validation, 'all_residual_four_conditions_verified': True,
                'literal_forward_composition_verified': True, 'backtracking_used': False,
                'pseudo_anosov_certified': False, 'all_stratum_states_certified': False}


def factor_paired(data, source, target, vertex_images, edge_images):
    """One-shot interface. Reuse PairedFactorizer for a batch on one graph."""
    return PairedFactorizer(data).factor(source, target, vertex_images, edge_images)
