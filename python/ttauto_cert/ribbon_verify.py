"""Exact ribbon automorphisms, fold inverses and surface-boundary checks.

This audits exported finite graph maps, not completeness of pA normal forms.
"""
from collections import deque


def _require(condition, *context):
    if not condition:
        raise ValueError("Ribbon certificate verification failed" + (": " + repr(context) if context else ""))


def inverse(word):
    return [-x for x in reversed(word)]


def reduce_word(word):
    result = []
    for d in word:
        if result and result[-1] == -d:
            result.pop()
        else:
            result.append(d)
    return result


def cyclic_word(word):
    result = reduce_word(word)
    while len(result) >= 2 and result[0] == -result[-1]:
        result = result[1:-1]
    return min(tuple(result[k:]+result[:k]) for k in range(len(result))) if result else ()


def substitute(words, word):
    return [d for e in word for d in (words[e-1] if e > 0 else inverse(words[-e-1]))]


def boundary_cycles(graph):
    previous = {}
    for directions in graph['cyclic_directions']:
        for i, d in enumerate(directions):
            previous[d] = directions[i-1]
    todo, cycles = set(previous), []
    while todo:
        d, word = min(todo), []
        while d in todo:
            todo.remove(d)
            word.append(d)
            d = previous[-d]
        _require(d == word[0])
        cycles.append(word)
    return cycles


def automorphisms(graph):
    rotation, colors = {}, {}
    for directions in graph['cyclic_directions']:
        for i, d in enumerate(directions):
            rotation[d] = directions[(i+1) % len(directions)]
    for p, polygon in enumerate(graph['polygons']):
        for d in polygon:
            colors[d] = (graph['punctured'][p], graph['labels'][p])
    result = []
    m = graph['main_edges']
    for target in [sign*e for e in range(1, m+1) for sign in (1, -1)]:
        mapping, reverse, queue, valid = {}, {}, deque([(1, target)]), True
        while queue and valid:
            a, b = queue.popleft()
            if a in mapping:
                valid = mapping[a] == b
                continue
            if b in reverse or (abs(a) <= m) != (abs(b) <= m):
                valid = False
                break
            if a > m and (b <= m or colors[a] != colors[b]):
                valid = False
                break
            mapping[a], reverse[b] = b, a
            queue.extend([(-a, -b), (rotation[a], rotation[b])])
        if valid and len(mapping) == len(rotation):
            result.append([mapping[e] for e in range(1, len(graph['endpoints'])+1)])
    return result


def inverse_elementary(words):
    back = [None for _ in words]
    long_words = []
    for e, word in enumerate(words, 1):
        if len(word) == 1:
            d = word[0]
            _require(back[abs(d)-1] is None)
            back[abs(d)-1] = [e if d > 0 else -e]
        else:
            long_words.append((e, word))
    _require(len(long_words) == 1)
    e, word = long_words[0]
    missing = [i for i, d in enumerate(word) if back[abs(d)-1] is None]
    _require(len(missing) == 1)
    k = missing[0]
    prefix = substitute(back, word[:k])
    suffix = substitute(back, word[k+1:])
    candidate = inverse(prefix)+[e]+inverse(suffix)
    back[abs(word[k])-1] = candidate if word[k] > 0 else inverse(candidate)
    _require(all(w is not None for w in back))
    for e in range(1, len(words)+1):
        _require(reduce_word(substitute(back, words[e-1])) == [e])
        _require(reduce_word(substitute(words, back[e-1])) == [e])
    return back



def automorphisms_ordered(graph):
    """Same complete group as C++, identity first then vertex/edge lex order."""
    origins = {d: v for v, ds in enumerate(graph['cyclic_directions']) for d in ds}
    records = []
    for permutation in automorphisms(graph):
        def image(d):
            return permutation[d-1] if d > 0 else -permutation[-d-1]
        vertices = [origins[image(ds[0])] for ds in graph['cyclic_directions']]
        edges = [[d] for d in permutation]
        m = graph['main_edges']
        matrix = [[0]*m for _ in range(m)]
        for j in range(m):
            matrix[abs(permutation[j])-1][j] = 1
        records.append({'vertex_images': vertices, 'edge_images': edges, 'matrix': matrix})
    records.sort(key=lambda a: (a['vertex_images'], a['edge_images']))
    identity = [a for a in records if a['edge_images'] == [[e+1] for e in range(len(graph['endpoints']))]]
    _require(len(identity) == 1, 'identity automorphism missing')
    return identity + [a for a in records if a is not identity[0]]


def validate_graph(graph):
    endpoints = graph['endpoints']
    directions = graph['cyclic_directions']
    m = graph['main_edges']
    _require(0 <= m <= len(endpoints), 'main edge count')
    origin = {}
    for v, star in enumerate(directions):
        _require(bool(star), 'isolated vertex')
        for d in star:
            _require(type(d) is int and d != 0 and abs(d) <= len(endpoints), 'invalid dart')
            _require(d not in origin, 'repeated dart')
            origin[d] = v
            _require(endpoints[abs(d)-1][0 if d > 0 else 1] == v, 'dart endpoint')
    _require(len(origin) == 2*len(endpoints), 'missing dart')
    _require(len(graph['polygons']) == len(graph['punctured']) == len(graph['labels']), 'polygon colors')
    cycles = boundary_cycles(graph)
    _require(len(cycles) == len(graph['polygons'])+1, 'not a polygon tree surface')
    _require(len(directions)-len(endpoints)+len(cycles) == 2, 'not genus zero')
    actual = {cyclic_word(w) for w in cycles}
    for polygon in graph['polygons']:
        _require(cyclic_word(polygon) in actual, 'polygon is not an oriented face')
    return origin, cycles


def verify_graph(data, require_exported_automorphisms=True):
    """Verify one exported automaton and its complete terminal symmetry groups.

    Raises ValueError on mathematical check failure, including under python -O.
    Does not certify a complete train-track normal-form theorem or pA status.
    """
    vertices = data['vertices']
    _require([v['id'] for v in vertices] == list(range(len(vertices))), 'vertex indices')
    prepared = [validate_graph(v['ribbon_graph']) for v in vertices]
    nontrivial, arcs, group_checks = [], 0, 0
    for vertex in vertices:
        graph = vertex['ribbon_graph']
        group = automorphisms_ordered(graph)
        if 'automorphisms' in vertex:
            exported = [{key: a[key] for key in ('vertex_images', 'edge_images', 'matrix')}
                        for a in vertex['automorphisms']]
            _require(exported == group, 'incomplete or misordered terminal automorphisms', vertex['id'])
            group_checks += 1
        elif require_exported_automorphisms:
            raise ValueError('Export missing terminal automorphisms at vertex ' + str(vertex['id']))
        if len(group) > 1:
            nontrivial.append({'vertex': vertex['id'], 'order': len(group), 'automorphisms': group})
        source_origin, source_boundaries = prepared[vertex['id']]
        for arc in vertex['arcs']:
            words, vertex_images = arc['edge_images'], arc['vertex_images']
            target = vertices[arc['target']]['ribbon_graph']
            target_origin, target_boundaries = prepared[arc['target']]
            _require(len(words) == len(graph['endpoints']), 'edge map size')
            _require(sorted(vertex_images) == list(range(len(target['cyclic_directions']))), 'vertex map is not bijective')
            for e, word in enumerate(words, 1):
                _require(bool(word), 'empty edge image')
                _require(target_origin[word[0]] == vertex_images[source_origin[e]], 'image start')
                _require(target_origin[-word[-1]] == vertex_images[source_origin[-e]], 'image end')
                for a, b in zip(word, word[1:]):
                    _require(target_origin[-a] == target_origin[b], 'discontinuous edge image')
                # Each expanded switch has main darts on one tangent side
                # and its two infinitesimal darts on the other. These finite
                # checks therefore prove that smooth paths remain smooth
                # under every iterate, not only that endpoints match once.
                if e <= graph['main_edges']:
                    _require(abs(word[0]) <= target['main_edges']
                             and abs(word[-1]) <= target['main_edges'], 'main endpoint smoothing')
                else:
                    _require(len(word) == 1 and abs(word[0]) > target['main_edges'],
                             'infinitesimal edge transport')
                for a, b in zip(word, word[1:]):
                    _require((abs(a) <= target['main_edges']) != (abs(b) <= target['main_edges']),
                             'illegal internal train-track turn')
            inverse_elementary(words)
            images = sorted(cyclic_word(substitute(words, w)) for w in source_boundaries)
            boundaries = sorted(cyclic_word(w) for w in target_boundaries)
            _require(images == boundaries, 'oriented surface boundaries not preserved', vertex['id'], arc['branch'])
            m = graph['main_edges']
            matrix = [[0]*m for _ in range(m)]
            for j in range(m):
                for d in words[j]:
                    if abs(d) <= m:
                        matrix[abs(d)-1][j] += 1
            _require(matrix == arc['matrix'], 'main matrix mismatch')
            arcs += 1
    return {
        'n': data['n'], 'stratum_index': data['stratum_index'],
        'vertices': len(vertices), 'arcs': arcs,
        'all_arcs_have_explicit_graph_homotopy_inverses': True,
        'all_arcs_preserve_oriented_surface_boundary_cycles': True,
        'all_arcs_respect_train_track_smoothing': True,
        'exported_terminal_groups_checked': group_checks,
        'nontrivial_orientation_preserving_automorphisms': nontrivial,
        'pseudo_anosov_or_global_completeness_certified': False,
    }
