# SPDX-License-Identifier: GPL-3.0-or-later
from collections import Counter

def need(condition, message):
    if not condition:
        raise ValueError(message)


def verify(source, target, arc, orders):
    """Use chord separation, without the generator's adjacent-pair reduction."""
    words = arc['edge_images']
    vm = arc['vertex_images']
    ss, ts = source['cyclic_directions'], target['cyclic_directions']
    e, m = len(words), source['main_edges']
    need(len(ss) == len(ts) and sorted(vm) == list(range(len(ts))), 'vertex permutation')
    need(target['main_edges'] == m and len(orders) == e, 'dimensions')
    darts = {d for k in range(1, e + 1) for d in (k, -k)}
    for stars in (ss, ts):
        need(Counter(d for star in stars for d in star) == Counter(darts), 'dart partition')
    so = {d: v for v, star in enumerate(ss) for d in star}
    to = {d: v for v, star in enumerate(ts) for d in star}
    uses = [[] for _ in range(e)]
    for k, word in enumerate(words, 1):
        need(word and all(d in darts for d in word), 'nonempty edge path')
        need(to[word[0]] == vm[so[k]] and to[-word[-1]] == vm[so[-k]], 'endpoints')
        need(all(to[-a] == to[b] for a, b in zip(word, word[1:])), 'continuity')
        need((abs(word[0]) <= m) == (k <= m) == (abs(word[-1]) <= m), 'end types')
        need(all((abs(a) <= m) != (abs(b) <= m) for a, b in zip(word, word[1:])), 'turn types')
        for j, d in enumerate(word):
            uses[abs(d)-1].append((k, j))
    for supplied, expected in zip(orders, uses):
        need(Counter(map(tuple, supplied)) == Counter(expected), 'lane partition')
    for v, star in enumerate(ts):
        ports = []
        for d in star:
            lanes = orders[abs(d)-1]
            for k, j in (lanes if d > 0 else reversed(lanes)):
                word = words[k-1]
                if word[j] == d:
                    ports.append(('H', k) if j == 0 else ('T', k, j-1))
                else:
                    ports.append(('H', -k) if j+1 == len(word) else ('T', k, j))
        hubs = [(i, p[1]) for i, p in enumerate(ports) if p[0] == 'H']
        source_star = ss[vm.index(v)]
        labels = [d for _, d in hubs]
        need(len(labels) == len(source_star), 'hub degree')
        need(any(labels == source_star[i:] + source_star[:i]
                 for i in range(len(source_star))), 'hub cyclic order')
        turns = {}
        for i, p in enumerate(ports):
            if p[0] == 'T':
                turns.setdefault(p, []).append(i)
        need(all(len(positions) == 2 for positions in turns.values()), 'turn pair')
        chords = list(turns.values())
        for a, b in chords:
            inside_hubs = sum(a < h < b for h, _ in hubs)
            need(inside_hubs in (0, len(hubs)), 'chord separates source hub')
        for i, (a, b) in enumerate(chords):
            for c, d in chords[i+1:]:
                need((a < c < b) == (a < d < b), 'crossing chords')
    return True

