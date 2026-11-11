# SPDX-License-Identifier: GPL-3.0-or-later
"""Resumable Newton enumeration; no graph construction precedes this source."""
from copy import deepcopy
from ttauto_cert.spectral import trace_bounds, coefficient_range


class SpectralSource:
    def __init__(self, space, config, cursor=None):
        self.space, self.config = space, config
        self.ranks = sorted({r['rank'] for r in space['cases']})
        self.cursor = deepcopy(cursor) if cursor is not None else {'rank_index': 0, 'stack': [], 'catalog_index': 0}
        self.depth_truncated = False
        self.catalog_complete = space.get('residual_complete', 'polynomial_catalog' not in space)

    def reject_prefix(self):
        self.cursor['stack'].pop()

    def next_event(self):
        c = self.cursor
        if 'polynomial_catalog' in self.space or 'residual_candidates' in self.space:
            rows = self.space.get('polynomial_catalog', self.space.get('residual_candidates'))
            if c['catalog_index'] == len(rows):
                return None
            row = rows[c['catalog_index']]
            p = row['polynomial'] if isinstance(row, dict) else row
            c['catalog_index'] += 1
            return self.candidate(p)
        while True:
            if not c['stack']:
                if c['rank_index'] == len(self.ranks):
                    return None
                rank = self.ranks[c['rank_index']]
                c['rank_index'] += 1
                c['stack'].append(self.node(rank, [1], [rank]))
            node = c['stack'][-1]
            rank, coeff = node['rank'], node['coeff']
            if len(coeff) == rank//2+1:
                c['stack'].pop()
                return self.candidate(coeff+coeff[-2::-1])
            if node['next'] > node['hi']:
                c['stack'].pop()
                continue
            value = node['next']; node['next'] += 1
            k = len(coeff)
            powers = node['powers']+[node['known']-k*value]
            child = self.node(rank, coeff+[value], powers)
            c['stack'].append(child)
            return {'kind': 'prefix', 'candidate_kind': 'reciprocal_coefficient_prefix',
                    'rank': rank, 'coefficients': child['coeff'], 'traces': powers[1:]}

    def node(self, rank, coeff, powers):
        node = dict(rank=rank, coeff=coeff, powers=powers)
        if len(coeff) <= rank//2:
            bounds = trace_bounds(rank, self.config['upper'], rank//2)
            lo, hi, known = coefficient_range(coeff, powers, bounds)
            node.update(next=lo, hi=hi, known=known)
        return node

    def candidate(self, p):
        rank = len(p)-1
        if 'residual_candidates' in self.space:
            matched = next((row for row in self.space['residual_candidates'] if row['polynomial'] == p), None)
            if matched is None:
                raise ValueError('spectrum absent from verified parent survivors')
            return {'kind': 'candidate', 'candidate_kind': 'homological_spectrum',
                    'rank': rank, 'polynomial': p, 'strata': matched['strata']}
        return {'kind': 'candidate', 'candidate_kind': 'homological_spectrum',
                'rank': rank, 'polynomial': p,
                'strata': [i for i, case in enumerate(self.space['cases']) if case['rank'] == rank]}
