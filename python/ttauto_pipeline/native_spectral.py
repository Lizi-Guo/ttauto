# SPDX-License-Identifier: GPL-3.0-or-later
"""Optional batched native Newton source with durable unvisited prefixes."""
from copy import deepcopy
from functools import lru_cache
from pathlib import Path
import json
import subprocess
import tempfile
from ttauto_cert.lefschetz import partitions, traces
from ttauto_cert.hyperelliptic_mod2 import branch_permutation_charpoly, _bits
from ttauto_cert.spectral import trace_bounds
from .spectral_source import SpectralSource


@lru_cache(maxsize=24)
def parity_masks(rank):
    return tuple(sorted({_bits(branch_permutation_charpoly(part)) for part in partitions(rank+2)}))


class NativeSpectralSource(SpectralSource):
    def __init__(self, space, config, cursor=None):
        super().__init__(space, config)
        if any(rank > 48 for rank in self.ranks):
            raise ValueError('native source supports reciprocal ranks up to 48; use newton for arbitrary integers')
        self.cursor = deepcopy(cursor) if cursor is not None else {
            'frontier': [[rank, []] for rank in reversed(self.ranks)], 'pending': [], 'pending_index': 0}

    def next_event(self):
        c = self.cursor
        if c['pending_index'] < len(c['pending']):
            p = c['pending'][c['pending_index']]
            c['pending_index'] += 1
            return self.candidate(p)
        c['pending'], c['pending_index'] = [], 0
        if not c['frontier']:
            return None
        rank, prefix = c['frontier'].pop()
        half = rank//2
        h = max(half, self.config['source']['trace_horizon']) if self.config['pruning']['trace'] else half
        bounds = trace_bounds(rank, self.config['upper'], h)
        masks = parity_masks(rank) if self.config['surface']['kind'] == 'disc' and self.config['pruning']['mod2'] else ()
        budget = self.config['source']['batch_nodes']
        with tempfile.TemporaryDirectory(prefix='ttauto-newton-') as temporary:
            root = Path(temporary)
            content = (f'{half} {h} {len(masks)}\n'+' '.join(map(str, bounds))+'\n'+
                       ' '.join(map(str, masks))+'\n'+str(len(prefix))+' '+
                       ' '.join(map(str, prefix))+f'\n{budget}\n')
            (root/'input').write_text(content, encoding='ascii')
            process = subprocess.run([self.config['source']['executable'], str(root/'input'),
                                      str(root/'polynomials'), str(root/'report'), str(root/'frontier')],
                                     capture_output=True, text=True)
            if process.returncode:
                raise RuntimeError('native Newton source failed (not an exclusion): '+process.stderr.strip())
            stats = json.loads((root/'report').read_text())
            rows = [list(map(int, line.split())) for line in (root/'polynomials').read_text().splitlines()]
            frontier = [list(map(int, line.split())) for line in (root/'frontier').read_text().splitlines()]
        if (stats['retained'] != len(rows) or stats['deferred_prefixes'] != len(frontier) or
                stats['completed'] != (not frontier) or
                stats['parity_admissible_completions'] != stats['trace_pruned']+len(rows)):
            raise ValueError('native Newton accounting mismatch')
        for p in rows:
            if (len(p) != rank+1 or p[0] != 1 or p != p[::-1] or p[1:len(prefix)+1] != prefix or
                    any(abs(t) > b for t, b in zip(traces(p, h)[1:], bounds)) or
                    masks and _bits(p) not in masks):
                raise ValueError('invalid native reciprocal spectrum')
        if len({tuple(p) for p in rows}) != len(rows):
            raise ValueError('duplicate native polynomial')
        for row in frontier:
            if row[0] != len(row)-1 or row[1:len(prefix)+1] != prefix:
                raise ValueError('invalid deferred native prefix')
        c['frontier'].extend([rank, row[1:]] for row in reversed(frontier))
        c['pending'] = rows
        return {'kind': 'prefix', 'candidate_kind': 'native_newton_batch', 'rank': rank,
                'coefficient_prefix': prefix, 'native_statistics': stats,
                'early_trace_filter': self.config['pruning']['trace'],
                'early_uncolored_mod2_filter': bool(masks)}
