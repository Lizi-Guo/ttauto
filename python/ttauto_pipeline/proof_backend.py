# SPDX-License-Identifier: GPL-3.0-or-later
"""Spectrum backend for the same transactional engine used by folding paths."""
from copy import deepcopy
from pathlib import Path
import json
from ttauto_cert.exact_math import compare_largest_real_roots
from .proof_inventory import inventory
from .spectral_source import SpectralSource
from .proof_pruning import evaluate
from .storage import atomic_json, digest, file_digest


class ProofBackend:
    def __init__(self, config, folder, resume, sources, snapshot_only=False):
        self.config, self.folder = config, Path(folder)
        self.space = inventory(config)
        if config['source']['kind'] == 'residual':
            from .residual import residual_space
            self.space = residual_space(config, self.space, folder, resume, snapshot_only)
        if config['source']['kind'] == 'polynomial_catalog':
            if snapshot_only:
                rows = json.loads((self.folder/'candidate_space.json').read_text(encoding='utf-8'))['polynomial_catalog']
            else:
                packet = json.loads(Path(config['source']['path']).read_text(encoding='utf-8'))
                rows = packet.get('polynomials') if isinstance(packet, dict) else packet
            ranks = {r['rank'] for r in self.space['cases']}
            if not isinstance(rows, list) or any(not isinstance(p, list) or len(p)-1 not in ranks or
                    p[0] != 1 or p != p[::-1] or any(type(x) is not int for x in p) for p in rows):
                raise ValueError('catalogue must contain monic integral reciprocal polynomials of the declared ranks')
            if len({tuple(p) for p in rows}) != len(rows):
                raise ValueError('duplicate catalogue polynomials')
            self.space['polynomial_catalog'] = rows
            self.space['coverage']['catalogue_completeness_proved'] = False
        path = self.folder/'candidate_space.json'
        if resume:
            if not path.is_file() or json.loads(path.read_text(encoding='utf-8')) != self.space:
                raise ValueError('resume refused: candidate space or source catalogue changed')
        else:
            atomic_json(path, self.space)
        self.witness = self.check_witness(resume, snapshot_only)
        bound = deepcopy(config)
        bound['limits'].pop('events', None)
        if 'path' in bound['source']:
            bound['source']['path'] = 'content-bound-catalogue'
        if 'executable' in bound['source']:
            bound['source']['executable'] = file_digest(bound['source']['executable'])
        if bound.get('witness'):
            bound['witness']['path'] = 'content-bound-witness'
        self.identity = {'schema': 'ttauto-proof-identity-v1', 'config': bound,
                         'candidate_space_sha256': digest(self.space),
                         'witness_sha256': digest(self.witness), 'python_sources_sha256': sources}

    def source(self, cursor=None):
        if self.config['source']['kind'] == 'newton_native':
            from .native_spectral import NativeSpectralSource
            return NativeSpectralSource(self.space, self.config, cursor)
        return SpectralSource(self.space, self.config, cursor)

    def record_prefix(self, event, counts):
        if 'native_statistics' in event:
            stats = counts.setdefault('native_statistics', {})
            for key in ('nodes', 'parity_admissible_completions', 'trace_pruned', 'retained'):
                stats[key] = stats.get(key, 0)+event['native_statistics'][key]

    def validate_runtime(self):
        source = self.config['source']
        if 'executable' in source and file_digest(source['executable']) != self.identity['config']['source']['executable']:
            raise ValueError('native spectral backend changed during the run')

    def evaluate(self, candidate, counts):
        return evaluate(self.space, self.config, candidate, counts)

    def check_witness(self, resume, snapshot_only=False):
        spec = self.config.get('witness')
        if spec is None:
            return None
        # This constructs ONLY the explicitly supplied witness graph. A proof
        # without a witness never imports or constructs an Automaton.
        from .automaton import Automaton
        from .config import normalize
        from .engine import evaluate as evaluate_map
        from .ribbon import real_matrix
        snapshot = self.folder/'witness_automaton.json'
        data = json.loads((snapshot if snapshot_only else Path(spec['path'])).read_text(encoding='utf-8'))
        if resume:
            if not snapshot.is_file() or json.loads(snapshot.read_text(encoding='utf-8')) != data:
                raise ValueError('resume refused: witness graph changed')
        else:
            atomic_json(snapshot, data)
        a = Automaton(data, data['surface'])
        cfg = normalize(dict(schema='ttauto-experiment-v1', surface=data['surface'],
                             upper=self.config['upper'], geometry={'kind': 'builtin'}))
        words, _ = a.compose(spec['start'], spec['folds'], spec['terminal'])
        event = dict(start=spec['start'], path=spec['folds'], terminal=spec['terminal'],
                     matrix=real_matrix(words, a.dimension))
        record = evaluate_map(a, cfg, event, [], {})
        requested = self.config['surface']
        if requested['kind'] != cfg['surface']['kind'] or requested['genus'] != cfg['surface']['genus']:
            raise ValueError('witness belongs to another surface')
        if requested['kind'] == 'disc' and requested['punctures'] != cfg['surface']['punctures']:
            raise ValueError('witness puncture count differs')
        in_scope = False
        if record['status'] == 'attained':
            geo = record['geometry']
            if requested['kind'] == 'closed':
                kind = 'H' if geo['actual_foliation_orientable'] else 'Q'
                prongs = geo['actual_singular_prongs']
                orders = sorted((p-2)//2 if kind == 'H' else p-2 for p in prongs)
                in_scope = any(case['type'] == kind and case['orders'] == orders for case in self.space['cases'])
            else:
                in_scope = any(case['n'] == requested['punctures'] and
                               case['interior_prongs'] == geo['actual_interior_prongs'] and
                               case['outer_prongs'] == geo['actual_outer_prongs'] for case in self.space['cases'])
        matched = bool(in_scope and self.config.get('target_polynomial') and
                       compare_largest_real_roots(record['polynomial'], self.config['target_polynomial']) == 0)
        return {'graph_sha256': digest(data), 'record': record,
                'in_requested_scope': in_scope, 'matches_target': matched,
                'external_classifier_used': False}

    def report(self, source, checkpoint, status, error=None):
        counts = checkpoint['counts']
        exhausted = status == 'exhausted'
        complete = exhausted and source.catalog_complete
        unresolved = counts.get('unresolved_stratum_candidates', 0)
        lower = complete and unresolved == 0
        attained = bool(self.witness and self.witness['record']['status'] == 'attained' and
                        self.witness['in_requested_scope'])
        report = {'schema': 'ttauto-proof-result-v1', 'name': self.config.get('name', 'proof'),
                  'status': status, 'enumeration_completed': exhausted, 'candidate_superset_complete': complete,
                  'surface': self.config['surface'], 'upper': self.config['upper'],
                  'target_polynomial': self.config.get('target_polynomial'),
                  'lower_bound': {'polynomial': self.config['target_polynomial'], 'root': 'largest_real'}
                                 if self.config.get('target_polynomial') else {'rational': self.config['upper']},
                  'coverage': self.space['coverage'], 'strata': len(self.space['cases']),
                  'candidate_kind': 'homological_spectrum', 'automaton_required': False,
                  'counts': counts, 'unresolved_candidates': unresolved,
                  'lower_bound_proved': lower, 'upper_bound_attained': attained,
                  'minimum_proved': bool(lower and attained and self.witness['matches_target']),
                  'best_attained': self.witness['record'] if attained else None,
                  'witness': self.witness, 'run_identity': self.identity,
                  'checkpoint_sha256': digest(checkpoint),
                  'proof_scope': 'strict lower bound in the stated candidate-space coverage, conditional on recorded mathematical theorem dependencies; no realization inferred from surviving necessary data'}
        if error:
            report['error'] = error
        return report
