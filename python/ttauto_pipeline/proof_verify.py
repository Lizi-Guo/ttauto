# SPDX-License-Identifier: GPL-3.0-or-later
"""Replay spectral decisions and optionally the full Newton frontier."""
from copy import deepcopy
from pathlib import Path
import json
import tempfile
from ttauto_cert.spectral import trace_bounds
from ttauto_cert.lefschetz import traces
from .config import normalize
from .proof_backend import ProofBackend
from .storage import Journal, digest, source_digest, atomic_json


def verify(folder, replay_search=False):
    from .engine import run
    folder = Path(folder).resolve()
    config = normalize(json.loads((folder/'config.json').read_text(encoding='utf-8')))
    saved = json.loads((folder/'result.json').read_text(encoding='utf-8'))
    backend = ProofBackend(config, folder, True, source_digest(), snapshot_only=True)
    journal = Journal(folder/'run.sqlite3', backend.identity, resume=True)
    try:
        checkpoint = journal.get('checkpoint')
        source = backend.source(checkpoint['cursor'])
        expected = backend.report(source, checkpoint, checkpoint['status'], checkpoint.get('error'))
        if digest(saved) != digest(expected):
            raise ValueError('proof report or mathematical claims differ from committed evidence')
        counts = {'events': checkpoint['counts']['events'], 'prefixes': checkpoint['counts']['prefixes'], 'candidates': 0}
        if 'native_statistics' in checkpoint['counts']:
            counts['native_statistics'] = checkpoint['counts']['native_statistics']
        catalog = {tuple(p) for p in backend.space.get('polynomial_catalog', [])}
        for number, record in journal.records():
            if number != counts['candidates']+1:
                raise ValueError('candidate journal has a gap or duplicate sequence')
            p = record.get('polynomial')
            if not isinstance(p, list) or not p or p[0] != 1 or p != p[::-1] or any(type(x) is not int for x in p):
                raise ValueError('invalid recorded spectrum')
            candidate = source.candidate(p)
            if not candidate['strata'] or any(record.get(k) != v for k, v in candidate.items()):
                raise ValueError('candidate type, rank or stratum association changed')
            if config['source']['kind'] == 'polynomial_catalog':
                if tuple(p) not in catalog:
                    raise ValueError('candidate not in supplied catalogue')
            else:
                half = (len(p)-1)//2
                if any(abs(t) > b for t, b in zip(traces(p, half)[1:], trace_bounds(len(p)-1, config['upper'], half))):
                    raise ValueError('candidate outside Newton containing bounds')
            counts['candidates'] += 1
            actual = backend.evaluate(candidate, counts)
            if digest(actual) != digest(record):
                raise ValueError(f'candidate {number}: proof decision or evidence replay differs')
        if counts != checkpoint['counts'] or counts['events'] != counts['prefixes']+counts['candidates']:
            raise ValueError('proof stage accounting differs')
        if checkpoint['best'] is not None:
            raise ValueError('spectral survivors cannot be recorded as attained maps')
        if replay_search:
            cfg = deepcopy(config)
            with tempfile.TemporaryDirectory() as temporary:
                temp = Path(temporary)
                if cfg['source']['kind'] == 'polynomial_catalog':
                    atomic_json(temp/'catalog.json', backend.space['polynomial_catalog'])
                    cfg['source']['path'] = str(temp/'catalog.json')
                elif cfg['source']['kind'] == 'residual':
                    cfg['source']['path'] = str(folder/'parent_evidence')
                if cfg.get('witness'):
                    cfg['witness']['path'] = str(folder/'witness_automaton.json')
                replay = run(cfg, temp/'run')
                if digest(replay) != digest(saved):
                    raise ValueError('full spectral search replay differs from saved proof')
        return {'verified': True, 'schema': saved['schema'], 'records': counts['candidates'],
                'search_replayed': replay_search, 'candidate_superset_complete': saved['candidate_superset_complete'],
                'lower_bound_proved': saved['lower_bound_proved'], 'minimum_proved': saved['minimum_proved'],
                'automaton_required': False}
    finally:
        journal.close()
