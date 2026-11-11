# SPDX-License-Identifier: GPL-3.0-or-later
"""Continue verified necessary-data survivors, preserving the parent proof."""
from pathlib import Path
import json
import shutil
import sqlite3
from .storage import digest


def snapshot(source, destination):
    """Copy only proof evidence, so a child inside the parent cannot recurse."""
    source, destination = Path(source).resolve(), Path(destination).resolve()
    if source == destination or destination.exists():
        raise ValueError('parent evidence destination must be new')
    destination.mkdir()
    for name in ('config.json', 'result.json', 'candidate_space.json', 'witness_automaton.json'):
        if (source/name).is_file():
            shutil.copy2(source/name, destination/name)
    src = sqlite3.connect(f'{(source/"run.sqlite3").as_uri()}?mode=ro', uri=True)
    dst = sqlite3.connect(destination/'run.sqlite3')
    try:
        src.backup(dst)
    finally:
        dst.close(); src.close()
    if (source/'parent_evidence').is_dir():
        snapshot(source/'parent_evidence', destination/'parent_evidence')


def residual_space(config, space, folder, resume, snapshot_only):
    from .proof_verify import verify
    parent = Path(folder)/'parent_evidence'
    if not resume:
        snapshot(config['source']['path'], parent)
    verified = verify(parent)
    report = json.loads((parent/'result.json').read_text(encoding='utf-8'))
    parent_cfg = json.loads((parent/'config.json').read_text(encoding='utf-8'))
    parent_space = json.loads((parent/'candidate_space.json').read_text(encoding='utf-8'))
    for key in ('surface', 'strata', 'upper', 'target_polynomial'):
        if config.get(key) != parent_cfg.get(key):
            raise ValueError(f'residual proof must preserve parent {key}')
    if space['cases'] != parent_space['cases']:
        raise ValueError('parent stratum inventory differs')
    if resume and not snapshot_only:
        original = json.loads((Path(config['source']['path'])/'result.json').read_text(encoding='utf-8'))
        if digest(original) != digest(report):
            raise ValueError('resume refused: parent proof changed')
    indices = {case['id']: i for i, case in enumerate(space['cases'])}
    pending = []
    connection = sqlite3.connect(parent/'run.sqlite3')
    try:
        for payload, in connection.execute('SELECT payload FROM candidates ORDER BY id'):
            row = json.loads(payload)
            if row['status'] == 'unresolved':
                targets = [indices[c['stratum']] for c in row['case_results'] if c['status'] == 'unresolved']
                pending.append({'polynomial': row['polynomial'], 'strata': targets})
    finally:
        connection.close()
    space['residual_candidates'] = pending
    space['residual_complete'] = verified['candidate_superset_complete']
    space['parent_report_sha256'] = digest(report)
    space['coverage']['parent_exclusions_verified'] = True
    space['coverage']['parent_superset_complete'] = space['residual_complete']
    return space
