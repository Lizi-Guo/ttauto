# SPDX-License-Identifier: GPL-3.0-or-later
"""Import old integer spectra as data, without inheriting historical claims."""
from pathlib import Path
import json
from .storage import atomic_json, file_digest


def import_catalog(path, output, field=None):
    path, output = Path(path).resolve(), Path(output).resolve()
    if output.exists():
        raise ValueError('catalogue output already exists')
    if path.suffix.lower() == '.json':
        data = json.loads(path.read_text(encoding='utf-8'))
        if field is not None:
            for key in field.split('.'):
                data = data[key]
        elif isinstance(data, dict):
            data = data['polynomials']
    else:
        if field is not None:
            raise ValueError('--field applies to JSON only')
        data = [list(map(int, line.split())) for line in path.read_text(encoding='utf-8').splitlines() if line.strip()]
    if not isinstance(data, list) or any(not isinstance(p, list) or len(p) < 3 or len(p) % 2 == 0 or
            p[0] != 1 or p != p[::-1] or any(type(c) is not int for c in p) for p in data):
        raise ValueError('expected a list of monic integral reciprocal coefficient lists')
    rows = sorted(set(map(tuple, data)), key=lambda p: (len(p), p))
    result = {'schema': 'ttauto-polynomial-catalog-v1', 'polynomials': rows,
              'provenance': {'input_name': path.name, 'input_sha256': file_digest(path), 'field': field},
              'duplicates_removed': len(data)-len(rows), 'completeness_proved': False,
              'historical_minimum_or_classifier_claims_imported': False}
    output.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(output, result)
    return {'output': str(output), 'polynomials': len(rows), 'completeness_proved': False}
