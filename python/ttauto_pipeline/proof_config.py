# SPDX-License-Identifier: GPL-3.0-or-later
"""Proof specifications do not require a train track or an automaton."""
from copy import deepcopy
from fractions import Fraction
from .config import keys, positive
from ttauto_cert.exact_math import count_real_roots, compare_largest_real_root

MODULES = ('trace', 'schur', 'perron', 'mod2', 'peripheral', 'lefschetz',
           'deck', 'base_coupling', 'cycle')


def normalize(data):
    cfg = deepcopy(data)
    keys(cfg, {'schema', 'name', 'surface', 'strata', 'upper', 'target_polynomial',
               'source', 'pruning', 'limits', 'witness'}, 'proof configuration')
    surface = cfg.get('surface')
    keys(surface, {'kind', 'genus', 'punctures', 'one_pronged'}, 'proof surface')
    if surface.get('kind') == 'disc':
        positive(surface.get('punctures'), 'punctures', 3)
        if 'genus' in surface and surface['genus'] != 0:
            raise ValueError('disc base genus must be zero')
        surface['genus'] = 0
        if type(surface.setdefault('one_pronged', False)) is not bool:
            raise ValueError('one_pronged must be boolean')
    elif surface.get('kind') == 'closed':
        positive(surface.get('genus'), 'genus', 2)
        if 'punctures' in surface or 'one_pronged' in surface:
            raise ValueError('closed proof has no puncture options')
    else:
        raise ValueError('proof surface.kind must be disc or closed')
    strata = cfg.setdefault('strata', 'all')
    if strata != 'all' and (not isinstance(strata, list) or not strata):
        raise ValueError('strata must be all or a nonempty explicit list')
    if isinstance(cfg.get('upper'), (bool, float)):
        raise ValueError('upper must be an exact rational')
    upper = Fraction(cfg['upper'])
    if upper <= 1:
        raise ValueError('upper must exceed one')
    cfg['upper'] = str(upper)
    target = cfg.get('target_polynomial')
    if target is not None:
        if (not isinstance(target, list) or len(target) < 2 or target[0] != 1 or
                any(type(c) is not int for c in target) or not count_real_roots(target) or
                compare_largest_real_root(target, 1) <= 0 or compare_largest_real_root(target, upper) >= 0):
            raise ValueError('target must have a largest real root strictly between 1 and upper')
    source = cfg.setdefault('source', {'kind': 'newton'})
    keys(source, {'kind', 'path', 'executable', 'batch_nodes', 'trace_horizon'}, 'proof source')
    if source.get('kind') not in ('newton', 'newton_native', 'polynomial_catalog', 'residual'):
        raise ValueError('proof source must be newton, newton_native, polynomial_catalog or residual')
    if source['kind'] in ('polynomial_catalog', 'residual') and not isinstance(source.get('path'), str):
        raise ValueError('catalogue or residual source requires a path')
    if source['kind'] not in ('polynomial_catalog', 'residual') and 'path' in source:
        raise ValueError('newton source takes no path')
    if source['kind'] == 'newton_native':
        if not isinstance(source.get('executable'), str):
            raise ValueError('newton_native requires the local trace-prefix executable')
        positive(source.setdefault('batch_nodes', 10000), 'batch_nodes')
        positive(source.setdefault('trace_horizon', 32), 'native trace_horizon')
    elif 'executable' in source or 'batch_nodes' in source or 'trace_horizon' in source:
        raise ValueError('native settings require newton_native')
    pruning = cfg.setdefault('pruning', {})
    keys(pruning, {*MODULES, 'trace_horizon'}, 'proof pruning')
    for name in MODULES:
        if type(pruning.setdefault(name, name != 'cycle')) is not bool:
            raise ValueError(f'pruning.{name} must be boolean')
    positive(pruning.setdefault('trace_horizon', 32), 'trace_horizon', 2)
    limits = cfg.setdefault('limits', {})
    keys(limits, {'events', 'models', 'cycle_nodes', 'base_catalog_nodes'}, 'proof limits')
    for name, value in limits.items():
        positive(value, f'limits.{name}', 0)
    witness = cfg.get('witness')
    if witness is not None:
        keys(witness, {'kind', 'path', 'start', 'folds', 'terminal'}, 'witness')
        if witness.get('kind') != 'fold_path' or not isinstance(witness.get('path'), str):
            raise ValueError('only a local fold_path witness is supported; external classifiers are forbidden')
        for field in ('start', 'terminal'):
            positive(witness.get(field), f'witness.{field}', 0)
        if not isinstance(witness.get('folds'), list) or not witness['folds']:
            raise ValueError('witness requires a nonempty folds list')
        for value in witness['folds']:
            positive(value, 'fold index', 0)
    return cfg
