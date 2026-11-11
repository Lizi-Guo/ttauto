# SPDX-License-Identifier: GPL-3.0-or-later
"""Optional proof modules on spectra, strata, and orbit data, before maps exist."""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from ttauto_cert.exact_math import (count_real_roots, compare_largest_real_root,
                                    compare_largest_real_roots, isolate_largest_real_root, _gcd, _trim)
from ttauto_cert.lefschetz import compatible, traces
from ttauto_cert.hyperelliptic_mod2 import branch_cycle_types, _bits
from ttauto_cert.perron_obstructions import strict_perron_obstruction
from ttauto_cert.spectral import (schur_inside, trace_bounds, opposite,
                                  peripheral_obstruction, reciprocal_polynomials, multiply)
from ttauto_cert.cover_orbits import (Budget, ModelBudgetExceeded, disc_patterns,
    disc_deck_feasible, joint_models, regular_from_anti, coupled_orbits)
from ttauto_cert.cycle_obstruction import signed_peripheral_polynomials, cycle_search
from .pruning import Decision
from .storage import digest


CONTRACTS = {
    'trace': 'reciprocal spectral trace bounds through the configured horizon',
    'schur': 'all roots strictly inside the rational enumeration envelope',
    'perron': 'positive lift or anti-invariant action: simple strictly dominant stretching root',
    'mod2': 'sphere cover only: colored branch permutation characteristic polynomial',
    'peripheral': 'one-pronged disc standard-track weight/homology exact sequence',
    'lefschetz': 'singular indices and nonnegative integral periodic orbit counts',
    'deck': 'sphere cover: both lifts, singular deck action, regular pairing and optional joint mod2',
    'base_coupling': 'closed Q: correlated singular data and invariant/anti-invariant homology',
    'cycle': 'disc residual: all signed peripheral factors and abstract cycle intersection connectivity',
}


def answer(excluded, reason, evidence=None):
    return Decision('reject' if excluded else 'pass', reason, evidence)


@lru_cache(maxsize=16)
def base_catalogue(rank, upper, limit):
    rows = []
    for index, p in enumerate(reciprocal_polynomials(rank, upper)):
        if limit is not None and index >= limit:
            raise ModelBudgetExceeded('base catalogue incomplete')
        if schur_inside(p, upper):
            rows.append(tuple(p))
    return tuple(rows)


class SpectrumContext:
    def __init__(self, config, polynomial, case):
        self.config, self.p, self.case = config, polynomial, case
        self.horizon = config['pruning']['trace_horizon']
        self.budget = Budget(config['limits'].get('models'))
        self._disc_patterns = None

    @property
    def disc_patterns(self):
        if self._disc_patterns is None:
            self._disc_patterns = disc_patterns(self.p, self.case, self.horizon, self.budget)
        return self._disc_patterns


def domain(p, config):
    if not count_real_roots(p) or compare_largest_real_root(p, 1) <= 0:
        return answer(True, 'no_positive_expanding_root')
    target = config.get('target_polynomial')
    cmp = compare_largest_real_roots(p, target) if target else compare_largest_real_root(p, config['upper'])
    return answer(cmp >= 0, 'strictly_below_target' if target else 'strictly_below_rational_target')


def algebraic(name, p, config):
    horizon = config['pruning']['trace_horizon']
    if name == 'trace':
        t = traces(p, horizon)[1:]
        bounds = trace_bounds(len(p)-1, config['upper'], horizon)
        failure = next((dict(period=k, trace=v, bound=b) for k, (v, b) in enumerate(zip(t, bounds), 1) if abs(v) > b), None)
        return answer(failure is not None, 'reciprocal_trace_envelope', failure)
    if name == 'schur':
        return answer(not schur_inside(p, config['upper']), 'exact_Schur_disk', {'radius': config['upper']})
    if name == 'perron':
        negative = opposite(p)
        if count_real_roots(negative) and compare_largest_real_roots(negative, p) >= 0:
            return answer(True, 'negative_eigenvalue_at_least_as_large')
        result = strict_perron_obstruction(p)
        return answer(result['excluded'], result['reason'], result)
    raise ValueError('unknown algebraic proof module')


def check(name, context):
    c = context
    p, case, horizon = c.p, c.case, c.horizon
    try:
        if name == 'mod2':
            if case['kind'] != 'disc':
                return Decision('skip', 'base_is_not_a_sphere')
            cv = case['cover']
            for row in branch_cycle_types(cv['branch_colors'], cv['fixed_branch_points']):
                c.budget.tick()
                if _bits(row['mod2_polynomial']) == _bits(p):
                    return answer(False, 'colored_branch_mod2', row)
            return answer(True, 'no_colored_branch_mod2_model')
        if name == 'peripheral':
            if case['kind'] != 'disc':
                return Decision('skip', 'disc_weight_homology_contract_only')
            e = sum(x % 2 == 0 for x in case['interior_prongs'])+int(case['outer_prongs'] % 2 == 0)
            failure = peripheral_obstruction(p, e, horizon)
            return answer(failure is not None, 'signed_peripheral_trace_orbits', failure)
        if name == 'lefschetz':
            if case['kind'] == 'disc':
                a, b = c.disc_patterns
                return answer(not a or not b, 'sphere_two_lift_orbits',
                              {'positive_patterns': len(a), 'negative_patterns': len(b), 'horizon': horizon})
            if case['type'] == 'H':
                groups = sorted(Counter(2*d+2 for d in case['orders']).items())
                sizes = []
                for negative in (False, True):
                    rows = compatible(p, groups, negative, stop=horizon, budget=c.budget)
                    sizes.append(len(rows))
                    if rows:
                        return answer(False, 'orientable_sign_and_orbits', {'negative': negative, 'model': rows[0]})
                return answer(True, 'neither_foliation_sign_has_orbits', {'pattern_counts': sizes})
            anti = traces(p, horizon)[1:]
            tested = 0
            for joint in joint_models(case['orders'], horizon, c.budget):
                tested += 1
                regular, failure = regular_from_anti(anti, joint)
                if failure is None:
                    return answer(False, 'anti_invariant_regular_orbits', {'joint_vector': joint, 'regular_orbits': regular[1]})
            return answer(True, 'no_joint_singular_orbit_model', {'models_tested': tested})
        if name == 'deck':
            if case['kind'] != 'disc':
                return Decision('skip', 'use_base_coupling_for_closed_Q')
            a, b = c.disc_patterns
            model = disc_deck_feasible(p, case, a, b, c.budget, c.config['pruning']['mod2'])
            return answer(model is None, 'sphere_deck_orbit_coupling', {'model': model, 'horizon': horizon})
        if name == 'base_coupling':
            if case['kind'] != 'closed' or case['type'] != 'Q':
                return Decision('skip', 'closed_nonorientable_foliation_only')
            bases = base_catalogue(2*case['genus'], c.config['upper'], c.config['limits'].get('base_catalog_nodes'))
            # Base eigenvalues have modulus <lambda. An upper endpoint gives a
            # conservative test, so no actual base action can be omitted.
            lo, hi = isolate_largest_real_root(p, bits=48)
            eligible = []
            for q in bases:
                if not schur_inside(q, hi):
                    continue
                shared = False
                for signed in (q, opposite(q)):
                    common = _gcd(_trim(p), _trim(signed))
                    if len(common) > 1 and count_real_roots(common) and compare_largest_real_roots(common, p) == 0:
                        shared = True
                        break
                if not shared:
                    eligible.append(q)
            anti = traces(p, horizon)[1:]
            base_traces = [(q, traces(q, horizon)[1:]) for q in eligible]
            tested = 0
            for joint in joint_models(case['orders'], horizon, c.budget):
                regular, failure = regular_from_anti(anti, joint)
                if failure is not None:
                    continue
                for q, base in base_traces:
                    c.budget.tick(); tested += 1
                    model, failure = coupled_orbits(anti, base, joint, regular)
                    if failure is None:
                        return answer(False, 'invariant_anti_invariant_coupling',
                                      {'base_polynomial': q, 'joint_vector': joint, 'model': model})
            return answer(True, 'all_base_and_joint_models_excluded',
                          {'complete_base_catalogue': len(bases), 'base_catalogue_sha256': digest(bases),
                           'eligible_base_spectra': len(eligible), 'coupled_pairs_tested': tested, 'horizon': horizon})
        if name == 'cycle':
            if case['kind'] != 'disc':
                return Decision('skip', 'no_closed_anti_invariant_to_real_matrix_identity_assumed')
            results, seen = [], set()
            for factor in signed_peripheral_polynomials(case, c.budget):
                key = tuple(factor)
                if key in seen:
                    continue
                seen.add(key)
                result = cycle_search(multiply(p, factor), c.config['limits'].get('cycle_nodes', 200000))
                results.append({'peripheral_factor': factor, 'search': result})
                if not result['excluded']:
                    return Decision('skip' if not result['completed'] and result.get('survivor') is None else 'pass',
                                    result['reason'], {'peripheral_searches': results})
            return answer(True, 'all_peripheral_cycle_graphs_excluded', {'peripheral_searches': results})
        raise ValueError('unknown stratum proof module')
    except ModelBudgetExceeded as error:
        return Decision('skip', 'resource_limit_unresolved', {'detail': str(error), 'models_consumed': c.budget.used})


def evaluate(space, config, candidate, counts):
    from .engine import stage
    record = {**candidate, 'decisions': [], 'case_results': []}
    def decide(name, decision, destination, scale):
        stage(counts, scale+name, decision.outcome)
        destination.append({'module': name, 'outcome': decision.outcome, 'reason': decision.reason,
                            'evidence': decision.evidence})
        return decision.outcome == 'reject'
    p = candidate['polynomial']
    # Cheap exact trace/disk filters can precede root isolation. They are
    # optional; the strict target interval defines the question itself.
    for name in ('trace', 'schur'):
        if config['pruning'][name] and decide(name, algebraic(name, p, config), record['decisions'], 'polynomial_'):
            return {**record, 'status': 'rejected', 'rejected_by': name}
    if decide('domain', domain(p, config), record['decisions'], 'polynomial_'):
        return {**record, 'status': 'rejected', 'rejected_by': 'domain'}
    if config['pruning']['perron'] and decide('perron', algebraic('perron', p, config), record['decisions'], 'polynomial_'):
        return {**record, 'status': 'rejected', 'rejected_by': 'perron'}
    for index in candidate['strata']:
        case = space['cases'][index]
        row = {'stratum': case['id'], 'decisions': [], 'status': 'unresolved'}
        context = SpectrumContext(config, p, case)
        counts['stratum_candidates'] = counts.get('stratum_candidates', 0)+1
        for name in ('mod2', 'peripheral', 'lefschetz', 'deck', 'base_coupling', 'cycle'):
            if config['pruning'][name] and decide(name, check(name, context), row['decisions'], 'stratum_'):
                row.update(status='rejected', rejected_by=name)
                break
        record['case_results'].append(row)
    unresolved = sum(row['status'] == 'unresolved' for row in record['case_results'])
    counts['unresolved_stratum_candidates'] = counts.get('unresolved_stratum_candidates', 0)+unresolved
    record['status'] = 'unresolved' if unresolved else 'rejected'
    return record
