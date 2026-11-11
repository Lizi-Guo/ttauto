# SPDX-License-Identifier: GPL-3.0-or-later
"""Proof coverage, historical arithmetic controls and non-Hall certification."""
from copy import deepcopy
from fractions import Fraction
from itertools import product, permutations
import json
from pathlib import Path
import sqlite3
import os
import tempfile
import unittest
from unittest.mock import patch
from ttauto_cert.exact_math import identity, matrix_multiply, charpoly
from ttauto_cert.spectral import reciprocal_polynomials, schur_inside, multiply, ordinary_orbits
from ttauto_cert.lefschetz import traces, compatible, compatible_traces
from ttauto_cert.cover_orbits import (Budget, joint_group, regular_from_anti,
                                     coupled_orbits, disc_patterns, disc_deck_feasible)
from ttauto_cert.cycle_obstruction import cycle_search
from ttauto_cert.perron_obstructions import strict_perron_obstruction
from ttauto_pipeline.config import normalize
from ttauto_pipeline.engine import run
from ttauto_pipeline.verify import verify
from ttauto_pipeline.automaton import Automaton
from ttauto_pipeline.pruning import Context, LefschetzModule
from ttauto_pipeline.proof_inventory import inventory
from ttauto_pipeline.proof_pruning import SpectrumContext, check, base_catalogue
from ttauto_pipeline.storage import atomic_json, canonical, digest

FIXTURES = Path(__file__).parent/'fixtures'
C = [1, 0, -1, -1, -1, 0, 1]


def disc(**changes):
    cfg = dict(schema='ttauto-proof-v1', surface=dict(kind='disc', punctures=3),
               upper='8/3', target_polynomial=[1, -3, 1])
    cfg.update(changes)
    return normalize(cfg)


def closed(kind='H', orders=(4,), **changes):
    cfg = dict(schema='ttauto-proof-v1', surface=dict(kind='closed', genus=3),
               strata=[dict(type=kind, orders=list(orders))], upper='701/500', target_polynomial=C)
    cfg.update(changes)
    return normalize(cfg)


class LiftRegressionTests(unittest.TestCase):
    def test_actual_genus_three_deck_action_is_not_minus_identity(self):
        data = json.loads((FIXTURES/'g3_q8_witness.json').read_text())
        a = Automaton(data, data['surface'])
        cfg = normalize(dict(schema='ttauto-experiment-v1', surface=data['surface'], upper='21/10'))
        context = Context(a, cfg, dict(start=0, path=[0]*5, terminal=0))
        cv = context.cover
        A, J, B = cv['positive_lift_matrix'], cv['deck_matrix'], cv['other_lift_matrix']
        self.assertEqual(cv['genus'], 5)
        self.assertEqual(sum(J[i][i] for i in range(10)), 2)
        self.assertEqual(matrix_multiply(J, J), identity(10))
        self.assertEqual(matrix_multiply(J, A), B)
        self.assertEqual(matrix_multiply(A, J), B)
        self.assertEqual(traces(charpoly(A), 4)[1:], [0, 4, 6, 20])
        self.assertEqual(traces(charpoly(B), 4)[1:], [-2, 4, -14, 20])
        d = LefschetzModule(cfg).check(context)
        self.assertEqual(d.outcome, 'pass')
        self.assertEqual(d.evidence['actual_lift_polynomials'], [charpoly(A), charpoly(B)])

    def test_foliation_sign_does_not_silently_negate_actual_traces(self):
        p = [1, -3, 1]
        t = traces(p, 12)[1:]
        self.assertEqual(compatible(p, [], True, stop=12),
                         compatible_traces([(-1)**k*v for k, v in enumerate(t, 1)], 2, [], True))
        self.assertFalse(compatible_traces(t, 2, [], True))
        self.assertTrue(compatible_traces(t, 2, [], False))


class HistoricalArithmeticTests(unittest.TestCase):
    def test_rank_six_catalogue_reproduces_original_75(self):
        rows = base_catalogue(6, '701/500', None)
        self.assertEqual(len(rows), 75)
        self.assertIn(tuple(C), rows)

    def test_newton_containing_catalogue_matches_direct_coefficients(self):
        bound = Fraction(3, 2)
        actual = {tuple(p) for p in reciprocal_polynomials(4, bound)}
        expected = set()
        for a, b in product(range(-8, 9), range(-24, 25)):
            p = [1, a, b, a, 1]
            t = traces(p, 2)
            if all(abs(t[k]) <= 2*(bound**k+bound**(-k)) for k in (1, 2)):
                expected.add(tuple(p))
        self.assertEqual(actual, expected)

    def test_schur_strict_boundary_and_hidden_rotating_factor(self):
        self.assertTrue(schur_inside([1, 0, 1], '3/2'))
        self.assertFalse(schur_inside([1, -2], 2))
        self.assertFalse(schur_inside([1, -3, 1], 2))
        p = multiply([1, 0, 0, -3, 0, 0, 1], [1, 1, 1])
        r = strict_perron_obstruction(p)
        self.assertTrue(r['excluded'])
        self.assertEqual(r['rotation_order'], 3)

    def test_joint_signatures_match_literal_prong_and_sheet_actions(self):
        horizon = 8
        for p, n in ((3, 2), (4, 2), (5, 1), (6, 1)):
            literal = set()
            for sigma in permutations(range(n)):
                for rotations in product(range(p), repeat=n):
                    lf, lt, z, e = ([0]*horizon for _ in range(4))
                    for start in range(n):
                        at, prong, parity = start, 0, 0
                        for k in range(horizon):
                            r = rotations[at]
                            prong = (prong+(p+1)*r) % (2*p) if p % 2 else (prong+r) % p
                            if not p % 2:
                                parity ^= r % 2
                            at = sigma[at]
                            if at != start:
                                continue
                            z[k] += 1
                            if p % 2:
                                lf[k] += 1-2*p if prong == 0 else 1
                                lt[k] += 1
                            elif parity == 0:
                                lf[k] += 2*(1-p if prong == 0 else 1); e[k] += 1
                            else:
                                lt[k] += 2; e[k] -= 1
                    literal.add(tuple((-(a-b)//2, c, d) for a, b, c, d in zip(lf, lt, z, e)))
            self.assertEqual({tuple(row) for row in joint_group(p, n, horizon)}, literal)

    def test_negative_orbits_use_odd_returns(self):
        horizon = 12
        O, N = {1: 3, 2: 5, 3: 1, 4: 2, 6: 3}, {1: 2, 2: 4, 3: 0, 4: 1, 6: 3}
        joint = list(joint_group(6, 1, horizon))[1]
        anti, base = [], []
        for k, (s, z, e) in enumerate(joint, 1):
            fixed = sum(d*v for d, v in O.items() if k % d == 0)
            negative = sum(d*v for d, v in N.items() if k % d == 0 and (k//d) % 2)
            anti.append(fixed+s); base.append(2-z+fixed+s+e-2*negative)
        regular, failure = regular_from_anti(anti, joint)
        self.assertIsNone(failure)
        model, failure = coupled_orbits(anti, base, joint, regular)
        self.assertIsNone(failure)
        self.assertEqual(model['negative_orbits'], [N.get(k, 0) for k in range(1, horizon+1)])
        base[1] += 1
        self.assertIsNotNone(coupled_orbits(anti, base, joint, regular)[1])

    def test_disc_branch_coupling_keeps_actual_torus_control(self):
        case = inventory(disc())['cases'][0]
        a, b = disc_patterns([1, -3, 1], case, 32)
        model = disc_deck_feasible([1, -3, 1], case, a, b, Budget(None))
        self.assertIsNotNone(model)
        self.assertEqual(sorted(model['regular_branch_cycles']), [1, 3])

    def test_original_D14_residual_cycle_obstruction(self):
        p = [1]
        factors = ([1, -1], [1, -1], [1, 0, 1], [1, 0, 1], [1, -1, 1],
                   [1, 1, -1, -2, 0, 1, 0, -2, -1, 1, 1])
        for f in factors:
            p = multiply(p, f)
        r = cycle_search(multiply(p, [1, 1]), 5000)
        self.assertTrue(r['completed'] and r['excluded'])
        self.assertEqual((r['nodes'], r['capacity_prunes'], r['matching_graphs']), (1936, 3008, 10))
        incomplete = cycle_search(multiply(p, [1, 1]), 1)
        self.assertFalse(incomplete['completed'] or incomplete['excluded'])
        control = cycle_search([1, -3, 1])
        self.assertFalse(control['excluded'])
        self.assertIsNotNone(control['survivor'])


class ProofEngineTests(unittest.TestCase):
    def test_complete_lower_without_automaton_or_external_process(self):
        with tempfile.TemporaryDirectory() as folder, patch('ttauto_pipeline.automaton.Automaton.construct',
                side_effect=AssertionError('automaton forbidden')), patch('subprocess.Popen', side_effect=AssertionError('process forbidden')):
            report = run(disc(), folder)
            self.assertTrue(report['lower_bound_proved'])
            self.assertFalse(report['minimum_proved'])
            self.assertFalse((Path(folder)/'automaton.json').exists())
            self.assertTrue(verify(folder, True)['verified'])

    def test_resume_exactly_replays_frontier_and_module_counts(self):
        cfg = closed()
        with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
            expected = run(cfg, a)
            limited = deepcopy(cfg); limited['limits']['events'] = 19
            first = run(limited, b)
            self.assertEqual(first['status'], 'paused_budget')
            self.assertFalse(first['lower_bound_proved'])
            actual = run(cfg, b, resume=True)
            self.assertEqual(actual, expected)
            self.assertTrue(verify(b, True)['lower_bound_proved'])
            changed = deepcopy(cfg); changed['pruning']['lefschetz'] = False
            with self.assertRaisesRegex(ValueError, 'resume refused'):
                run(changed, b, resume=True)

    def test_disabling_an_obstruction_retains_unresolved_candidates(self):
        cfg = closed(pruning={'lefschetz': False})
        with tempfile.TemporaryDirectory() as folder:
            report = run(cfg, folder)
            self.assertFalse(report['lower_bound_proved'])
            self.assertEqual(report['unresolved_candidates'], 1)

    def test_budget_is_not_a_mathematical_exclusion(self):
        with tempfile.TemporaryDirectory() as folder:
            report = run(closed(limits={'models': 0}), folder)
            self.assertFalse(report['lower_bound_proved'])
            self.assertEqual(report['counts']['stages']['stratum_lefschetz']['skipped'], 1)

    def test_catalogue_cannot_claim_completeness_and_is_portable(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            atomic_json(root/'catalog.json', {'polynomials': [[1, -3, 1]], 'completed': True})
            cfg = disc(source={'kind': 'polynomial_catalog', 'path': str(root/'catalog.json')})
            result = run(cfg, root/'run')
            self.assertTrue(result['enumeration_completed'])
            self.assertFalse(result['candidate_superset_complete'] or result['lower_bound_proved'])
            (root/'catalog.json').unlink()
            self.assertTrue(verify(root/'run', True)['verified'])

    def test_residual_continuation_keeps_parent_exclusions_and_coverage(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            parent = run(closed(pruning={'lefschetz': False}), root/'parent')
            self.assertEqual(parent['unresolved_candidates'], 1)
            cfg = closed(source={'kind': 'residual', 'path': str(root/'parent')})
            child = run(cfg, root/'child')
            self.assertEqual(child['counts']['candidates'], 1)
            self.assertTrue(child['lower_bound_proved'])
            self.assertTrue(verify(root/'child', True)['lower_bound_proved'])
            # The continuation carries its parent's evidence, so validation
            # does not depend on an external original run directory.
            (root/'parent'/'result.json').rename(root/'parent'/'hidden.json')
            self.assertTrue(verify(root/'child')['verified'])

    def test_partial_parent_cannot_be_upgraded_by_empty_residuals(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            run(closed(limits={'events': 1}), root/'parent')
            child = run(closed(source={'kind': 'residual', 'path': str(root/'parent')}), root/'child')
            self.assertTrue(child['enumeration_completed'])
            self.assertFalse(child['candidate_superset_complete'] or child['lower_bound_proved'])

    def test_mutated_claim_and_evidence_are_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            result = run(closed(), folder)
            result['minimum_proved'] = True
            atomic_json(Path(folder)/'result.json', result)
            with self.assertRaisesRegex(ValueError, 'claims'):
                verify(folder)
        with tempfile.TemporaryDirectory() as folder:
            run(closed(), folder)
            with sqlite3.connect(Path(folder)/'run.sqlite3') as db:
                number, payload = db.execute('SELECT id,payload FROM candidates LIMIT 1').fetchone()
                record = json.loads(payload); record['decisions'][0]['reason'] = 'forged'
                db.execute('UPDATE candidates SET payload=?,digest=? WHERE id=?', (canonical(record), digest(record), number))
            db.close()
            with self.assertRaisesRegex(ValueError, 'replay differs'):
                verify(folder)

    def test_lower_and_actual_map_are_combined_without_Hall(self):
        fold = normalize(dict(schema='ttauto-experiment-v1', surface=dict(kind='disc', inner_prongs=[1]*3,
                              outer_prongs=[1]), upper='8/3'))
        graph = Automaton.construct(fold).data
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); atomic_json(root/'witness.json', graph)
            cfg = disc(witness=dict(kind='fold_path', path=str(root/'witness.json'), start=0, folds=[0, 1], terminal=0))
            with patch('subprocess.Popen', side_effect=AssertionError('process forbidden')):
                result = run(cfg, root/'run')
                self.assertTrue(result['minimum_proved'])
                self.assertTrue(result['upper_bound_attained'])
                (root/'witness.json').unlink()
                self.assertTrue(verify(root/'run', True)['minimum_proved'])

    def test_inventory_covers_original_disc_reduction_and_closed_types(self):
        self.assertEqual(len(inventory(disc(surface=dict(kind='disc', punctures=10)))['cases']), 120)
        self.assertEqual(len(inventory(disc(surface=dict(kind='disc', punctures=12, one_pronged=True)))['cases']), 97)
        self.assertEqual(len(inventory(disc(surface=dict(kind='disc', punctures=14, one_pronged=True)))['cases']), 195)
        cases = inventory(closed(strata='all'))['cases']
        self.assertEqual(len(cases), 27)
        q8 = next(c for c in cases if c['type'] == 'Q' and c['orders'] == [8])
        self.assertEqual(q8['rank'], 4)

    def test_external_classifiers_and_unsupported_options_are_rejected(self):
        with self.assertRaises(ValueError):
            disc(witness={'kind': 'hall', 'path': 'unused'})
        with self.assertRaises(ValueError):
            disc(upper=1.5)
        with self.assertRaises(ValueError):
            inventory(closed(strata=[dict(type='Q', orders=[6])]))


class NativeProofSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[2]
        cls.executable = Path(os.environ.get('TTAUTO_BUILD', root/'build'))/'bin'/('ttauto_trace_prefix.exe' if os.name == 'nt' else 'ttauto_trace_prefix')
        if not cls.executable.is_file():
            raise unittest.SkipTest('native Newton source not built')

    def test_small_batches_resume_and_match_python_proof(self):
        cfg = closed(source={'kind': 'newton_native', 'executable': str(self.executable), 'batch_nodes': 37})
        with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
            expected = run(cfg, a)
            partial = deepcopy(cfg); partial['limits']['events'] = 9
            report = run(partial, b)
            self.assertFalse(report['lower_bound_proved'])
            resumed = run(cfg, b, resume=True)
            self.assertEqual(resumed, expected)
            self.assertTrue(verify(b, True)['lower_bound_proved'])
            self.assertEqual(resumed['counts']['native_statistics']['retained'], resumed['counts']['candidates'])
        with tempfile.TemporaryDirectory() as a:
            python = run(closed(), a)
            self.assertEqual(python['unresolved_candidates'], resumed['unresolved_candidates'])
            self.assertEqual(python['lower_bound_proved'], resumed['lower_bound_proved'])

    def test_native_candidate_set_matches_independent_python_recursion(self):
        from ttauto_pipeline.native_spectral import NativeSpectralSource
        from ttauto_cert.spectral import trace_bounds
        for enabled in (False, True):
            cfg = closed(source={'kind': 'newton_native', 'executable': str(self.executable), 'batch_nodes': 10000},
                         pruning={'trace': enabled})
            source = NativeSpectralSource(inventory(cfg), cfg)
            actual = []
            while (event := source.next_event()) is not None:
                if event['kind'] == 'candidate':
                    actual.append(tuple(event['polynomial']))
            expected = set(map(tuple, reciprocal_polynomials(6, cfg['upper'])))
            if enabled:
                bounds = trace_bounds(6, cfg['upper'], cfg['source']['trace_horizon'])
                expected = {p for p in expected if all(abs(t) <= b for t, b in zip(traces(p, len(bounds))[1:], bounds))}
            self.assertEqual(len(actual), len(set(actual)))
            self.assertEqual(set(actual), expected)

    def test_early_filters_are_optional_and_scope_is_preserved(self):
        for enabled in (False, True):
            cfg = disc(source={'kind': 'newton_native', 'executable': str(self.executable), 'batch_nodes': 5},
                       pruning={'trace': enabled, 'mod2': enabled})
            with tempfile.TemporaryDirectory() as folder:
                report = run(cfg, folder)
                self.assertTrue(report['lower_bound_proved'])
                self.assertTrue(verify(folder, True)['verified'])


if __name__ == '__main__':
    unittest.main()
