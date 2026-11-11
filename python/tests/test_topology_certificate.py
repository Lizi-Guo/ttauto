"""Geometric reconstruction and independent affine positive/negative examples."""
import copy
import json
from pathlib import Path
import unittest

from ttauto_cert.topology_certificate import (
    TopologyCertificateError, UnsupportedTopology,
    reconstruct_prong_certificate, validate_reported_prongs,
    certify_pillowcase_affine,
)

FIXTURES = Path(__file__).parent / 'fixtures'


def load(name):
    return json.loads((FIXTURES / name).read_text(encoding='utf-8'))


def _subdivide_peripheral_loops(graph):
    """Insert bivalent vertices, with an explicit conjugate edge map."""
    original = copy.deepcopy(graph)
    result = copy.deepcopy(graph)
    by_label = {e['label']: e for e in result['edges']}
    old_vertex = {v['label']: v for v in result['vertices']}
    peripheral = [e for e in original['edges'] if e['type'] == 1]
    first_edge = max(by_label) + 1
    first_vertex = max(old_vertex) + 1
    halves = {e['label']: first_edge + i for i, e in enumerate(peripheral)}
    middles = {e['label']: first_vertex + i for i, e in enumerate(peripheral)}

    def expand(word):
        out = []
        for d in word:
            if abs(d) not in halves:
                out.append(d)
            elif d > 0:
                out.extend([d, halves[d]])
            else:
                out.extend([-halves[-d], d])
        return out

    for v in result['vertices']:
        v['cyclic_directions'] = [-halves[-d] if d < 0 and -d in halves else d
                                  for d in v['cyclic_directions']]
    for edge in original['edges']:
        e = by_label[edge['label']]
        if e['label'] not in halves:
            e['image'] = expand(e['image'])
            continue
        if not (edge['start'] == edge['end'] and len(edge['image']) == 1 and edge['image'][0] > 0):
            raise AssertionError('fixture is not a positive peripheral loop permutation')
        label, target = e['label'], edge['image'][0]
        e['end'] = middles[label]
        e['image'] = [target]
        result['edges'].append({'label': halves[label], 'start': middles[label],
                                'end': edge['end'], 'type': 1, 'puncture': e['puncture'],
                                'image': [halves[target]]})
        result['vertices'].append({'label': middles[label], 'image': middles[target],
                                   'cyclic_directions': [-label, halves[label]]})
    return result


class GateBlowupTests(unittest.TestCase):
    def test_reconstructs_regions_without_reading_backend_type_or_prongs(self):
        data = load('bh_n6_s3_topology.json')
        data.pop('type')
        data.pop('singularity_orbits')
        result = reconstruct_prong_certificate(data['final_graph'], 6)
        self.assertEqual(result['reconstructed_interior_prongs'], [5])
        self.assertEqual(result['reconstructed_outer_prongs'], 1)
        self.assertEqual(result['reconstructed_puncture_prongs'], [1] * 6)
        self.assertFalse(result['pseudo_anosov_certified'])
        self.assertFalse(result['invariant_foliation_realization_certified'])
        interior = next(r for r in result['regions'] if r['kind'] == 'interior')
        self.assertEqual(len(interior['boundary']), 5)
        self.assertEqual(len(interior['cusps']), 5)
        self.assertEqual(interior['projected_boundary'], [])

    def test_real_reports_and_marked_regular_outer(self):
        examples = [('bh_n3_s0_0.json', [], 1, [1] * 3),
                    ('bh_n5_s1_0.json', [3], 2, [1] * 5),
                    ('bh_n6_s3_topology.json', [5], 1, [1] * 6),
                    ('bh_fully_punctured_n6_s2.json', [], 2, [1, 4, 1, 1, 1, 1, 1])]
        for name, interior, outer, punctures in examples:
            with self.subTest(name=name):
                data = load(name)
                before = copy.deepcopy(data)
                result = validate_reported_prongs(data)
                self.assertEqual(data, before)
                self.assertEqual(result['reconstructed_interior_prongs'], interior)
                self.assertEqual(result['reconstructed_outer_prongs'], outer)
                self.assertEqual(result['reconstructed_puncture_prongs'], punctures)

    def test_euler_preserving_forged_prongs_rejected_by_region_checker(self):
        data = load('bh_n6_s3_topology.json')
        for orbit in data['singularity_orbits']:
            for s in orbit['singularities']:
                if s['interior']:
                    s['prongs'] = 4
                elif len(s['location']) > 1:
                    s['prongs'] = 2
        # This is a malformed REPORT, not a mathematical exclusion or witness.
        with self.assertRaisesRegex(TopologyCertificateError, 'reconstructed cusp count'):
            validate_reported_prongs(data)

    def test_bad_orbit_and_reversed_boundary_are_rejected(self):
        data = load('bh_n6_s3_topology.json')
        orbit = next(o['singularities'] for o in data['singularity_orbits'] if len(o['singularities']) > 1)
        orbit[0], orbit[1] = orbit[1], orbit[0]
        with self.assertRaisesRegex(TopologyCertificateError, 'orbit is incorrect'):
            validate_reported_prongs(data)
        data = load('bh_n6_s3_topology.json')
        singularity = data['singularity_orbits'][1]['singularities'][0]
        singularity['location'] = [-d for d in reversed(singularity['location'])]
        with self.assertRaisesRegex(TopologyCertificateError, 'does not match'):
            validate_reported_prongs(data)

    def test_duplicate_omitted_and_wrong_location_rejected(self):
        data = load('bh_n6_s3_topology.json')
        data['singularity_orbits'].pop(0)
        with self.assertRaisesRegex(TopologyCertificateError, 'omits'):
            validate_reported_prongs(data)
        data = load('bh_n6_s3_topology.json')
        data['singularity_orbits'].append(copy.deepcopy(data['singularity_orbits'][0]))
        with self.assertRaisesRegex(TopologyCertificateError, 'reported twice'):
            validate_reported_prongs(data)
        data = load('bh_n6_s3_topology.json')
        data['singularity_orbits'][0]['singularities'][0]['location'] = [1]
        with self.assertRaisesRegex(TopologyCertificateError, 'does not match'):
            validate_reported_prongs(data)

    def test_two_gate_subdivision_does_not_manufacture_singularities(self):
        data = load('bh_n6_s3_topology.json')
        graph = _subdivide_peripheral_loops(data['final_graph'])
        result = reconstruct_prong_certificate(graph, 6)
        self.assertEqual(result['reconstructed_interior_prongs'], [5])
        self.assertEqual(result['reconstructed_outer_prongs'], 1)
        self.assertEqual(result['reconstructed_puncture_prongs'], [1] * 6)
        self.assertEqual(sum(len(g['gates']) == 2 for g in result['gate_data']), 6)

    def test_cyclic_start_changes_do_not_change_regions(self):
        data = load('bh_n6_s3_topology.json')
        for v in data['final_graph']['vertices']:
            ds = v['cyclic_directions']
            v['cyclic_directions'] = ds[1:] + ds[:1]
        # In particular this splits the old [11,12] gate across a list endpoint.
        result = validate_reported_prongs(data)
        self.assertEqual(result['reconstructed_interior_prongs'], [5])

    def test_known_reducible_example_is_outside_model(self):
        data = load('bh_n4_s0_0.json')
        data['type'] = 'pseudo_anosov'
        with self.assertRaises(UnsupportedTopology):
            reconstruct_prong_certificate(data['final_graph'], 4)

    def test_corrupt_graph_and_marking_rejected(self):
        for kind in ('unknown_edge', 'duplicate_label', 'wrong_marking', 'backtrack', 'bad_rotation'):
            with self.subTest(kind=kind):
                graph = load('bh_n6_s3_topology.json')['final_graph']
                if kind == 'unknown_edge':
                    graph['edges'][0]['image'] = [999]
                elif kind == 'duplicate_label':
                    graph['edges'].append(copy.deepcopy(graph['edges'][0]))
                elif kind == 'wrong_marking':
                    graph['edges'][0]['puncture'] = 2
                elif kind == 'backtrack':
                    graph['edges'][0]['image'] = [6, -6, 6]
                else:
                    graph['vertices'][0]['cyclic_directions'].append(1)
                with self.assertRaises(TopologyCertificateError):
                    reconstruct_prong_certificate(graph, 6)


ORIGIN = [[0, 1], [0, 1]]
BRANCHES = [[[i, 2], [j, 2]] for i in (0, 1) for j in (0, 1)]
MATRIX = [[2, 1], [1, 1]]


class AffinePillowcaseTests(unittest.TestCase):
    def test_three_puncture_attainer_is_a_complete_independent_pa_certificate(self):
        result = certify_pillowcase_affine(MATRIX, BRANCHES[1:], ORIGIN)
        self.assertTrue(result['pseudo_anosov_certified'])
        self.assertTrue(result['invariant_foliation_realization_certified'])
        self.assertFalse(result['uses_backend_type_or_prong_claims'])
        self.assertEqual(result['dilatation_minimal_polynomial'], [1, -3, 1])
        self.assertEqual(result['puncture_prongs'], [1, 1, 1])
        self.assertEqual(result['outer_prongs'], 1)
        self.assertEqual(sorted(map(len, result['puncture_orbits'])), [3])

    def test_four_puncture_attainer_has_a_fixed_regular_outer_mark(self):
        result = certify_pillowcase_affine(MATRIX, BRANCHES, [[2, 5], [-1, 5]])
        self.assertEqual(result['puncture_prongs'], [1] * 4)
        self.assertEqual(result['outer_prongs'], 2)
        self.assertEqual(result['dilatation_minimal_polynomial'], [1, -3, 1])

    def test_negative_trace_and_integral_coordinate_translates(self):
        positive = certify_pillowcase_affine(MATRIX, BRANCHES[1:], ORIGIN)
        negative = certify_pillowcase_affine([[-2, -1], [-1, -1]], BRANCHES[1:], [[1, 1], [-3, 1]])
        self.assertEqual(positive['dilatation_minimal_polynomial'], negative['dilatation_minimal_polynomial'])
        self.assertEqual(positive['puncture_permutation'], negative['puncture_permutation'])

    def test_nonhyperbolic_wrong_determinant_and_noninteger_matrix_rejected(self):
        for matrix in ([[1, 1], [0, 1]], [[0, -1], [1, 0]], [[3, 1], [1, 1]], [[2.0, 1], [1, 1]]):
            with self.subTest(matrix=matrix):
                with self.assertRaises(TopologyCertificateError):
                    certify_pillowcase_affine(matrix, BRANCHES[1:], ORIGIN)

    def test_missing_duplicate_noninvariant_and_unfixed_marks_rejected(self):
        bad_marks = [(BRANCHES[1:3], ORIGIN),
                     (BRANCHES + [ORIGIN], [[2, 5], [-1, 5]]),
                     (BRANCHES + [[[1, 7], [0, 1]]], [[2, 5], [-1, 5]]),
                     (BRANCHES, [[1, 7], [0, 1]])]
        for interior, outer in bad_marks:
            with self.subTest(interior=interior, outer=outer):
                with self.assertRaises(TopologyCertificateError):
                    certify_pillowcase_affine(MATRIX, interior, outer)

    def test_float_and_zero_denominator_marks_rejected(self):
        for outer in ([[0.0, 1], [0, 1]], [[0, 0], [0, 1]]):
            with self.assertRaises(TopologyCertificateError):
                certify_pillowcase_affine(MATRIX, BRANCHES[1:], outer)


if __name__ == '__main__':
    unittest.main()
