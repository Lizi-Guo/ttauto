import copy
import json
from pathlib import Path
import unittest

from ttauto_cert.ribbon_verify import automorphisms_ordered, verify_graph

FIXTURES = Path(__file__).parent/'fixtures'


class RibbonVerification(unittest.TestCase):
    def load(self, name):
        return json.loads((FIXTURES/name).read_text())

    def test_actual_three_puncture_fold_maps(self):
        result = verify_graph(self.load('ribbon_n3_graph.json'))
        self.assertEqual(result['arcs'], 2)
        self.assertEqual(result['exported_terminal_groups_checked'], 1)

    def test_wrong_infinitesimal_orientation_is_rejected(self):
        # The main matrix and every endpoint remain unchanged. Only the
        # surface-boundary check sees this wrong winding around a puncture.
        data = self.load('ribbon_n3_graph.json')
        data['vertices'][0]['arcs'][0]['edge_images'][0][1] = 4
        with self.assertRaisesRegex(ValueError, 'oriented surface boundaries'):
            verify_graph(data)

    def test_matrix_and_transport_corruption_are_rejected(self):
        data = self.load('ribbon_n3_graph.json')
        bad = copy.deepcopy(data)
        bad['vertices'][0]['arcs'][0]['matrix'][0][0] = 9
        with self.assertRaisesRegex(ValueError, 'matrix mismatch'):
            verify_graph(bad)
        data['vertices'][0]['arcs'][0]['vertex_images'][0] = 99
        with self.assertRaises(ValueError):
            verify_graph(data)

    def test_terminal_group_must_be_complete_and_identity_first(self):
        graph = self.load('ribbon_n4_symmetric_state.json')
        group = automorphisms_ordered(graph)
        self.assertEqual(len(group), 2)
        self.assertEqual(group[0]['edge_images'], [[e+1] for e in range(7)])
        self.assertEqual(group[1]['edge_images'], [[-3],[-2],[-1],[7],[6],[5],[4]])
        data = {'n': 4, 'stratum_index': 0, 'vertices': [
            {'id': 0, 'ribbon_graph': graph, 'arcs': [], 'automorphisms': group}]}
        self.assertEqual(verify_graph(data)['exported_terminal_groups_checked'], 1)
        data['vertices'][0]['automorphisms'] = group[:1]
        with self.assertRaisesRegex(ValueError, 'incomplete or misordered'):
            verify_graph(data)
        data['vertices'][0]['automorphisms'] = group[::-1]
        with self.assertRaisesRegex(ValueError, 'incomplete or misordered'):
            verify_graph(data)


if __name__ == '__main__':
    unittest.main()
