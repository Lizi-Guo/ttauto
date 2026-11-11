import copy
from fractions import Fraction
from itertools import product
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

from ttauto_cert.cycle_pruning import CycleGrowthPruner, verify_cycle_witness
from ttauto_cert.exact_math import compare_perron_to_rational
from ttauto_cert.standard_maps import StandardMapSearch

FIXTURES = Path(__file__).with_name('fixtures')


class CyclePruningTests(unittest.TestCase):
    def test_all_small_matrices_against_independent_sturm_comparison(self):
        for upper in (Fraction(3, 2), Fraction(2), Fraction(5, 2)):
            pruner = CycleGrowthPruner(upper)
            rejected = 0
            for entries in product(range(4), repeat=4):
                matrix = [entries[:2], entries[2:]]
                witness = pruner.witness(matrix)
                if witness:
                    rejected += 1
                    self.assertTrue(verify_cycle_witness(witness))
                    self.assertGreater(compare_perron_to_rational(matrix, upper), 0)
            self.assertGreater(rejected, 0)

    def test_three_vertex_cycles_and_entrywise_completions(self):
        rng = random.Random(20261028)
        pruner = CycleGrowthPruner(Fraction(3, 2))
        for _ in range(80):
            partial = [[rng.randrange(3) for _ in range(3)] for _ in range(3)]
            witness = pruner.witness(partial)
            if witness:
                self.assertTrue(verify_cycle_witness(witness))
                complete = [[v+rng.randrange(2) for v in row] for row in partial]
                self.assertGreater(compare_perron_to_rational(complete, Fraction(3, 2)), 0)
            columns = list(zip(*partial))
            support = pruner.excluded_leading_columns(columns, 3)
            if support is not None:
                block = [[partial[i][j] for j in support] for i in support]
                self.assertGreater(compare_perron_to_rational(block, Fraction(3, 2)), 0)
                self.assertTrue(verify_cycle_witness(pruner.witness(block)))

    def test_equality_imprimitive_and_truncation_are_safe(self):
        for matrix in ([[2]], [[0, 2], [2, 0]], [[2, 0], [0, 2]]):
            self.assertIsNone(CycleGrowthPruner(2).witness(matrix))
        # A prefix with no cycles may still have an expanding completion.
        self.assertIsNone(CycleGrowthPruner(2).witness([[0, 0], [3, 0]]))
        # Bounded discovery can miss an obstruction but cannot create one.
        pruner = CycleGrowthPruner(2, max_length=1, max_cycles=1, max_steps=1)
        witness = pruner.witness([[0, 3], [3, 0]])
        self.assertIsNone(witness)

    def test_clique_bound_is_stronger_than_a_single_anchor(self):
        matrix = [[1, 1, 0, 0, 1], [0, 0, 1, 1, 0], [1, 0, 1, 0, 0],
                  [0, 1, 0, 0, 1], [1, 1, 0, 0, 1]]
        witness = CycleGrowthPruner(2).witness(matrix)
        self.assertEqual(witness['rule'], 'induced_cycle_clique')
        self.assertTrue(verify_cycle_witness(witness))
        self.assertGreater(compare_perron_to_rational(matrix, 2), 0)

    def test_witness_tampering_and_bad_inputs(self):
        witness = CycleGrowthPruner(Fraction(3, 2)).witness([[2]])
        for mutate in (lambda w: w['cycles'][0].update(multiplicity=3),
                       lambda w: w['cycles'].append(copy.deepcopy(w['cycles'][0])),
                       lambda w: w.update(negative_integer_value=0)):
            bad = copy.deepcopy(witness)
            mutate(bad)
            with self.assertRaises(ValueError):
                verify_cycle_witness(bad)
        with self.assertRaises(TypeError):
            CycleGrowthPruner(1.5)
        with self.assertRaises(ValueError):
            CycleGrowthPruner(2).witness([[1, -1], [0, 1]])

    def test_full_standard_search_preserves_complete_signed_maps(self):
        graph = json.loads((FIXTURES/'ribbon_n3_graph.json').read_text())['vertices'][0]['ribbon_graph']
        baseline = StandardMapSearch(graph, Fraction(2619, 1000))
        enhanced = StandardMapSearch(graph, Fraction(2619, 1000), cycle_pruning=True)
        self.assertEqual(list(baseline.run()), list(enhanced.run()))
        self.assertEqual(enhanced.status, 'completed')
        self.assertGreater(enhanced.counts['cycle_growth_prunes'], 0)
        self.assertLess(enhanced.counts['endpoint_compatible_maps'], baseline.counts['endpoint_compatible_maps'])
        for witness in enhanced.cycle_witnesses:
            self.assertTrue(verify_cycle_witness(witness))



if __name__ == '__main__':
    unittest.main()
