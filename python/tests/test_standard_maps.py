import copy
from fractions import Fraction
import json
from pathlib import Path
from types import SimpleNamespace
import unittest

from ttauto_cert.standard_maps import (StandardMapSearch, StandardState,
                                      exact_sum_bound, primitive,
                                      surjective_free_group, lazy_product)

FIXTURES = Path(__file__).resolve().parent/'fixtures'


def graph():
    return json.loads((FIXTURES/'ribbon_n3_graph.json').read_text())['vertices'][0]['ribbon_graph']


class StandardMapTests(unittest.TestCase):
    def test_long_word_no_longer_depends_on_python_recursion_limit(self):
        search = StandardMapSearch(graph(), 50, max_prefixes=2000)
        word = next(search._words(0, 0, 1500))
        self.assertEqual(len(word), 2999)
        self.assertEqual(sum(abs(d) <= search.state.m for d in word), 1500)
        self.assertEqual(search.state.origins[word[0]], 0)
        self.assertEqual(search.state.origins[-word[-1]], 0)
        self.assertTrue(all(search.state.origins[-a] == search.state.origins[b]
                            for a, b in zip(word, word[1:])))
        self.assertEqual(search.counts['word_prefixes'], 1500)

    def test_assignment_stack_handles_more_edges_than_recursion_limit(self):
        search = StandardMapSearch.__new__(StandardMapSearch)
        size = 1500
        search.state = SimpleNamespace(m=size, graph={'endpoints': [(i, i+1) for i in range(size)]})
        search._words = lambda start, finish, length: iter([[start+1]] if length == 1 else [])
        result = next(search._assign(list(range(size+1)), [None]*size, size))
        self.assertEqual(result, [[i+1] for i in range(size)])

    def test_cartesian_factors_are_lazy_and_stack_safe(self):
        observed = []
        def large_factor():
            for i in range(1000000):
                observed.append(i)
                yield i
        values = lazy_product([large_factor, lambda: iter((3, 4))])
        self.assertEqual(next(values), (0, 3))
        self.assertEqual(observed, [0])
        values.close()
        self.assertEqual(next(lazy_product([lambda: iter((1,))]*1500)), (1,)*1500)
        self.assertEqual(list(lazy_product([])), [()])
        self.assertEqual(list(lazy_product([lambda: iter(())])), [])

    def test_full_n3_direct_enumeration(self):
        search = StandardMapSearch(graph(), Fraction(2619, 1000))
        maps = list(search.run())
        self.assertEqual(search.status, 'completed')
        self.assertEqual(search.counts['endpoint_compatible_maps'], 2090)
        self.assertEqual(search.counts['primitive_maps'], 1340)
        self.assertEqual(search.counts['primitive_below_upper'], 76)
        self.assertEqual(search.counts['free_group_automorphisms'], 28)
        self.assertEqual(len(maps), 2)
        self.assertTrue(all(m['charpoly'] == [1, -3, 1] for m in maps))
        self.assertTrue(all(not m['pseudo_anosov_certified'] for m in maps))
        self.assertFalse(search.summary()['all_stratum_states_certified'])

    def test_prefix_limit_and_early_close_are_not_completion(self):
        search = StandardMapSearch(graph(), 3, max_prefixes=2)
        list(search.run())
        self.assertEqual(search.status, 'truncated')
        self.assertEqual(search.counts['word_prefixes'], 2)
        self.assertFalse(search.summary()['completed'])
        stopped = StandardMapSearch(graph(), Fraction(2619, 1000))
        generator = stopped.run()
        next(generator)
        generator.close()
        self.assertEqual(stopped.status, 'aborted')

    def test_zero_map_budget_and_exact_bound(self):
        search = StandardMapSearch(graph(), 3, max_maps=0)
        self.assertEqual(list(search.run()), [])
        self.assertEqual(search.status, 'truncated')
        self.assertEqual(exact_sum_bound(Fraction(2619, 1000), 2), 7)
        with self.assertRaises(TypeError):
            exact_sum_bound(2.619, 2)
        with self.assertRaises(ValueError):
            exact_sum_bound(1, 2)

    def test_stallings_detects_nonsurjective_image(self):
        self.assertTrue(surjective_free_group([[1, 2], [2]], 2))
        self.assertTrue(surjective_free_group([[2, 1, -2], [2]], 2))
        self.assertFalse(surjective_free_group([[1, 1], [2]], 2))
        self.assertFalse(surjective_free_group([[1, 2, -1], [2]], 2))
        self.assertFalse(surjective_free_group([[], [2]], 2))
        # Identity on abelianization, and an image subgroup of rank two,
        # but not the whole F_2. Neither rank nor determinant detects this.
        self.assertFalse(surjective_free_group([[1], [2, 1, 2, -1, -2]], 2))

    def test_primitive_period_and_nonstandard_switch(self):
        self.assertTrue(primitive([[2, 1], [1, 1]]))
        self.assertFalse(primitive([[0, 1], [1, 0]]))
        self.assertFalse(primitive([[2, 0], [0, 3]]))
        bad = copy.deepcopy(graph())
        bad['cyclic_directions'][1] = [4, -1, -4, 2]
        with self.assertRaises(ValueError):
            StandardState(bad)

    def test_colored_polygon_transport_is_not_terminal_symmetry(self):
        state = StandardState(graph())
        transports = list(state.peripheral_transports())
        self.assertEqual(len(transports), 6)
        self.assertEqual({tuple(v) for v, _, _ in transports},
                         {(0,1,2), (0,2,1), (1,0,2), (1,2,0), (2,0,1), (2,1,0)})
        marked = copy.deepcopy(graph())
        marked['labels'] = [1, 2, 3]
        self.assertEqual(len(list(StandardState(marked).peripheral_transports())), 1)

    def test_trigon_rotations_transport_all_endpoints_and_free_basis(self):
        graph = json.loads((FIXTURES/'ribbon_n4_trigon_state.json').read_text())
        state = StandardState(graph)
        transports = list(state.peripheral_transports())
        self.assertEqual(len(transports), 24*3)
        rotations = set()
        for vertices, words, polygon_images in transports:
            self.assertEqual(polygon_images[1], 1)  # The zero is a separate color.
            rotations.add(tuple(words[e-1][0] for e in (6, 7, 8)))
            for e in range(state.m+1, state.e+1):
                image = words[e-1][0]
                self.assertEqual(state.origins[image], vertices[state.origins[e]])
                self.assertEqual(state.origins[-image], vertices[state.origins[-e]])
        self.assertEqual(rotations, {(6,7,8), (7,8,6), (8,6,7)})
        identity = [[e] for e in range(1, state.e+1)]
        images, rejection = state.surface_automorphism(identity, list(range(state.q)))
        self.assertIsNone(rejection)
        self.assertEqual(images, [[i] for i in range(1, state.q+1)])
        search = StandardMapSearch(graph, Fraction(23,10), max_prefixes=100)
        list(search.run())
        self.assertEqual(search.status, 'truncated')


if __name__ == '__main__':
    unittest.main()
