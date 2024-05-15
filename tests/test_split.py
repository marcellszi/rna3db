import unittest

from rna3db.split import find_optimal_components

import unittest


class TestFindOptimalComponents(unittest.TestCase):
    def test_small_lengths(self):
        lengths_dict = {"component1": 1, "component2": 2, "component3": 3}
        capacity = 5
        expected_components = set(["component2", "component3"])
        expected_length = sum(
            lengths_dict[component] for component in expected_components
        )
        result_components = find_optimal_components(lengths_dict, capacity)
        result_length = sum(lengths_dict[component] for component in result_components)
        self.assertEqual(result_length, expected_length)
        self.assertCountEqual(result_components, expected_components)

    def test_exact_fit(self):
        lengths_dict = {"component1": 4, "component2": 3, "component3": 2}
        capacity = 7
        expected_components = set(["component1", "component2"])
        expected_length = sum(
            lengths_dict[component] for component in expected_components
        )
        result_components = find_optimal_components(lengths_dict, capacity)
        result_length = sum(lengths_dict[component] for component in result_components)
        self.assertEqual(result_length, expected_length)
        self.assertCountEqual(result_components, expected_components)

    def test_over_capacity(self):
        lengths_dict = {"component1": 8, "component2": 9, "component3": 10}
        capacity = 5
        expected_components = set([])
        expected_length = sum(
            lengths_dict[component] for component in expected_components
        )
        result_components = find_optimal_components(lengths_dict, capacity)
        result_length = sum(lengths_dict[component] for component in result_components)
        self.assertEqual(result_length, expected_length)
        self.assertEqual(result_components, expected_components)

    def test_zero_capacity(self):
        lengths_dict = {"component1": 1, "component2": 2, "component3": 3}
        capacity = 0
        expected_components = set([])
        expected_length = sum(
            lengths_dict[component] for component in expected_components
        )
        result_components = find_optimal_components(lengths_dict, capacity)
        result_length = sum(lengths_dict[component] for component in result_components)
        self.assertEqual(result_length, expected_length)
        self.assertEqual(result_components, expected_components)

    def test_large_numbers(self):
        lengths_dict = {"component1": 100, "component2": 300, "component3": 400}
        capacity = 800
        expected_components = set(["component1", "component2", "component3"])
        expected_length = sum(
            lengths_dict[component] for component in expected_components
        )
        result_components = find_optimal_components(lengths_dict, capacity)
        result_length = sum(lengths_dict[component] for component in result_components)
        self.assertEqual(result_length, expected_length)
        self.assertCountEqual(result_components, expected_components)


if __name__ == "__main__":
    unittest.main()
