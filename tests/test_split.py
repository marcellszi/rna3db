import unittest

from rna3db.split import find_optimal_components, split
from rna3db.utils import write_json

import unittest
import tempfile


def total_deviation(sol, components, bins):
    weights = [sum([components[i] for i in s]) for s in sol]
    deviations = [abs(w - c) for w, c in zip(weights, bins)]
    return sum(deviations)


class TestFindOptimalComponents(unittest.TestCase):
    def test_solver_small(self):
        components, bins = [3, 4, 1, 2], [1, 5, 4]
        sol = find_optimal_components(components, bins)
        self.assertEqual(total_deviation(sol, components, bins), 0)

        components, bins = [2, 2, 1, 2, 3], [1, 6, 3]
        sol = find_optimal_components(components, bins)
        self.assertEqual(total_deviation(sol, components, bins), 0)

        components, bins = [1, 4, 2, 3], [5, 1, 4]
        sol = find_optimal_components(components, bins)
        self.assertEqual(total_deviation(sol, components, bins), 0)

        components, bins = [3, 4, 3], [4, 6]
        sol = find_optimal_components(components, bins)
        self.assertEqual(total_deviation(sol, components, bins), 0)

        components, bins = [2, 1, 4, 3], [4, 5, 1]
        sol = find_optimal_components(components, bins)
        self.assertEqual(total_deviation(sol, components, bins), 0)

    def test_solver_large(self):
        # optimal sol of all of these has been manually verified
        # there are multiple solutions in many cases (all?)
        # but the minimum total deviation is guaranteed correct

        components, bins = [85, 6, 781, 1, 14, 63, 50], [128, 872]
        sol = find_optimal_components(components, bins)
        self.assertEqual(total_deviation(sol, components, bins), 0)

        components, bins = [14, 44, 404, 70, 417, 51], [893, 107]
        sol = find_optimal_components(components, bins)
        self.assertEqual(total_deviation(sol, components, bins), 4)

        components, bins = [8, 230, 150, 5, 124, 483], [13, 694, 293]
        sol = find_optimal_components(components, bins)
        self.assertEqual(total_deviation(sol, components, bins), 38)

        components, bins = [1, 164, 274, 41, 176, 1, 168, 175], [41, 716, 243]
        sol = find_optimal_components(components, bins)
        self.assertEqual(total_deviation(sol, components, bins), 62)

        components, bins = [4, 228, 2, 1, 477, 282, 6], [5, 594, 401]
        sol = find_optimal_components(components, bins)
        self.assertEqual(total_deviation(sol, components, bins), 152)


class TestSplit(unittest.TestCase):
    def test_real_case(self):
        example = {
            "component_1": {
                "8g9s_O": {
                    "8g9s_O": {
                        "release_date": "2024-03-06",
                        "structure_method": "electron microscopy",
                        "resolution": 3.4,
                        "length": 42,
                        "sequence": "AUUGAAACAGGGUCAGCUUGCCGUAGGUGGCAUCGCCCUCGU",
                    },
                },
                "8g9t_O": {
                    "8g9t_O": {
                        "release_date": "2024-03-06",
                        "structure_method": "electron microscopy",
                        "resolution": 3.6,
                        "length": 43,
                        "sequence": "GAAACAGGGUCAGCUUGCCGUAGGUGGCAUCGCCCUCGUAAAA",
                    }
                },
            },
            "component_2": {
                "6gov_R": {
                    "6gov_R": {
                        "release_date": "2019-02-13",
                        "structure_method": "electron microscopy",
                        "resolution": 3.7,
                        "length": 66,
                        "sequence": "GGCGCUCUUUAACAUUAAGCCCUGAAGAAGGGCAAAAAUCAAAUUAAACCACACCUGGCGUGUGGC",
                    },
                },
            },
            "component_3": {
                "5une_B": {
                    "5une_B": {
                        "release_date": "2017-03-29",
                        "structure_method": "x-ray diffraction",
                        "resolution": 2.90009915777,
                        "length": 47,
                        "sequence": "GGGAAAUGAUGGGCGUAGACGCACGUCAGCGGCGGAAAUGGUUUCCC",
                    },
                    "5une_A": {
                        "release_date": "2017-03-29",
                        "structure_method": "x-ray diffraction",
                        "resolution": 2.90009915777,
                        "length": 47,
                        "sequence": "GGGAAAUGAUGGGCGUAGACGCACGUCAGCGGCGGAAAUGGUUUCCC",
                    },
                },
            },
            "component_4": {
                "7wii_V": {
                    "7wii_V": {
                        "release_date": "2023-01-18",
                        "structure_method": "x-ray diffraction",
                        "resolution": 2.75,
                        "length": 50,
                        "sequence": "GGCGUGGUCCGUUCAACUCGUUCCUCGAAAGAGGAACUACGGGAGACGCC",
                    },
                    "7wi9_V": {
                        "release_date": "2023-01-18",
                        "structure_method": "x-ray diffraction",
                        "resolution": 2.98,
                        "length": 50,
                        "sequence": "GGCGUGGUCCGUUCAACUCGUUCCUCGAAAGAGGAACUACGGGAGACGCC",
                    },
                    "7wie_V": {
                        "release_date": "2023-01-18",
                        "structure_method": "x-ray diffraction",
                        "resolution": 2.9,
                        "length": 50,
                        "sequence": "GGCGUGGUCCGUUCAACUCGUUCCUCGAAAGAGGAACUACGGGAGACGCC",
                    },
                    "7wib_V": {
                        "release_date": "2023-01-18",
                        "structure_method": "x-ray diffraction",
                        "resolution": 2.83,
                        "length": 50,
                        "sequence": "GGCGUGGUCCGUUCAACUCGUUCCUCGAAAGAGGAACUACGGGAGACGCC",
                    },
                    "7wif_V": {
                        "release_date": "2023-01-18",
                        "structure_method": "x-ray diffraction",
                        "resolution": 2.86,
                        "length": 50,
                        "sequence": "GGCGUGGUCCGUUCAACUCGUUCCUCGAAAGAGGAACUACGGGAGACGCC",
                    },
                },
                "7wia_V": {
                    "7wia_V": {
                        "release_date": "2023-01-18",
                        "structure_method": "x-ray diffraction",
                        "resolution": 3.22,
                        "length": 50,
                        "sequence": "GGCGUGGUCCGUUCAAGUCGUUCCUCGAAAGAGGAACUACGGGAGACGCC",
                    }
                },
            },
        }

        splits_n = [3, 1, 2]
        splits = [i / sum(splits_n) for i in splits_n]

        with tempfile.NamedTemporaryFile() as f:
            write_json(example, f.name)
            f.seek(0)
            ans = split(f.name, splits=splits)
            for i, (K, V) in enumerate(ans.items()):
                curr_length = sum(len(v) for v in V.values())
                # check if we have desired length
                self.assertEqual(curr_length, splits_n[i])
                for k, v in V.items():
                    # check that we have the data we expect
                    # (i.e. release_date, structure_method, etc.)
                    self.assertEqual(v, example[k])

    def test_force_zero_last(self):
        example = {
            "component_0": "0",
            "component_1": "01",
            "component_2": "01",
        }
        splits_n = [1, 0, 4]
        splits = [i / sum(splits_n) for i in splits_n]

        with tempfile.NamedTemporaryFile() as f:
            write_json(example, f.name)
            f.seek(0)

            ans = split(f.name, splits=splits)
            res = [sum([len(v) for v in V.values()]) for V in ans.values()]
            self.assertEqual(sum([abs(w - c) for w, c in zip(splits_n, res)]), 0)

            ans = split(f.name, splits=splits, force_zero_last=True)
            res = [sum([len(v) for v in V.values()]) for V in ans.values()]
            self.assertEqual(sum([abs(w - c) for w, c in zip(splits_n, res)]), 2)


if __name__ == "__main__":
    unittest.main()
