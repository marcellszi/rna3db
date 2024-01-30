from rna3db.cluster import InfernalGraph
import unittest


class TestInfernalGraph(unittest.TestCase):
    def setUp(self):
        self.graph = InfernalGraph()

    def test_add_family(self):
        self.graph.add_family("A")
        self.assertTrue("A" in self.graph.graph)
        self.assertEqual(self.graph.graph["A"]["is_family"], True)

    def test_add_chain(self):
        self.graph.add_chain("A")
        self.assertTrue("A" in self.graph.graph)
        self.assertEqual(self.graph.graph["A"]["is_family"], False)

    def test_add_edge(self):
        self.graph.add_family("A")
        self.graph.add_chain("B")
        self.graph.add_edge("A", "B")
        self.assertTrue("B" in self.graph.graph["A"]["neighbours"])
        self.assertTrue("A" in self.graph.graph["B"]["neighbours"])

    def test_dfs(self):
        self.graph.add_family("RF_A")
        self.graph.add_family("RF_B")

        self.graph.add_chain("A")
        self.graph.add_chain("B")
        self.graph.add_chain("C")
        self.graph.add_chain("D")
        self.graph.add_chain("E")

        self.graph.add_edge("RF_A", "A")
        self.graph.add_edge("RF_A", "B")

        self.graph.add_edge("RF_B", "C")
        self.graph.add_edge("RF_B", "D")

        components = self.graph.components()
        self.assertEqual(len(components), 3)
        self.assertIn(set(["A", "B"]), components)
        self.assertIn(set(["C", "D"]), components)
        self.assertIn(set(["E"]), components)


if __name__ == "__main__":
    unittest.main()
