import os
import unittest

import pandas as pd
from ndex2.cx2 import RawCX2NetworkFactory

from cellmaps_hierarchyeval.perturb import PerturbSeqAnalysis


class TestPerturbSeqAnalysis(unittest.TestCase):

    def setUp(self):
        factory = RawCX2NetworkFactory()
        hier_net = factory.get_cx2network(os.path.join(os.path.dirname(__file__), 'data', 'hierarchy_perturb_test.cx2'))
        self.perturb_table = pd.read_table(os.path.join(os.path.dirname(__file__), 'data', 'sample_perturb_data.csv'),
                                           sep=',', index_col=0)
        self.analysis_obj = PerturbSeqAnalysis(hier_net)

    def test_get_heatmap_for_given_hierarchy_system(self):
        r_data = self.analysis_obj.get_heatmap_for_given_hierarchy_system(72, self.perturb_table)
        self.assertEqual(10, len(r_data))
        self.assertAlmostEqual(r_data.iloc[0, 0], -1.57, delta=0.01)

    def test_get_root_gene_pair_similarities(self):
        r_value = self.analysis_obj.get_root_gene_pair_similarities()
        self.assertEqual(5147, len(r_value))
        self.assertEqual(1, r_value.iloc[1, 2])

    def test_get_root_overlapping_pair_similarities(self):
        root_pairs = self.analysis_obj.get_root_gene_pair_similarities()
        r_val1, r_val2 = self.analysis_obj.get_root_overlapping_pair_similarities(root_pairs, self.perturb_table)
        self.assertEqual(len(r_val1), 1233)
        self.assertEqual(len(r_val2), 1233)
        self.assertAlmostEqual(r_val1.loc['ESF1', 'TMA16'], 0.65, delta=0.01)

    def test_get_root_functional_data_similarity(self):
        root_pairs = self.analysis_obj.get_root_gene_pair_similarities()
        r_val1, r_val2 = self.analysis_obj.get_root_overlapping_pair_similarities(root_pairs, self.perturb_table)
        result = self.analysis_obj.get_root_functional_data_similarity(r_val1, r_val2)
        self.assertEqual(len(result), 426402)

    def test_get_cluster_similarity(self):
        root_pairs = self.analysis_obj.get_root_gene_pair_similarities()
        r_val1, r_val2 = self.analysis_obj.get_root_overlapping_pair_similarities(root_pairs, self.perturb_table)
        result = self.analysis_obj.get_cluster_similarity(r_val1, 72)
        self.assertEqual(len(result), 45)
        self.assertAlmostEqual(result[0], 0.83, delta=0.01)

    def test_compare_cluster_root_similarities(self):
        root_pairs = self.analysis_obj.get_root_gene_pair_similarities()
        r_val1, r_val2 = self.analysis_obj.get_root_overlapping_pair_similarities(root_pairs, self.perturb_table)
        sim_root = self.analysis_obj.get_root_functional_data_similarity(r_val1, r_val2)
        sim_cluster = self.analysis_obj.get_cluster_similarity(r_val1, 72)
        stat, p_value = self.analysis_obj.compare_cluster_root_similarities(sim_cluster, sim_root)
        self.assertAlmostEqual(stat, 7.49, delta=0.01)
        self.assertAlmostEqual(p_value, 3.35 * 10 ** -14, delta=0.01)




if __name__ == '__main__':
    unittest.main()
