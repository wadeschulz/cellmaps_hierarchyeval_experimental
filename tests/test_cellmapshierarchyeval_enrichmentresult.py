import unittest

from cellmaps_hierarchyeval.runner import EnrichmentResult


class TestEnrichmentResult(unittest.TestCase):

    def setUp(self):
        self.enrichment_result = EnrichmentResult(
            term="Pathway A",
            pval=0.05,
            jaccard_index=0.3,
            overlap_genes=["gene1", "gene2", "gene3"]
        )

    def test_initialization(self):
        self.assertEqual(self.enrichment_result.term, "Pathway A")
        self.assertEqual(self.enrichment_result.pval, 0.05)
        self.assertEqual(self.enrichment_result.jaccard_index, 0.3)
        self.assertEqual(self.enrichment_result.overlap_genes, ["gene1", "gene2", "gene3"])
        self.assertIsNone(self.enrichment_result.description)
        self.assertIsNone(self.enrichment_result.adjusted_pval)
        self.assertFalse(self.enrichment_result.accepted)

    def test_set_adjusted_pval(self):
        self.enrichment_result.set_adjusted_pval(0.01)
        self.assertEqual(self.enrichment_result.adjusted_pval, 0.01)

    def test_set_accepted_conditions_fulfilled(self):
        self.enrichment_result.set_adjusted_pval(0.01)
        self.enrichment_result.set_accepted(min_jaccard_index=0.25, max_fdr=0.05)
        self.assertTrue(self.enrichment_result.accepted)

    def test_set_accepted_jaccard_below_threshold(self):
        self.enrichment_result.set_adjusted_pval(0.01)
        self.enrichment_result.set_accepted(min_jaccard_index=0.4, max_fdr=0.05)
        self.assertFalse(self.enrichment_result.accepted)

    def test_set_accepted_fdr_above_threshold(self):
        self.enrichment_result.set_adjusted_pval(0.5)
        self.enrichment_result.set_accepted(min_jaccard_index=0.25, max_fdr=0.02)
        self.assertFalse(self.enrichment_result.accepted)

    def test_set_accepted_jaccard_is_one(self):
        self.enrichment_result.set_adjusted_pval(0.01)
        self.enrichment_result.jaccard_index = 1
        self.enrichment_result.set_accepted(min_jaccard_index=0.25, max_fdr=0.05)
        self.assertTrue(self.enrichment_result.accepted)

    def test_set_description(self):
        self.enrichment_result.set_description("This is a test description.")
        self.assertEqual(self.enrichment_result.description, "This is a test description.")


if __name__ == '__main__':
    unittest.main()
