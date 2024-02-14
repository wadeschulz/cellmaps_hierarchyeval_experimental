import unittest
from unittest.mock import MagicMock

from cellmaps_hierarchyeval.runner import GO_EnrichmentTerms, CORUM_EnrichmentTerms, HPA_EnrichmentTerms, \
    HiDeF_EnrichmentTerms


class TestEnrichmentTerms(unittest.TestCase):

    def setUp(self):
        self.hierarchy_genes = ['gene1', 'gene2', 'gene3']
        self.terms = MagicMock()
        self.terms.get_node_attribute_value.side_effect = self.mock_get_node_attribute_value

    @staticmethod
    def mock_get_node_attribute_value(_, attribute):
        if attribute == 'description':
            return 'Description for Term1'
        elif attribute == 'genes':
            return 'gene1,gene2'
        elif attribute == 'subunits(Gene name)':
            return ['gene1', 'gene2']
        elif attribute == 'Main location':
            return ['Location1', 'Location2']
        return None

    def test_go_enrichment_terms(self):
        self.terms.get_nodes.return_value = [
            (1, {'n': 'Term1'})
        ]
        go_terms = GO_EnrichmentTerms(self.terms, 'GO Term', self.hierarchy_genes, min_comp_size=2)
        self.assertTrue('gene1' in go_terms.term_genes['Term1'])
        self.assertTrue('gene2' in go_terms.term_genes['Term1'])
        self.assertEqual(go_terms.term_description, {'Term1': 'Description for Term1'})

    def test_corum_enrichment_terms(self):
        self.terms.get_nodes.return_value = [
            (1, {'n': 'Term1'})
        ]
        corum_terms = CORUM_EnrichmentTerms(self.terms, 'CORUM Term', self.hierarchy_genes, min_comp_size=2)
        self.assertTrue('gene1' in corum_terms.term_genes['Term1'])
        self.assertTrue('gene2' in corum_terms.term_genes['Term1'])

    def test_hpa_enrichment_terms(self):
        self.terms.get_nodes.return_value = [
            (1, {'n': 'gene1'}),
            (2, {'n': 'gene2'}),
            (3, {'n': 'gene4'})
        ]
        hpa_terms = HPA_EnrichmentTerms(self.terms, 'HPA Term', self.hierarchy_genes, min_comp_size=2)
        self.assertTrue('gene1' in hpa_terms.term_genes['Location1'])
        self.assertTrue('gene2' in hpa_terms.term_genes['Location1'])
        self.assertTrue('gene1' in hpa_terms.term_genes['Location2'])
        self.assertTrue('gene2' in hpa_terms.term_genes['Location2'])
        self.assertTrue('gene4' not in hpa_terms.term_genes.get('Location1', []))
        self.assertTrue('gene4' not in hpa_terms.term_genes.get('Location2', []))
        self.assertEqual(len(hpa_terms.term_genes['Location1']), 2)
        self.assertEqual(len(hpa_terms.term_genes['Location2']), 2)

    def test_hidef_enrichment_terms(self):
        self.terms.get_nodes.return_value = [
            (1, {'n': 'Term1'}),
            (2, {'n': 'Term2'}),
            (3, {'n': 'Term3'})
        ]
        self.terms.get_node_attribute_value.side_effect = lambda node, attribute: {
            ('Term1', 'CD_MemberList'): 'gene1 gene2 gene3',
            ('Term2', 'CD_MemberList'): 'gene3 gene4',
            ('Term3', 'CD_MemberList'): 'gene5',
        }.get((node['n'], attribute), None)

        hidef_terms = HiDeF_EnrichmentTerms(self.terms, 'HiDeF Term', self.hierarchy_genes, min_comp_size=3)

        self.assertIn('Term1', hidef_terms.term_genes)
        self.assertNotIn('Term2', hidef_terms.term_genes)
        self.assertNotIn('Term3', hidef_terms.term_genes)

        expected_genes = set(['gene1', 'gene2', 'gene3'])
        self.assertEqual(hidef_terms.all_term_genes, expected_genes)


if __name__ == '__main__':
    unittest.main()
