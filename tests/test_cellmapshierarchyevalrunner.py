#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `cellmaps_hierarchyeval` package."""

import os
import tempfile
import shutil
import unittest
from unittest.mock import patch, Mock, MagicMock

import ndex2
from ndex2.cx2 import CX2Network
from cellmaps_utils import constants
from cellmaps_utils.provenance import ProvenanceUtil
from requests import RequestException

from cellmaps_hierarchyeval.analysis import FakeGeneSetAgent
from cellmaps_hierarchyeval.exceptions import CellmapshierarchyevalError
from cellmaps_hierarchyeval.runner import CellmapshierarchyevalRunner, NiceCXNetworkHelper, CX2NetworkHelper

@unittest.skipIf(os.getenv('CELLMAPS_HIERARCHYEVAL_BAD_INTERNET') is not None, 'Too slow internet')
class TestCellmapshierarchyevalrunner(unittest.TestCase):
    """Tests for `cellmaps_hierarchyeval` package."""

    def setUp(self):
        """Set up test fixtures, if any."""
        self.runner = CellmapshierarchyevalRunner('outdir')
        self.runner._hierarchy_dir = "mock_hierarchy_dir"

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def get_4nodehierarchy(self):
        return os.path.join(os.path.dirname(__file__), 'data',
                            '4nodehierarchy.cx')

    def test_constructor(self):
        """Tests constructor"""
        runner = CellmapshierarchyevalRunner('outdir')
        self.assertIsNotNone(runner)

    def test_4nodehierarchy(self):
        temp_dir = tempfile.mkdtemp()
        try:
            hier_dir = os.path.join(temp_dir, '4.hierarchy')
            os.makedirs(hier_dir, mode=0o755)
            prov = ProvenanceUtil()
            prov.register_rocrate(hier_dir, name='hierarchy1',
                                  organization_name='hierarchy org',
                                  project_name='hierarchy project')

            srchier_file = os.path.join(hier_dir,
                                        constants.HIERARCHY_NETWORK_PREFIX +
                                        constants.CX_SUFFIX)
            shutil.copy(self.get_4nodehierarchy(),
                        srchier_file)

            outdir = os.path.join(temp_dir, 'outdir')
            runner = CellmapshierarchyevalRunner(outdir,
                                                 hierarchy_dir=hier_dir,
                                                 skip_logging=False,
                                                 input_data_dict={})
            self.assertEqual(0, runner.run())

            hierarchy_file = runner.get_annotated_hierarchy_dest_file()
            hier_net = ndex2.create_nice_cx_from_file(hierarchy_file)
            self.assertEqual('test 4node hierarchy', hier_net.get_name())

            # verify node attributes have been added
            for node_id, node_obj in hier_net.get_nodes():
                for db_name in ['GO_CC', 'HPA', 'CORUM']:
                    for db_attr in ['_terms', '_descriptions', '_jaccard_indexes',
                                    '_overlap_genes']:
                        if (db_name == 'HPA' or db_name == 'CORUM') and db_attr == '_descriptions':
                            # CORUM & HPA do not have descriptions
                            continue
                        n_attr = hier_net.get_node_attribute(node_id, db_name + db_attr)
                        self.assertIsNotNone(n_attr, str(node_obj) + ' is none for ' + db_name + db_attr)

        finally:
            shutil.rmtree(temp_dir)

    @patch('os.path.exists')
    def test_initialize_hierarchy_helper_with_cx(self, mock_exists):
        mock_exists.side_effect = lambda path: path.endswith(constants.CX_SUFFIX)
        self.runner.initialize_hierarchy_helper()
        self.assertIsInstance(self.runner._hierarchy_helper, NiceCXNetworkHelper)

    @patch('os.path.exists')
    def test_initialize_hierarchy_helper_with_cx2(self, mock_exists):
        mock_exists.side_effect = lambda path: path.endswith(constants.CX2_SUFFIX)
        self.runner.initialize_hierarchy_helper()
        self.assertIsInstance(self.runner._hierarchy_helper, CX2NetworkHelper)

    @patch('os.path.exists')
    def test_initialize_hierarchy_helper_with_no_files(self, mock_exists):
        mock_exists.return_value = False
        with self.assertRaises(CellmapshierarchyevalError) as err:
            self.runner.initialize_hierarchy_helper()
        self.assertTrue(f"Input directory '{self.runner._hierarchy_dir}' does not contain "
                        f"neither cx nor cx2 files." in str(err.exception))

    @patch('ndex2.create_nice_cx_from_server')
    def test_get_network_failure(self, mock_create_nice_cx):
        mock_response = Mock()
        mock_response.text = "Server error"
        mock_create_nice_cx.side_effect = [RequestException(response=mock_response),
                                           RequestException(response=mock_response),
                                           RequestException(response=mock_response)]

        your_instance = CellmapshierarchyevalRunner('foo')

        with self.assertRaises(CellmapshierarchyevalError) as context:
            your_instance._get_network_from_server(uuid="some_uuid", retry_wait=0)

        self.assertIn('3 attempts to get network some_uuid failed', str(context.exception))

    def test_add_empty_attr_to_hierarchy(self):
        mock_hierarchy = MagicMock()
        mock_terms = MagicMock()
        mock_terms.term_name = "TestTerm"
        node_ids = [1, 2, 3]
        self.runner._add_empty_attr_to_hierarchy(mock_hierarchy, mock_terms, node_ids=node_ids)
        expected_call_count = len(node_ids) * 5

        actual_call_count = mock_hierarchy.set_node_attribute.call_count
        self.assertEqual(expected_call_count, actual_call_count,
                         f"Expected set_node_attribute to be called {expected_call_count} times, got {actual_call_count}")

    def test_annotate_hierarchy_with_geneset_annotators(self):
        gsai = MagicMock()
        gsai.get_attribute_name_prefix = MagicMock(return_val='foo::')
        gsai.annotate_gene_set = MagicMock()
        temp_dir = tempfile.mkdtemp()
        try:
            mockannotator = MagicMock()
            mockannotator.annotate_hierarchy = MagicMock()
            mockannotator.set_hierarchy_helper = MagicMock()
            runner = CellmapshierarchyevalRunner(os.path.join(temp_dir, 'foo'),
                                                 geneset_agents=[gsai],
                                                 geneset_annotator=mockannotator)

            hierarchy = CX2Network()
            runner._annotate_hierarchy_with_geneset_annotators(hierarchy=hierarchy)
            mockannotator.annotate_hierarchy.assert_called()
            mockannotator.set_hierarchy_helper.assert_called()
        finally:
            shutil.rmtree(temp_dir)

    def test_four_node_hierarchy(self):
        temp_dir = tempfile.mkdtemp()
        try:
            agent = FakeGeneSetAgent()
            runner = CellmapshierarchyevalRunner(os.path.join(temp_dir, 'foo'),
                                                 geneset_agents=[agent])

            hierhelper = CX2NetworkHelper(os.path.join(os.path.dirname(__file__),
                                                       'data', 'hierarchy.cx2'))
            hierarchy = hierhelper.get_hierarchy()
            runner._hierarchy_helper = hierhelper
            runner._annotate_hierarchy_with_geneset_annotators(hierarchy=hierarchy)
            attributes_with_process = []
            attributes_with_confidence = []
            attributes_with_raw = []

            for node, attrs in hierarchy.get_nodes().items():
                for attr in attrs['v']:
                    if attr.startswith('fake') and attr.endswith('_process'):
                        attributes_with_process.append(attr)
                    elif attr.startswith('fake') and attr.endswith('_confidence'):
                        attributes_with_confidence.append(attr)
                    elif attr.startswith('fake') and attr.endswith('_raw'):
                        attributes_with_raw.append(attr)

            self.assertEqual(len(attributes_with_process), len(hierarchy.get_nodes()))
            self.assertEqual(len(attributes_with_confidence), len(hierarchy.get_nodes()))
            self.assertEqual(len(attributes_with_raw), len(hierarchy.get_nodes()))
        finally:
            shutil.rmtree(temp_dir)
