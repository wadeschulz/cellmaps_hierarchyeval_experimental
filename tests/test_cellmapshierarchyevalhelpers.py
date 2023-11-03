#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `cellmaps_hierarchyeval` package."""

import os
import json
import unittest
from unittest.mock import MagicMock

from ndex2 import NiceCXNetwork
from ndex2.cx2 import CX2Network
from cellmaps_utils import constants

from cellmaps_hierarchyeval.runner import BaseNetworkHelper, CX2NetworkHelper, \
    NiceCXNetworkHelper


class TestCellmapshierarchyevalHelpers(unittest.TestCase):
    """Tests for `cellmaps_hierarchyeval` package."""

    def setUp(self):
        """Set up test fixtures, if any."""
        self.hierarchy_cx_path = os.path.join(os.path.dirname(__file__), 'data',
                                              'hierarchy.cx')
        self.hierarchy_cx2_path = os.path.join(os.path.dirname(__file__), 'data',
                                               'hierarchy.cx2')
        self.cx_hierarchy_helper = NiceCXNetworkHelper(self.hierarchy_cx_path)
        self.cx2_hierarchy_helper = CX2NetworkHelper(self.hierarchy_cx2_path)
        self.test_file = os.path.join(os.path.dirname(__file__), 'data',
                                      'test_hierarchy_output.cx2')

    def tearDown(self):
        """Tear down test fixtures, if any."""
        # This will run after each test method
        if os.path.exists(self.test_file):
            os.remove(self.test_file)

    def test_get_hierarchy_input_file(self):
        helper = BaseNetworkHelper(self.hierarchy_cx_path)
        self.assertEqual(helper.get_hierarchy_input_file(), self.hierarchy_cx_path)

    def test_get_hierarchy_cx(self):
        self.assertTrue(isinstance(self.cx_hierarchy_helper.get_hierarchy(), NiceCXNetwork))

    def test_get_hierarchy_cx2(self):
        self.assertTrue(isinstance(self.cx2_hierarchy_helper.get_hierarchy(), CX2Network))

    def test_get_suffix_cx(self):
        self.assertEqual(self.cx_hierarchy_helper.get_suffix(), constants.CX_SUFFIX)

    def test_get_suffix_cx2(self):
        self.assertEqual(self.cx2_hierarchy_helper.get_suffix(), constants.CX2_SUFFIX)

    def test_get_format_cx(self):
        self.assertEqual(self.cx_hierarchy_helper.get_format(), 'CX')

    def test_get_format_cx2(self):
        self.assertEqual(self.cx2_hierarchy_helper.get_format(), 'CX2')

    def test_dump_to_file_cx(self):
        hierarchy = self.cx_hierarchy_helper.get_hierarchy()
        self.cx_hierarchy_helper.dump_to_file(hierarchy, self.test_file)
        self.assertTrue(os.path.exists(self.test_file))

        with open(self.test_file, 'r') as f:
            file_contents = json.load(f)
            self.assertEqual(file_contents, hierarchy.to_cx())

    def test_dump_to_file_cx2(self):
        hierarchy = self.cx2_hierarchy_helper.get_hierarchy()
        self.cx2_hierarchy_helper.dump_to_file(hierarchy, self.test_file)
        self.assertTrue(os.path.exists(self.test_file))

        with open(self.test_file, 'r') as f:
            file_contents = json.load(f)
            self.assertEqual(file_contents, hierarchy.to_cx2())

    def test_get_hierarchy_real_ids_cx(self):
        hierarchy = self.cx_hierarchy_helper.get_hierarchy()
        hierarchy_size = len(hierarchy.get_nodes())
        real_ids = self.cx_hierarchy_helper.get_hierarchy_real_ids(hierarchy, hierarchy_size)
        self.assertEqual(list(range(hierarchy_size)), real_ids)

    def test_get_hierarchy_real_ids_cx2(self):
        hierarchy = self.cx2_hierarchy_helper.get_hierarchy()
        hierarchy_size = len(hierarchy.get_nodes())
        real_ids = self.cx2_hierarchy_helper.get_hierarchy_real_ids(hierarchy, hierarchy_size)
        self.assertEqual(list(hierarchy.get_nodes().keys()), real_ids)

    def test_get_node_genes_cx(self):
        mock_hierarchy = MagicMock()
        mock_hierarchy.get_node_attribute.return_value = {'v': 'gene1 gene2 gene3'}
        genes = self.cx_hierarchy_helper.get_node_genes(mock_hierarchy, 1)
        self.assertEqual(genes, {'gene1', 'gene2', 'gene3'})
        mock_hierarchy.get_node_attribute.assert_called_with(1, 'CD_MemberList')

    def test_get_node_genes_cx2(self):
        mock_node = {
            'v': {
                'CD_MemberList': 'gene1 gene2 gene3'
            }
        }
        genes = self.cx2_hierarchy_helper.get_node_genes(None, mock_node)
        self.assertEqual(genes, ['gene1', 'gene2', 'gene3'])

    def test_get_nodes_cx(self):
        mock_hierarchy = MagicMock()
        mock_hierarchy.nodes = {'node1': {}, 'node2': {}}
        nodes = self.cx_hierarchy_helper.get_nodes(mock_hierarchy)
        self.assertEqual(nodes, {'node1': {}, 'node2': {}})

    def test_get_nodes_cx2(self):
        mock_hierarchy = MagicMock()
        mock_hierarchy.get_nodes.return_value = {'node1': {}, 'node2': {}}
        nodes = self.cx2_hierarchy_helper.get_nodes(mock_hierarchy)
        self.assertEqual(nodes, {'node1': {}, 'node2': {}})
        mock_hierarchy.get_nodes.assert_called_once()

    def test_write_as_nodelist_cx(self):
        mock_hierarchy = MagicMock()
        mock_hierarchy.get_node_attributes.return_value = [{'n': 'attribute1', 'v': 'value1'},
                                                           {'n': 'attribute2', 'v': 'value2'}]
        mock_hierarchy.get_nodes.return_value = iter([(1, {'n': 'node1'}), (2, {'n': 'node2'})])
        test_output_path = os.path.join(os.path.dirname(__file__), 'data', 'test_nodelist_output.txt')
        try:
            self.cx_hierarchy_helper.write_as_nodelist(mock_hierarchy, test_output_path)
            self.assertTrue(os.path.exists(test_output_path))
            with open(test_output_path, 'r') as file:
                lines = file.readlines()
                self.assertEqual(lines[0].strip(), 'Name\tattribute1\tattribute2')
                self.assertEqual(lines[1].strip(), 'node1\tvalue1\tvalue2')
                self.assertEqual(lines[2].strip(), 'node2\tvalue1\tvalue2')
        finally:
            if os.path.exists(test_output_path):
                os.remove(test_output_path)

    def test_write_as_nodelist_cx2(self):
        mock_hierarchy = MagicMock()
        mock_hierarchy.get_attribute_declarations.return_value = {"nodes": {
            "attribute1": {
                "d": "string"
            },
            "attribute2": {
                "d": "string"
            }}}
        mock_hierarchy.get_nodes.return_value = {
            '1': {'v': {'attribute1': 'value1', 'attribute2': 'value2'}},
            '2': {'v': {'attribute1': 'value3', 'attribute2': 'value4'}}
        }
        test_output_path = os.path.join(os.path.dirname(__file__), 'data', 'test_nodelist_output.txt')

        try:
            self.cx2_hierarchy_helper.write_as_nodelist(mock_hierarchy, test_output_path)
            self.assertTrue(os.path.exists(test_output_path))
            with open(test_output_path, 'r') as file:
                lines = file.readlines()
                self.assertEqual(lines[0].strip(), 'attribute1\tattribute2')
                self.assertEqual(lines[1].strip(), 'value1\tvalue2')
                self.assertEqual(lines[2].strip(), 'value3\tvalue4')
        finally:
            if os.path.exists(test_output_path):
                os.remove(test_output_path)


