#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `cellmaps_hierarchyeval` package."""

import os
import tempfile
import shutil
import unittest

import ndex2
from cellmaps_utils import constants
from cellmaps_utils.provenance import ProvenanceUtil

from cellmaps_hierarchyeval.runner import CellmapshierarchyevalRunner


class TestCellmapshierarchyevalrunner(unittest.TestCase):
    """Tests for `cellmaps_hierarchyeval` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

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

            error_log_file = os.path.join(outdir, 'error.log')
            self.assertEqual(0, os.path.getsize(error_log_file))
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
