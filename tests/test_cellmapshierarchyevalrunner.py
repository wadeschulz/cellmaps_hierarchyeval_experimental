#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `cellmaps_hierarchyeval` package."""

import os
import tempfile
import shutil
import unittest
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

    def test_get_hierarchy_file(self):
        runner = CellmapshierarchyevalRunner('outdir',
                                             hierarchy_dir='/hier')
        self.assertIsNotNone('/hier/hierarchy.cx',
                             runner.get_hierarchy_input_file())

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
                                                 input_data_dict={})
            self.assertEqual(0, runner.run())

            # todo finish test
            self.assertEqual(1, 2)
        finally:
            shutil.rmtree(temp_dir)
