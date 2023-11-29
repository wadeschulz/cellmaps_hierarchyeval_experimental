#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `cellmaps_hierarchyeval` package."""

import os
import tempfile
import shutil

import unittest
from cellmaps_hierarchyeval import cellmaps_hierarchyevalcmd


class TestCellmapshierarchyeval(unittest.TestCase):
    """Tests for `cellmaps_hierarchyeval` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_parse_arguments(self):
        """Tests parse arguments"""
        res = cellmaps_hierarchyevalcmd._parse_arguments('hi',
                                                         ['outdir',
                                                          cellmaps_hierarchyevalcmd.HIERARCHYDIR,
                                                          'foox'])

        self.assertEqual(res.verbose, 1)
        self.assertEqual(res.logconf, None)
        self.assertEqual(res.outdir, 'outdir')
        self.assertEqual(res.hierarchy_dir, 'foox')

        someargs = ['-vv', '--logconf', 'hi', 'resdir',
                    cellmaps_hierarchyevalcmd.HIERARCHYDIR,
                    'foo']
        res = cellmaps_hierarchyevalcmd._parse_arguments('hi', someargs)

        self.assertEqual(res.verbose, 3)
        self.assertEqual(res.logconf, 'hi')
        self.assertEqual(res.hierarchy_dir, 'foo')
        self.assertEqual(res.outdir, 'resdir')

    def test_main(self):
        """Tests main function"""

        # try where loading config is successful
        try:
            temp_dir = tempfile.mkdtemp()
            res = cellmaps_hierarchyevalcmd.main(['myprog.py',
                                                  'outdir',
                                                  cellmaps_hierarchyevalcmd.HIERARCHYDIR,
                                                  'foox'])
            self.assertEqual(res, 2)
        finally:
            shutil.rmtree(temp_dir)
