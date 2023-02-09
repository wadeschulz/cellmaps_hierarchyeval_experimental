#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Integration Tests for `cellmaps_hierarchyeval` package."""

import os

import unittest
from cellmaps_hierarchyeval import cellmaps_hierarchyevalcmd

SKIP_REASON = 'CELLMAPS_HIERARCHYEVAL_INTEGRATION_TEST ' \
              'environment variable not set, cannot run integration ' \
              'tests'

@unittest.skipUnless(os.getenv('CELLMAPS_HIERARCHYEVAL_INTEGRATION_TEST') is not None, SKIP_REASON)
class TestIntegrationCellmaps_hierarchyeval(unittest.TestCase):
    """Tests for `cellmaps_hierarchyeval` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_something(self):
        """Tests parse arguments"""
        self.assertEqual(1, 1)
