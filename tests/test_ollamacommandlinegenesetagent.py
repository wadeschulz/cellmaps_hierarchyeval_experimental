#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `OllamaCommandLineGeneSetAgent` ."""

import os
import tempfile
import shutil
import unittest

from cellmaps_hierarchyeval.analysis import OllamaCommandLineGeneSetAgent
from cellmaps_hierarchyeval.exceptions import CellmapshierarchyevalError


class TestOllamaCommandLineGeneSetAgent(unittest.TestCase):
    """Tests for `OllamaCommandLineGeneSetAgent` ."""

    def test_update_prompt_with_gene_set(self):
        agent = OllamaCommandLineGeneSetAgent(prompt='Hello @@GENE_SET@@\n'
                                                     'well\n')
        new_prompt = agent._update_prompt_with_gene_set(gene_names=['gene1',
                                                                    'gene2'])
        self.assertEqual('Hello gene1,gene2\nwell\n', new_prompt)

    def test_prompt_not_set(self):
        agent = OllamaCommandLineGeneSetAgent(prompt=None)

        self.assertTrue('@@GENE_SET@@' in agent.get_prompt())



