#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `OllamaCommandLineGeneSetAgent` ."""

import os
import tempfile
import shutil
import unittest
from unittest.mock import patch


from cellmaps_hierarchyeval.analysis import GenesetAgent
from cellmaps_hierarchyeval.analysis import OllamaCommandLineGeneSetAgent
from cellmaps_hierarchyeval.exceptions import CellmapshierarchyevalError


class TestOllamaCommandLineGeneSetAgent(unittest.TestCase):
    """Tests for `OllamaCommandLineGeneSetAgent` ."""

    def test_update_prompt_with_gene_set(self):
        agent = OllamaCommandLineGeneSetAgent(prompt='Hello {' +
                                                     GenesetAgent.GENE_SET_TOKEN +
                                                     '}\n'
                                                     'well\n')
        new_prompt = agent._update_prompt_with_gene_set(gene_names=['gene1',
                                                                    'gene2'])
        self.assertEqual('Hello gene1,gene2\nwell\n', new_prompt)

    def test_get_attribute_name_prefix(self):
        agent = OllamaCommandLineGeneSetAgent(prompt=None)
        self.assertEqual('ollama_llama2:latest::', agent.get_attribute_name_prefix())

    def test_prompt_not_set(self):
        agent = OllamaCommandLineGeneSetAgent(prompt=None)

        self.assertTrue('{GENE_SET}' in agent.get_prompt())

    def test_annotate_gene_set_success(self):
        with patch.object(OllamaCommandLineGeneSetAgent, '_run_cmd',
                          return_value=(0, 'Process: someproc\n'
                                           'Confidence Score: 0.50\n'
                                           'some output', '')):
            agent = OllamaCommandLineGeneSetAgent(prompt=None)
            res = agent.annotate_gene_set(['gene1', 'gene2'])
            self.assertEqual(res, ('someproc', '0.50',
                                   'Process: someproc\nConfidence Score: '
                                   '0.50\nsome output'))





