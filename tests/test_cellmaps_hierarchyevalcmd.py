#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `cellmaps_hierarchyeval` package."""

import os
import tempfile
import shutil

import unittest
from cellmaps_hierarchyeval import cellmaps_hierarchyevalcmd
from cellmaps_hierarchyeval.analysis import FakeGeneSetAgent
from cellmaps_hierarchyeval.analysis import OllamaRestServiceGenesetAgent
from cellmaps_hierarchyeval.analysis import OllamaCommandLineGeneSetAgent


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

    def test_get_ollama_geneset_agents_no_prompts(self):

        res = cellmaps_hierarchyevalcmd.get_ollama_geneset_agents(ollama_prompts=None)
        self.assertIsNone(res)

    def test_get_ollama_geneset_agents_service_prompts(self):
        temp_dir = tempfile.mkdtemp()
        try:
            prompt_file = os.path.join(temp_dir, 'myprompt.txt')
            with open(prompt_file, 'w') as f:
                f.write('my prompt')
            o_prompts = ['fake', 'modela,' + prompt_file,
                         'modelb,a prompt']
            res = cellmaps_hierarchyevalcmd.get_ollama_geneset_agents(ollama='http://foo',
                                                                      ollama_prompts=o_prompts)
            self.assertEqual(3, len(res))
            self.assertTrue(isinstance(res[0], FakeGeneSetAgent))
            self.assertTrue(isinstance(res[1], OllamaRestServiceGenesetAgent))
            self.assertEqual(res[1].get_prompt(), 'my prompt')
            self.assertTrue(isinstance(res[2], OllamaRestServiceGenesetAgent))
            self.assertEqual(res[2].get_prompt(), 'a prompt')

        finally:
            shutil.rmtree(temp_dir)

    def test_get_ollama_geneset_agents_commandline_prompts(self):
        temp_dir = tempfile.mkdtemp()
        try:
            prompt_file = os.path.join(temp_dir, 'myprompt.txt')
            with open(prompt_file, 'w') as f:
                f.write('my prompt')
            o_prompts = ['fake', 'modela,' + prompt_file,
                         'modelb,a prompt']
            res = cellmaps_hierarchyevalcmd.get_ollama_geneset_agents(ollama='/bin/ollama',
                                                                      ollama_prompts=o_prompts)
            self.assertEqual(3, len(res))
            self.assertTrue(isinstance(res[0], FakeGeneSetAgent))
            self.assertTrue(isinstance(res[1], OllamaCommandLineGeneSetAgent))
            self.assertEqual(res[1].get_prompt(), 'my prompt')
            self.assertTrue(isinstance(res[2], OllamaCommandLineGeneSetAgent))
            self.assertEqual(res[2].get_prompt(), 'a prompt')

        finally:
            shutil.rmtree(temp_dir)
