#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `OllamaCommandLineGeneSetAgent` ."""

import os
import tempfile
import shutil
import unittest
from unittest.mock import patch, Mock, MagicMock

import ndex2
from cellmaps_utils import constants
from cellmaps_utils.provenance import ProvenanceUtil
from requests import RequestException

from cellmaps_hierarchyeval.exceptions import CellmapshierarchyevalError
from cellmaps_hierarchyeval.runner import CellmapshierarchyevalRunner, NiceCXNetworkHelper, CX2NetworkHelper


class TestOllamaCommandLineGeneSetAgent(unittest.TestCase):
    """Tests for `OllamaCommandLineGeneSetAgent` ."""

    def test_update_promt_with_gene_set(self):




