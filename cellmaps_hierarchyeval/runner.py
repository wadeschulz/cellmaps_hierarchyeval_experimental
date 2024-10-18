#! /usr/bin/env python
import os
import logging
import shutil
import time
import json
import numpy as np
from datetime import date
from tqdm import tqdm

from requests import RequestException, JSONDecodeError
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import warnings
import ndex2
from ndex2.cx2 import CX2Network
from cellmaps_utils import constants
from cellmaps_utils import logutils
from cellmaps_utils.provenance import ProvenanceUtil
import cellmaps_hierarchyeval
from cellmaps_hierarchyeval.exceptions import CellmapshierarchyevalError

logger = logging.getLogger(__name__)


class EnrichmentTerms(object):
    """
    Base class for implementations that generate
    term databases for enrichment (i.e., HPA, CORUM, GO)
    """

    def __init__(self, terms=None, term_name=None,
                 hierarchy_genes=None, min_comp_size=4):
        """
        Constructor

        :param terms: The terms to be processed.
        :type terms: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork` or None
        :param term_name: Name of the term.
        :type term_name: str or None
        :param hierarchy_genes: Genes in the hierarchy.
        :type hierarchy_genes: list or None
        :param min_comp_size: Minimum number of genes in a term for it to be considered.
        :type min_comp_size: int
        """
        self.terms = terms
        self.term_name = term_name
        self.hierarchy_genes = hierarchy_genes
        self.min_comp_size = min_comp_size
        self.term_genes = None
        self.term_description = None


class GO_EnrichmentTerms(EnrichmentTerms):
    """
    This class extends the EnrichmentTerms class to handle terms specific to Gene Ontology (GO).
    """

    def __init__(self, terms=None, term_name=None,
                 hierarchy_genes=None, min_comp_size=4):
        """
        Constructor. Sets the parameters and initializes the term genes and term description.

        :param terms: The terms to be processed.
        :type terms: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork` or None
        :param term_name: Name of the term.
        :type term_name: str or None
        :param hierarchy_genes: Genes in the hierarchy.
        :type hierarchy_genes: list or None
        :param min_comp_size: Minimum number of genes in a term for it to be considered.
        :type min_comp_size: int
        """
        super().__init__(terms=terms, term_name=term_name,
                         hierarchy_genes=hierarchy_genes,
                         min_comp_size=min_comp_size)
        self.term_genes, self.all_term_genes = self._get_term_genes(terms)
        self.term_description = self._get_term_description(terms)

    def _get_term_genes(self, terms):
        """
        Retrieves genes for given GO terms and filters out terms with genes
        below the minimum term size.

        :param terms: The GO terms for which genes are to be retrieved.
        :type terms: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork`
        :return: Dictionary of GO terms mapped to their respective genes.
        :rtype: dict
        """
        term_genes_dict = {}
        all_term_genes = set()

        for node_id, node in terms.get_nodes():
            term = node.get('n')
            genes = terms.get_node_attribute_value(node, 'genes')
            if genes is None:
                continue
            genes = genes.split(',')
            genes = list(set(genes).intersection(set(self.hierarchy_genes)))
            all_term_genes.update(genes)
            if len(genes) < self.min_comp_size:
                continue
            term_genes_dict[term] = genes
        return term_genes_dict, all_term_genes

    def _get_term_description(self, terms):
        """
        Retrieves descriptions for the given GO terms.

        :param terms: The GO terms for which descriptions are to be retrieved.
        :type terms: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork`
        :return: Dictionary of GO terms mapped to their respective descriptions.
        :rtype: dict
        """
        term_description = {}
        for node_id, node in terms.get_nodes():
            term = node.get('n')
            term_description[term] = terms.get_node_attribute_value(node,
                                                                    'description')
        return term_description


class HiDeF_EnrichmentTerms(EnrichmentTerms):
    """
    This class extends the EnrichmentTerms class to handle terms specific to HiDeF output.
    """

    def __init__(self, terms=None, term_name=None, hierarchy_genes=None,
                 min_comp_size=4):
        """
        Constructor. Sets the parameters and initializes the term genes.

        :param terms: The terms to be processed.
        :type terms: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork` or None
        :param term_name: Name of the term.
        :type term_name: str or None
        :param hierarchy_genes: Genes in the hierarchy.
        :type hierarchy_genes: list or None
        :param min_comp_size: Minimum number of genes in a term for it to be considered.
        :type min_comp_size: int
        """
        super().__init__(terms=terms, term_name=term_name,
                         hierarchy_genes=hierarchy_genes,
                         min_comp_size=min_comp_size)
        self.term_genes, self.all_term_genes = self._get_term_genes(terms)
        self.term_description = None

    def _get_term_genes(self, terms):
        """
        Retrieves genes for given terms, and filters out terms with genes below the minimum term size.

        :param terms: The terms for which genes are to be retrieved.
        :type terms: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork`
        :return: Dictionary of terms mapped to their respective genes.
        :rtype: dict
        """
        term_genes_dict = {}
        all_term_genes = set()

        for node_id, node in terms.get_nodes():
            genes = terms.get_node_attribute_value(node, 'CD_MemberList')
            if genes is None:
                continue
            genes = genes.split(' ')
            genes = list(set(genes).intersection(set(self.hierarchy_genes)))
            all_term_genes.update(genes)
            if len(genes) < self.min_comp_size:
                continue
            term = node.get('n')
            term_genes_dict[term] = genes
        return term_genes_dict, all_term_genes


class CORUM_EnrichmentTerms(EnrichmentTerms):
    """
    This class extends the EnrichmentTerms class to handle terms specific to CORUM.
    """

    def __init__(self, terms=None, term_name=None, hierarchy_genes=None,
                 min_comp_size=4):
        """
        Constructor. Sets the parameters and initializes the term genes.

        :param terms: The terms to be processed.
        :type terms: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork` or None
        :param term_name: Name of the term.
        :type term_name: str or None
        :param hierarchy_genes: Genes in the hierarchy.
        :type hierarchy_genes: list or None
        :param min_comp_size: Minimum number of genes in a term for it to be considered.
        :type min_comp_size: int
        """
        super().__init__(terms=terms, term_name=term_name,
                         hierarchy_genes=hierarchy_genes,
                         min_comp_size=min_comp_size)
        self.term_genes, self.all_term_genes = self._get_term_genes(terms)
        self.term_description = None

    def _get_term_genes(self, terms):
        """
        Retrieves genes for given terms, and filters out terms with genes below the minimum term size.

        :param terms: The terms for which genes are to be retrieved.
        :type terms: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork`
        :return: Dictionary of terms mapped to their respective genes.
        :rtype: dict
        """
        term_genes_dict = {}
        all_term_genes = set()

        for node_id, node in terms.get_nodes():
            genes = terms.get_node_attribute_value(node, 'subunits(Gene name)')
            if genes is None:
                continue
            genes = list(set(genes).intersection(set(self.hierarchy_genes)))
            all_term_genes.update(genes)
            if len(genes) < self.min_comp_size:
                continue
            term = node.get('n')
            term_genes_dict[term] = genes
        return term_genes_dict, all_term_genes


class HPA_EnrichmentTerms(EnrichmentTerms):
    """
    This class extends the EnrichmentTerms class to handle terms specific to the Human Protein Atlas (HPA).
    """

    def __init__(self, terms=None, term_name=None, hierarchy_genes=None,
                 min_comp_size=4):
        super().__init__(terms=terms, term_name=term_name,
                         hierarchy_genes=hierarchy_genes,
                         min_comp_size=min_comp_size)
        self.term_genes, self.all_term_genes = self._get_term_genes(terms)
        self.term_description = None

    def _get_term_genes(self, terms):
        """
        Retrieves genes for given HPA terms and builds a dictionary where keys are annotations
        and values are lists of genes corresponding to those annotations.

        :param terms: The terms for which genes are to be retrieved.
        :type terms: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork`
        :return: Dictionary of terms mapped to their respective genes.
        :rtype: dict
        """
        term_genes_dict = {}
        all_term_genes = set()

        for node_id, node in terms.get_nodes():
            node_name = node.get('n')
            if node_name not in self.hierarchy_genes:
                continue
            all_term_genes.add(node_name)
            for a in ['Main location', 'Additional location']:
                annotations = terms.get_node_attribute_value(node, a)
                if annotations is None:
                    continue
                for c in annotations:
                    if c in term_genes_dict:
                        term_genes_dict[c].append(node_name)
                    else:
                        term_genes_dict[c] = [node_name]
        return term_genes_dict, all_term_genes


class EnrichmentResult(object):
    """
    Base class for representing the results of enrichment analysis.
    It generates a hierarchy that is output in the CX format following the CDAPS style.
    """

    def __init__(self,
                 term=None,
                 pval=None,
                 jaccard_index=None,
                 overlap_genes=None):
        """
        Constructor

        :param term: The term name.
        :type term: str
        :param pval: P-value of the enrichment result.
        :type pval: float
        :param jaccard_index: Jaccard index of the enrichment result.
        :type jaccard_index: float
        :param overlap_genes: List of overlapping genes.
        :type overlap_genes: list
        """
        self.term = term
        self.description = None
        self.pval = pval
        self.jaccard_index = jaccard_index
        self.overlap_genes = overlap_genes
        self.adjusted_pval = None
        self.accepted = False

    def set_adjusted_pval(self, adjusted_pval):
        """
        Sets the adjusted p-value for the enrichment result.

        :param adjusted_pval: Adjusted p-value.
        :type adjusted_pval: float
        """
        self.adjusted_pval = adjusted_pval

    def set_accepted(self, min_jaccard_index, max_fdr):
        """
        Sets the accepted status of the enrichment result based on Jaccard index and FDR criteria.

        :param min_jaccard_index: Minimum required Jaccard index for the result to be accepted.
        :type min_jaccard_index: float
        :param max_fdr: Maximum allowed adjusted p-value (FDR) for the result to be accepted.
        :type max_fdr: float
        """
        if self.jaccard_index == 1:
            self.accepted = True
        elif (self.jaccard_index >= min_jaccard_index) & (self.adjusted_pval < max_fdr):
            self.accepted = True

    def set_description(self, description):
        """
        Sets the description of the enrichment term results.

        :param description: Description for the term results.
        :type description: str
        """
        self.description = description


class BaseNetworkHelper(object):
    """
    Base class for network helpers.
    """

    def __init__(self, hierarchy_path):
        """
        Constructor.

        :param hierarchy_path: File system path where the hierarchy network data is stored.
        :type hierarchy_path: str
        """
        self._hierarchy_path = hierarchy_path

    def get_hierarchy_input_file(self):
        """
        Creates file path prefix for hierarchy

        Example path: ``/tmp/foo/hierarchy``

        :return: Prefix path on filesystem where Hierarchy Network resides
        :rtype: str
        """
        return self._hierarchy_path


class CX2NetworkHelper(BaseNetworkHelper):
    """
    Helper class for CX2 network data manipulation that extends the
    BaseNetworkHelper class with CX2-specific logic.
    """

    def __init__(self, hierarchy_path):
        """
        Constructor.

        :param hierarchy_path: File system path where the CX2 hierarchy network data is stored.
        :type hierarchy_path: str
        """
        super().__init__(hierarchy_path)

    def get_hierarchy(self):
        """
        Create and return a CX2 network object from the hierarchy path.

        :return: An instance of the CX2Network class.
        :rtype: CX2Network
        """
        cx2_obj = CX2Network()
        cx2_obj.create_from_raw_cx2(self._hierarchy_path)
        return cx2_obj

    @staticmethod
    def get_suffix():
        """
        Get the file suffix associated with CX2 files.

        :return: The suffix for CX2 file types.
        :rtype: str
        """
        return constants.CX2_SUFFIX

    @staticmethod
    def get_format():
        """
        Get string format identifier for CX2 network data.

        :return: The format identifier for CX2.
        :rtype: str
        """
        return 'CX2'

    @staticmethod
    def dump_to_file(hierarchy, hierarchy_out_file):
        """
        Save the hierarchy to a CX2 formatted JSON file.

        :param hierarchy: The hierarchy to save.
        :type hierarchy: CX2Network
        :param hierarchy_out_file: The file path where the hierarchy should be written.
        :type hierarchy_out_file: str
        """
        with open(hierarchy_out_file, 'w') as f:
            json.dump(hierarchy.to_cx2(), f)

    @staticmethod
    def get_hierarchy_real_ids(hierarchy=None, hierarchy_size=None):
        """
        Retrieve the real identifiers of nodes within the hierarchy.

        :param hierarchy: The hierarchy from which to extract node IDs.
        :type hierarchy: CX2Network
        :param hierarchy_size: Not used, but specified for compatibility.
        :return: A list of node identifiers.
        :rtype: list
        """
        return [node_id for node_id in hierarchy.get_nodes().keys()]

    @staticmethod
    def get_node_genes(_, node=None):
        """
        Extract the gene identifiers from a given node.

        :param _: Placeholder, not used.
        :param node: The node from which to extract gene identifiers.
        :type node: dict
        :return: A list of gene identifiers.
        :rtype: list
        """
        return node.get('v', {}).get('CD_MemberList', '').split(' ')

    @staticmethod
    def get_nodes(hierarchy):
        """
        Retrieve the nodes from the hierarchy.

        :param hierarchy: The hierarchy from which to retrieve nodes.
        :type hierarchy: CX2Network
        :return: A dictionary of nodes.
        :rtype: dict
        """
        return hierarchy.get_nodes()

    @staticmethod
    def write_as_nodelist(hierarchy, dest_path):
        """
        Write the nodes of the hierarchy to a specified file path as a tab-delimited list.

        :param hierarchy: The hierarchy containing the nodes to write.
        :type hierarchy: CX2Network
        :param dest_path: The destination file path for the nodelist.
        :type dest_path: str
        """
        with open(dest_path, 'w') as f:
            # write headers
            attribute_declarations = hierarchy.get_attribute_declarations()
            node_attribute_names = attribute_declarations.get('nodes', [])
            for attribute_name in node_attribute_names:
                f.write(attribute_name + '\t')
            f.write('\n')

            # write node attributes
            nodes_dict = hierarchy.get_nodes()
            for node_id, node_data in nodes_dict.items():
                node_attributes = node_data["v"]
                for attribute_name in node_attribute_names:
                    value = node_attributes.get(attribute_name, "")
                    cleaned_value = str(value).replace('\n', ' ').replace('\t', ' ')
                    f.write(cleaned_value + '\t')
                f.write('\n')


class NiceCXNetworkHelper(BaseNetworkHelper):
    """
    Helper class for NiceCX network data manipulation that extends the
    BaseNetworkHelper class with CX-specific logic.
    """

    def __init__(self, hierarchy_path):
        """
        Constructor.

        :param hierarchy_path: File system path where the NiceCX hierarchy network data is stored.
        :type hierarchy_path: str
        """
        super().__init__(hierarchy_path)

    def get_hierarchy(self):
        """
        Create and return a NiceCXNetwork object from the hierarchy path.

        :return: An instance of the NiceCX network class.
        :rtype: ndex2.nice_cx_network.NiceCXNetwork
        """
        return ndex2.create_nice_cx_from_file(self._hierarchy_path)

    @staticmethod
    def get_suffix():
        """
        Get the file suffix associated with CX files.

        :return: The suffix for NiceCX file types.
        :rtype: str
        """
        return constants.CX_SUFFIX

    @staticmethod
    def get_format():
        """
        Get the string format identifier for CX data.

        :return: The format identifier for NiceCX.
        :rtype: str
        """
        return 'CX'

    @staticmethod
    def dump_to_file(hierarchy, hierarchy_out_file):
        """
        Save the hierarchy to a CX formatted JSON file.

        :param hierarchy: The hierarchy to save.
        :type hierarchy: ndex2.nice_cx_network.NiceCXNetwork
        :param hierarchy_out_file: The file path where the hierarchy should be written.
        :type hierarchy_out_file: str
        """
        with open(hierarchy_out_file, 'w') as f:
            json.dump(hierarchy.to_cx(), f)

    @staticmethod
    def get_hierarchy_real_ids(hierarchy=None, hierarchy_size=None):
        """
        Generate a list of real IDs for a given hierarchy size.

        :param hierarchy: Not used, provided for compatibility.
        :param hierarchy_size: The size of the hierarchy to generate IDs for.
        :type hierarchy_size: int
        :return: A list of sequential integers representing node IDs.
        :rtype: list
        """
        return list(range(hierarchy_size))

    @staticmethod
    def get_node_genes(hierarchy=None, node=None):
        """
        Extract the set of gene identifiers from a given node in the hierarchy.

        :param hierarchy: The hierarchy containing the node.
        :type hierarchy: ndex2.nice_cx_network.NiceCXNetwork
        :param node: The node from which to extract gene identifiers.
        :type node: int
        :return: A set of gene identifiers.
        :rtype: set
        """
        return set(hierarchy.get_node_attribute(node, 'CD_MemberList')['v'].split(' '))

    @staticmethod
    def get_nodes(hierarchy):
        """
        Retrieve the nodes from the hierarchy.

        :param hierarchy: The hierarchy from which to retrieve nodes.
        :type hierarchy: ndex2.nice_cx_network.NiceCXNetwork
        :return: A dictionary of nodes.
        :rtype: dict
        """
        return hierarchy.nodes

    @staticmethod
    def write_as_nodelist(hierarchy, dest_path):
        """
        Write the nodes of the hierarchy to a specified file path as a tab-delimited list.

        :param hierarchy: The hierarchy containing the nodes to write.
        :type hierarchy: ndex2.nice_cx_network.NiceCXNetwork
        :param dest_path: The destination file path for the nodelist.
        :type dest_path: str
        """
        with open(dest_path, 'w') as f:
            # write headers
            f.write('Name' + '\t')
            for a in hierarchy.get_node_attributes(0):
                f.write(str(a['n']) + '\t')
            f.write('\n')

            # write node attributes
            for node_id, node_obj in hierarchy.get_nodes():
                f.write(node_obj['n'] + '\t')
                for a in hierarchy.get_node_attributes(node_obj):
                    f.write(str(a['v']) + '\t')
                f.write('\n')


class GeneSetAgentAnnotator(object):
    """
    Annotates hierarchy with results from one or more
    :py:class:`~cellmaps_hierarchyeval.analysis.GeneSetAgent` objects
    """

    def __init__(self):
        """
        Constructor
        """
        self._hierarchy_helper = None
        self._min_comp_size = 4

    def set_hierarchy_helper(self, hierarchy_helper):
        """
        Sets HierarchyHelper

        :param hierarchy_helper:
        :return:
        """
        self._hierarchy_helper = hierarchy_helper

    def set_minimum_comparison_size(self, val):
        """
        Only examine genesets of size **val** or larger
        :param val:
        :type val: int
        """
        self._min_comp_size = val

    def annotate_hierarchy(self, geneset_agent=None,
                           hierarchy=None):
        """
        Annotates hierarchy with
        :py:class:`~cellmaps_hierarchyeval.analysis.GeneSetAgent`
        by adding new node attributes
        :param geneset_agent:
        :param hierarchy:
        :return:
        """
        for node_id, node in tqdm(self._hierarchy_helper.get_nodes(hierarchy).items(), desc='Assemblies'):
            gene_names = self._hierarchy_helper.get_node_genes(hierarchy, node)
            if gene_names is None or len(gene_names) == 0:
                logger.debug('No genes to analyze')
                hierarchy.set_node_attribute(node_id, f'{geneset_agent.get_attribute_name_prefix()}_process', '')
                hierarchy.set_node_attribute(node_id, f'{geneset_agent.get_attribute_name_prefix()}_confidence', '')
                hierarchy.set_node_attribute(node_id, f'{geneset_agent.get_attribute_name_prefix()}_raw', '')
                continue
            if len(gene_names) < self._min_comp_size:
                logger.debug('Skipping node: ' + str(node_id) +
                             ' has only ' + len(gene_names) +
                             '  which is below threshold of ' +
                             str(self._min_comp_size))
                continue
            proc_name, \
                confidence, \
                output = geneset_agent.annotate_gene_set(gene_names=gene_names)
            print('Proc name: ' + str(proc_name))
            print('confidence: ' + str(confidence))
            hierarchy.set_node_attribute(node_id, f'{geneset_agent.get_attribute_name_prefix()}_process', proc_name)
            hierarchy.set_node_attribute(node_id, f'{geneset_agent.get_attribute_name_prefix()}_confidence', confidence)
            hierarchy.set_node_attribute(node_id, f'{geneset_agent.get_attribute_name_prefix()}_raw', output)


class CellmapshierarchyevalRunner(object):
    """
    Class to run Hierarchy evaluation
    """
    MAX_FDR = 0.05
    MIN_JACCARD_INDEX = 0.1
    MIN_COMP_SIZE = 4
    CORUM = '633291aa-6e1d-11ef-a7fd-005056ae23aa'
    GO_CC = '6722d74d-6e20-11ef-a7fd-005056ae23aa'
    HPA = '68c2f2c0-6e20-11ef-a7fd-005056ae23aa'
    NDEX_SERVER = 'http://www.ndexbio.org'

    def __init__(self, outdir=None,
                 hierarchy_dir=None,
                 min_comp_size=MIN_COMP_SIZE,
                 max_fdr=MAX_FDR,
                 min_jaccard_index=MIN_JACCARD_INDEX,
                 corum=CORUM,
                 go_cc=GO_CC,
                 hpa=HPA,
                 ndex_server=NDEX_SERVER,
                 geneset_agents=None,
                 name=None,
                 organization_name=None,
                 project_name=None,
                 input_data_dict=None,
                 skip_term_enrichment=False,
                 skip_logging=True,
                 provenance_utils=ProvenanceUtil(),
                 geneset_annotator=GeneSetAgentAnnotator(),
                 provenance=None):
        """
        Constructor

        :param outdir:
        :param hierarchy_dir: Output directory from cellmaps_generate_hierarchy
        :type hierarchy_dir: str
        :param min_comp_size:
        :type min_comp_size: int
        :param max_fdr:
        :type max_fdr: float
        :param min_jaccard_index:
        :type min_jaccard_index: float
        :param corum:
        :type corum: str
        :param go_cc:
        :type go_cc: str
        :param hpa:
        :type hpa: str
        :param ndex_server:
        :type ndex_server: str
        :param name:
        :type name: str
        :param organization_name:
        :type organization_name: str
        :param project_name:
        :type project_name: str
        :param input_data_dict: Command line parameters
        :type input_data_dict: dict
        :param skip_logging: If ``True`` skip logging, if ``None`` or ``False`` do NOT skip logging
        :type skip_logging: bool
        :param provenance_utils: ProvenanceUtil object to use for
                                 FAIRSCAPE registration
        :type provenance_utils: py:class:`cellmaps_utils.provenance.ProvenanceUtil`
        """
        logger.debug('In constructor')
        if outdir is None:
            raise CellmapshierarchyevalError('outdir is None')
        self._outdir = os.path.abspath(outdir)
        self._hierarchy_dir = hierarchy_dir
        self._min_comp_size = min_comp_size
        self._max_fdr = max_fdr
        self._min_jaccard_index = min_jaccard_index
        self._corum = corum
        self._go_cc = go_cc
        self._hpa = hpa
        self._ndex_server = ndex_server
        self._start_time = int(time.time())
        self._geneset_agents = geneset_agents
        self._name = name
        self._project_name = project_name
        self._organization_name = organization_name
        self._input_data_dict = input_data_dict
        if skip_logging is None:
            self._skip_logging = False
        else:
            self._skip_logging = skip_logging

        self._skip_term_enrichment = skip_term_enrichment

        self._provenance_utils = provenance_utils
        self._geneset_annotator = geneset_annotator
        self._hierarchy_helper = None
        self._hierarchy_real_ids = []
        self._provenance = provenance

        if self._input_data_dict is None:
            self._input_data_dict = {'outdir': self._outdir,
                                     'hierarchy_dir': self._hierarchy_dir,
                                     'min_comp_size': self._min_comp_size,
                                     'max_fdr': self._max_fdr,
                                     'min_jaccard_index': self._min_jaccard_index,
                                     'corum': self._corum,
                                     'go_cc': self._go_cc,
                                     'hpa': self._hpa,
                                     'ndex_server': self._ndex_server,
                                     'geneset_agents': str(self._geneset_agents),
                                     'name': self._name,
                                     'project_name': self._project_name,
                                     'organization_name': self._organization_name,
                                     'skip_logging': self._skip_logging,
                                     'provenance': str(self._provenance)
                                     }

    def _term_enrichment_hierarchy(self, hierarchy):
        """
        Performs term enrichment on the given hierarchy.

        :param hierarchy: The hierarchy for which term enrichment is performed.
        :type hierarchy: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork` or :py:class:`~ndex2.cx2.CX2Network`
        :return: Updated hierarchy after term enrichment.
        :rtype: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork` or :py:class:`~ndex2.cx2.CX2Network`
        """
        # get genes in hierarchy
        hierarchy_genes = self._get_hierarchy_genes(hierarchy)

        term_definitions = [
            ('CORUM', CORUM_EnrichmentTerms, self._corum),
            ('GO_CC', GO_EnrichmentTerms, self._go_cc),
            ('HPA', HPA_EnrichmentTerms, self._hpa),
        ]

        for term_name, term_class, term_uuid in term_definitions:
            self._process_term(term_name, term_class, hierarchy, hierarchy_genes, term_uuid)

        return hierarchy

    def _get_network_from_server(self, uuid=None, max_retries=3, retry_wait=10):
        """
        :param uuid:
        :type uuid: str
        :param max_retries:
        :type max_retries: int
        :param retry_wait:
        :type retry_wait: int
        :return:
        """
        retry_num = 1
        while retry_num < max_retries:
            logger.debug('Getting term network try # ' + str(retry_num))
            try:
                return ndex2.create_nice_cx_from_server(self._ndex_server, uuid=uuid)
            except RequestException as e:
                logger.debug(f"RequestException: {str(e)}")
            except JSONDecodeError as e:
                logger.debug(f"Timeout error: {str(e)}")
            except Exception as e:
                logger.debug(f"Unexpected error: {str(e)}")
            finally:
                retry_num += 1
                time.sleep(retry_wait)
        raise CellmapshierarchyevalError(str(max_retries) + ' attempts to get network ' +
                                         str(uuid) + ' failed')

    def _process_term(self, term_name, term_class, hierarchy, hierarchy_genes, uuid):
        """
        Processes a given term by retrieving it from the server, performing enrichment testing,
        and adding the results to the hierarchy.

        :param term_name: The name of the term to be processed.
        :type term_name: str
        :param term_class: The class of the term.
        :type term_class: class
        :param hierarchy: The hierarchy in CX format.
        :type hierarchy: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork`
        :param hierarchy_genes: List of genes in the hierarchy.
        :type hierarchy_genes: list
        :param uuid: The UUID of the term.
        :type uuid: str
        """
        terms_cx = self._get_network_from_server(uuid)
        terms = term_class(terms_cx, term_name, hierarchy_genes, self._min_comp_size)
        if len(terms.term_genes) == 0:
            warnings.warn(f"Skipping {term_name} enrichment due to no genes present when "
                          f"min_comp_size set to {self._min_comp_size}")
            self._add_empty_attr_to_hierarchy(hierarchy, terms)
        else:
            enrichment_results = self._enrichment_test(hierarchy, terms, hierarchy_genes)
            self._add_results_to_hierarchy(hierarchy, terms, enrichment_results)

    def _enrichment_test(self, hierarchy, terms, hierarchy_genes):
        """
        Performs the enrichment test on the provided hierarchy and terms.

        :param hierarchy: The hierarchy in CX format.
        :type hierarchy: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork`
        :param terms: The terms for enrichment test.
        :type terms:
        :param hierarchy_genes: List of genes in the hierarchy.
        :type hierarchy_genes: list
        :return: Matrix of enrichment results.
        :rtype: np.ndarray
        """
        hierarchy_size = len(hierarchy.get_nodes())
        self._hierarchy_real_ids = self._hierarchy_helper.get_hierarchy_real_ids(hierarchy, hierarchy_size)
        term_genes_dict = terms.term_genes
        term_size = len(term_genes_dict)
        term_names = list(term_genes_dict.keys())
        enrichment_results = np.empty((hierarchy_size, term_size), dtype=object)

        # get overlap genes
        all_overlap_genes = list(set(hierarchy_genes).intersection(terms.all_term_genes))
        cap_m = len(all_overlap_genes)

        for hierarchy_index in np.arange(hierarchy_size):

            node = hierarchy.get_node(self._hierarchy_real_ids[hierarchy_index])
            node_genes = self._hierarchy_helper.get_node_genes(hierarchy, node)

            # intersection with genes in the term
            node_genes = set(node_genes).intersection(all_overlap_genes)
            n = len(node_genes)

            for term_index in np.arange(term_size):
                term = term_names[term_index]
                value = term_genes_dict[term]
                term_genes = set(value)
                cap_n = len(term_genes)
                overlap_genes = list(node_genes.intersection(term_genes))
                x = len(overlap_genes)
                pval = hypergeom.sf(x - 1, cap_m, n, cap_n)
                jaccard_index = len(overlap_genes) / len(node_genes.union(term_genes))
                result = EnrichmentResult(term, pval, jaccard_index, overlap_genes)

                if terms.term_description is not None:
                    result.set_description(terms.term_description[term])

                enrichment_results[hierarchy_index, term_index] = result

        try:
            pvals = np.array([[obj.pval for obj in row] for row in enrichment_results])
            fdr = multipletests(pvals.flatten(), method='fdr_bh')[1].reshape(pvals.shape)
        except ZeroDivisionError:
            raise CellmapshierarchyevalError(f"No genes were found with min_comp_size set to {self._min_comp_size}")

        # filter results by fdr and JI, sort by max JI
        for hierarchy_index in np.arange(enrichment_results.shape[0]):

            # set adjusted p-values and if the enrichment is accepted
            for term_index in np.arange(enrichment_results.shape[1]):
                enrichment_results[hierarchy_index, term_index].set_adjusted_pval(fdr[hierarchy_index, term_index])
                enrichment_results[hierarchy_index, term_index].set_accepted(self._min_jaccard_index, self._max_fdr)

        return enrichment_results

    def _add_results_to_hierarchy(self, hierarchy, terms, enrichment_results):
        """
        Incorporates the enrichment results into the hierarchy by adding relevant node attributes.

        :param hierarchy: The hierarchy in CX format.
        :type hierarchy: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork` or :py:class:`~ndex2.cx2.CX2Network`
        :param terms: The terms used for enrichment.
        :type terms:
        :param enrichment_results: Matrix of enrichment results.
        :type enrichment_results: np.ndarray
        """
        updated_node_ids = set()
        for hierarchy_index in np.arange(enrichment_results.shape[0]):
            node_id = self._hierarchy_real_ids[hierarchy_index]
            sorted_results = sorted(enrichment_results[hierarchy_index], key=lambda obj: obj.jaccard_index,
                                    reverse=True)
            sorted_results_threshold = [x for x in sorted_results if x.accepted]
            hierarchy.set_node_attribute(node_id, '{}_terms'.format(terms.term_name),
                                         '|'.join([x.term for x in sorted_results_threshold]))
            if terms.term_description is not None:
                hierarchy.set_node_attribute(node_id, '{}_descriptions'.format(terms.term_name),
                                             '|'.join([x.description for x in sorted_results_threshold]))
            hierarchy.set_node_attribute(node_id, '{}_FDRs'.format(terms.term_name),
                                         '|'.join(
                                             ['{:0.2e}'.format(x.adjusted_pval) for x in sorted_results_threshold]))
            hierarchy.set_node_attribute(node_id, '{}_jaccard_indexes'.format(terms.term_name),
                                         '|'.join(
                                             [str(np.round(x.jaccard_index, 2)) for x in sorted_results_threshold]))
            hierarchy.set_node_attribute(node_id, '{}_overlap_genes'.format(terms.term_name),
                                         '|'.join([','.join(x.overlap_genes) for x in sorted_results_threshold]))
            if len(sorted_results_threshold) > 0:
                hierarchy.set_node_attribute(node_id, '{}_max_jaccard_index'.format(terms.term_name),
                                             np.round(sorted_results_threshold[0].jaccard_index, 2))

            updated_node_ids.add(node_id)

        node_ids = list(set(self._hierarchy_helper.get_nodes(hierarchy)).difference(updated_node_ids))
        self._add_empty_attr_to_hierarchy(hierarchy, terms, node_ids=node_ids)

    def _add_empty_attr_to_hierarchy(self, hierarchy, terms, node_ids=None):
        """
        Adds empty attributes to nodes in the hierarchy for the given term.
        This is used when no genes are present for enrichment of the specific term.

        :param hierarchy: The hierarchy in CX format.
        :type hierarchy: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork`
        :param terms: The terms for which empty attributes should be added.
        :type terms:
        """
        if node_ids is None:
            node_ids = self._hierarchy_helper.get_nodes(hierarchy)

        for node_id in node_ids:
            hierarchy.set_node_attribute(node_id, '{}_terms'.format(terms.term_name), "")
            hierarchy.set_node_attribute(node_id, '{}_descriptions'.format(terms.term_name), "")
            hierarchy.set_node_attribute(node_id, '{}_FDRs'.format(terms.term_name), "")
            hierarchy.set_node_attribute(node_id, '{}_jaccard_indexes'.format(terms.term_name), "")
            hierarchy.set_node_attribute(node_id, '{}_overlap_genes'.format(terms.term_name), "")

    def _get_hierarchy_genes(self, hierarchy):
        """
        Extracts and returns all genes from the provided hierarchy.

        :param hierarchy: The hierarchy from which genes are extracted.
        :type hierarchy: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork` or :py:class:`~ndex2.cx2.CX2Network`
        :return: List of genes in the hierarchy.
        :rtype: list
        """
        hierarchy_genes = set()
        for _, node in self._hierarchy_helper.get_nodes(hierarchy).items():
            node_genes = self._hierarchy_helper.get_node_genes(hierarchy, node)
            hierarchy_genes.update(node_genes)
        return list(hierarchy_genes)

    def _create_rocrate(self):
        """
        Creates rocrate for output directory

        :raises CellMapsProvenanceError: If there is an error
        """
        logger.debug('Registering rocrate with FAIRSCAPE')
        name, proj_name, org_name = 'Hierarchy Enrichment', 'NA', 'NA'
        if os.path.exists(os.path.join(self._hierarchy_dir, constants.RO_CRATE_METADATA_FILE)):
            name, proj_name, org_name = self._provenance_utils.get_name_project_org_of_rocrate(self._hierarchy_dir)
        elif self._provenance is not None:
            if 'name' in self._provenance:
                self._name = self._provenance['name']
            if 'organization-name' in self._provenance:
                self._organization_name = self._provenance['organization-name']
            if 'project-name' in self._provenance:
                self._project_name = self._provenance['project-name']

        if self._name is not None:
            name = self._name

        if self._organization_name is not None:
            org_name = self._organization_name

        if self._project_name is not None:
            proj_name = self._project_name
        try:
            self._provenance_utils.register_rocrate(self._outdir,
                                                    name=name,
                                                    organization_name=org_name,
                                                    project_name=proj_name)
        except TypeError as te:
            raise CellmapshierarchyevalError('Invalid provenance: ' + str(te))
        except KeyError as ke:
            raise CellmapshierarchyevalError('Key missing in provenance: ' + str(ke))

    def _register_software(self):
        """
        Registers this tool

        :raises CellMapsImageEmbeddingError: If fairscape call fails
        """
        self._softwareid = self._provenance_utils.register_software(self._outdir,
                                                                    name=cellmaps_hierarchyeval.__name__,
                                                                    description=cellmaps_hierarchyeval.__description__,
                                                                    author=cellmaps_hierarchyeval.__author__,
                                                                    version=cellmaps_hierarchyeval.__version__,
                                                                    file_format='py',
                                                                    url=cellmaps_hierarchyeval.__repo_url__)

    def _register_computation(self, generated_dataset_ids=[]):
        """
        Registers the computation run details with the FAIRSCAPE platform.

        :param generated_dataset_ids: List of IDs for the datasets generated during the computation.
        :type generated_dataset_ids: list
        # Todo: added in used dataset, software and what is being generated
        :return:
        """
        logger.debug('Getting id of input rocrate')
        input_dataset_id = self._provenance_utils.get_id_of_rocrate(self._hierarchy_dir)
        self._provenance_utils.register_computation(self._outdir,
                                                    name=cellmaps_hierarchyeval.__name__ + ' computation',
                                                    run_by=str(self._provenance_utils.get_login()),
                                                    command=str(self._input_data_dict),
                                                    description='run of ' + cellmaps_hierarchyeval.__name__,
                                                    used_software=[self._softwareid],
                                                    used_dataset=[input_dataset_id],
                                                    generated=generated_dataset_ids)

    def _write_and_register_annotated_hierarchy_as_nodelist(self, hierarchy):
        """
        Writes out **hierarchy** passed in as node list file

        :param hierarchy:
        :type hierarchy: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork` or :py:class:`~ndex2.cx2.CX2Network`
        :return: (dataset id, path to output file)
        :rtype: tuple
        """
        logger.debug('Writing hierarchy nodelist')
        dest_path = self.get_annotated_hierarchy_as_nodelist_dest_file()

        # write node list to filesystem
        self._hierarchy_helper.write_as_nodelist(hierarchy, dest_path)

        # register node list file with fairscape
        data_dict = {'name': os.path.basename(dest_path) + ' PPI nodelist file',
                     'description': 'Annotated Nodelist file',
                     'data-format': 'tsv',
                     'author': cellmaps_hierarchyeval.__name__,
                     'version': cellmaps_hierarchyeval.__version__,
                     'date-published': date.today().strftime('%m-%d-%Y')}
        dataset_id = self._provenance_utils.register_dataset(self._outdir,
                                                             source_file=dest_path,
                                                             data_dict=data_dict)
        return dataset_id

    def _write_and_register_annotated_hierarchy(self, hierarchy):
        """
         Writes out the provided hierarchy in the CX format and registers it.

        :param hierarchy: The hierarchy to be written out.
        :type hierarchy: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork` or :py:class:`~ndex2.cx2.CX2Network`
        :return: Dataset ID for the registered hierarchy.
        :rtype: str
        """
        logger.debug('Writing hierarchy')
        hierarchy_out_file = self.get_annotated_hierarchy_dest_file()
        self._hierarchy_helper.dump_to_file(hierarchy, hierarchy_out_file)
        # register ppi network file with fairscape
        data_dict = {'name': os.path.basename(hierarchy_out_file) + ' Hierarchy network file',
                     'description': 'Hierarchy network file',
                     'data-format': self._hierarchy_helper.get_format(),
                     'author': cellmaps_hierarchyeval.__name__,
                     'version': cellmaps_hierarchyeval.__version__,
                     'date-published': date.today().strftime('%m-%d-%Y')}
        dataset_id = self._provenance_utils.register_dataset(self._outdir,
                                                             source_file=hierarchy_out_file,
                                                             data_dict=data_dict)
        return dataset_id

    def get_hierarchy_parent_network_dest_file(self):
        """
        Creates file path prefix for hierarchy parent network

        Example path: ``/tmp/foo/hierarchy_parent``
        :return:
        """
        return os.path.join(self._outdir, 'hierarchy_parent')

    def _write_and_register_hierarchy_parent_network(self, parent_input_file=None):
        """
        :param network:
        :return:
        """
        logger.debug('Writing hierarchy parent')
        suffix = self._hierarchy_helper.get_suffix()
        parent_out_file = self.get_hierarchy_parent_network_dest_file() + suffix
        shutil.copy(parent_input_file, parent_out_file)

        data_dict = {'name': 'Hierarchy parent network',
                     'description': 'Hierarchy parent network file',
                     'keywords': ['file', 'parent', 'interactome', 'ppi', 'network', 'CX2'],
                     'data-format': self._hierarchy_helper.get_format(),
                     'author': cellmaps_hierarchyeval.__name__,
                     'version': cellmaps_hierarchyeval.__version__,
                     'date-published': date.today().strftime(self._provenance_utils.get_default_date_format_str())}
        dataset_id = self._provenance_utils.register_dataset(self._outdir,
                                                             source_file=parent_out_file,
                                                             data_dict=data_dict)
        return dataset_id

    def get_annotated_hierarchy_dest_file(self):
        """
        Creates file path prefix for hierarchy

        Example path: ``/tmp/foo/hierarchy``

        :return: Prefix path on filesystem to write Hierarchy Network
        :rtype: str
        """
        return os.path.join(self._outdir, constants.HIERARCHY_NETWORK_PREFIX +
                            self._hierarchy_helper.get_suffix())

    def get_annotated_hierarchy_as_nodelist_dest_file(self):
        """
        Creates file path prefix for hierarchy

        Example path: ``/tmp/foo/hierarchy``

        :return: Prefix path on filesystem to write Hierarchy Network
        :rtype: str
        """
        return os.path.join(self._outdir, constants.HIERARCHY_NODES_FILE)

    def _update_annotate_hierarchy(self, network=None, path=None):
        """
        Adds HCX attributes to network as well as sets

        ``prov:wasGeneratedBy`` to the name and version of this tool

        ``prov:wasDerivedFrom`` to FAIRSCAPE dataset id of this rocrate

        :param network: Hierarchy
        :type network: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork` or :py:class:`~ndex2.cx2.CX2Network`
        :param path: Path to largest PPI network in CX or CX2 format
        :type path: str
        """
        network.add_network_attribute('prov:wasGeneratedBy',
                                      cellmaps_hierarchyeval.__name__ + ' ' + cellmaps_hierarchyeval.__version__)

        rocrate_id = self._provenance_utils.get_id_of_rocrate(path)
        network.add_network_attribute('prov:wasDerivedFrom', 'RO-crate: ' + str(rocrate_id))

    def initialize_hierarchy_helper(self):
        """
        Initializes hierarchy helper which will  be used to call custom methods
        depending on whether the input was in CX or CX2 format.

        :return:
        """
        cx_file_path = os.path.join(self._hierarchy_dir, constants.HIERARCHY_NETWORK_PREFIX + constants.CX_SUFFIX)
        if os.path.exists(cx_file_path):
            self._hierarchy_helper = NiceCXNetworkHelper(cx_file_path)
            return

        cx2_file_path = os.path.join(self._hierarchy_dir, constants.HIERARCHY_NETWORK_PREFIX + constants.CX2_SUFFIX)
        if os.path.exists(cx2_file_path):
            self._hierarchy_helper = CX2NetworkHelper(cx2_file_path)
            return

        raise CellmapshierarchyevalError(f"Input directory '{self._hierarchy_dir}' does not contain neither cx nor cx2 "
                                         f"files.")

    def _annotate_hierarchy_with_geneset_annotators(self, hierarchy=None):
        """
        Annotates hierarchy with each GeneSetAgent set in constructor.
        """
        if self._geneset_annotator is None:
            logger.debug('Skipping because geneset_annotator is None')
            return
        if self._geneset_agents is None:
            logger.debug('Skipping because there are no geneset agents')
            return
        self._geneset_annotator.set_hierarchy_helper(self._hierarchy_helper)
        logger.debug('Processing ' + str(len(self._geneset_agents)) + ' geneset agents')
        for a in tqdm(self._geneset_agents, desc='GeneSet Agents'):
            self._geneset_annotator.annotate_hierarchy(hierarchy=hierarchy,
                                                       geneset_agent=a)

    def generate_readme(self):
        description = getattr(cellmaps_hierarchyeval, '__description__', 'No description provided.')
        version = getattr(cellmaps_hierarchyeval, '__version__', '0.0.0')

        with open(os.path.join(os.path.dirname(__file__), 'readme_outputs.txt'), 'r') as f:
            readme_outputs = f.read()

        readme = readme_outputs.format(DESCRIPTION=description, VERSION=version)
        with open(os.path.join(self._outdir, 'README.txt'), 'w') as f:
            f.write(readme)

    def run(self):
        """
        Evaluates CM4AI Hierarchy

        :return:
        """
        exitcode = 99
        try:
            logger.debug('In run method')

            if os.path.isdir(self._outdir):
                raise CellmapshierarchyevalError(self._outdir + ' already exists')

            if not os.path.isdir(self._outdir):
                os.makedirs(self._outdir, mode=0o755)

            if self._skip_logging is False:
                logutils.setup_filelogger(outdir=self._outdir,
                                          handlerprefix='cellmaps_hierarchyeval')
            logutils.write_task_start_json(outdir=self._outdir,
                                           start_time=self._start_time,
                                           data={'commandlineargs': self._input_data_dict},
                                           version=cellmaps_hierarchyeval.__version__)

            self.generate_readme()

            self._create_rocrate()
            self._register_software()

            # Initializes hierarchy helper
            self.initialize_hierarchy_helper()

            generated_dataset_ids = []

            # annotate hierarchy
            hierarchy = self._hierarchy_helper.get_hierarchy()
            if self._skip_term_enrichment is None or self._skip_term_enrichment is False:
                hierarchy = self._term_enrichment_hierarchy(hierarchy)
            else:
                logger.info('Skipping term enrichment because '
                            'skip_term_enrichment flag is True')

            self._annotate_hierarchy_with_geneset_annotators(hierarchy=hierarchy)

            self._update_annotate_hierarchy(hierarchy, self._outdir)

            # write out annotated hierarchy
            dataset_id = self._write_and_register_annotated_hierarchy(hierarchy)
            generated_dataset_ids.append(dataset_id)

            # write out parent network
            parent_network_path = os.path.join(
                self._hierarchy_dir, 'hierarchy_parent' + self._hierarchy_helper.get_suffix())
            if os.path.exists(parent_network_path):
                dataset_id = self._write_and_register_hierarchy_parent_network(parent_network_path)
                generated_dataset_ids.append(dataset_id)

            # write out nodes file
            dataset_id = self._write_and_register_annotated_hierarchy_as_nodelist(hierarchy)
            generated_dataset_ids.append(dataset_id)

            # register generated datasets
            self._register_computation(generated_dataset_ids=generated_dataset_ids)
            exitcode = 0
        finally:
            logutils.write_task_finish_json(outdir=self._outdir,
                                            start_time=self._start_time,
                                            status=exitcode)

        return exitcode
