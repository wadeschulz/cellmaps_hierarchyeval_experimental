#! /usr/bin/env python
import os
import logging
import time
import json
import numpy as np
from datetime import date
from scipy.stats import hypergeom 
from statsmodels.stats.multitest import multipletests
import ndex2
from cellmaps_utils import constants
from cellmaps_utils import logutils
from cellmaps_utils.provenance import ProvenanceUtil
import cellmaps_hierarchyeval
from cellmaps_hierarchyeval.exceptions import CellmapshierarchyevalError



logger = logging.getLogger(__name__)

NDEX_UUID = {
    'HPA': 'a6a88e2d-9c0f-11ed-9a1f-005056ae23aa',
    'CORUM' : '764f7471-9b79-11ed-9a1f-005056ae23aa',
    'GO_CC' : 'f484e8ee-0b0f-11ee-aa50-005056ae23aa'}


class EnrichmentTerms(object):
    """
    Base class for implementations that generate
    term databases for enrichment (i.e., HPA, CORUM, GO)
    """
    def __init__(self, terms=None, term_name=None, hierarchy_genes=None, min_comp_size=4
                 ):
        """
        Constructor
        """
        self.terms = terms
        self.term_name = term_name
        self.hierarchy_genes = hierarchy_genes
        self.min_comp_size = min_comp_size
        self.terms = None
        self.term_genes = None
        self.term_description = None

class GO_EnrichmentTerms(EnrichmentTerms):
    def __init__(self, terms=None, term_name=None, hierarchy_genes=None, min_comp_size=4
                 ):
        super().__init__(terms=terms, term_name = term_name, hierarchy_genes = hierarchy_genes, min_comp_size=min_comp_size)
        self.term_genes = self._get_term_genes(terms)
        self.term_description = self._get_term_description(terms)

    def _get_term_genes(self, terms):
        term_genes_dict = {}
        for id, node in terms.get_nodes():
            term = node.get('n')
            genes =  terms.get_node_attribute_value(node, 'genes')
            if genes is None:
                continue
            genes = genes.split(',')
            genes = list(set(genes).intersection(set(self.hierarchy_genes)))
            if len(genes) < self.min_comp_size:
                continue
            term_genes_dict[term] = genes
        return term_genes_dict
            
    def _get_term_description(self, terms):
        term_description = {}
        for id, node in terms.get_nodes():
            term = node.get('n')
            term_description[term] = terms.get_node_attribute_value(node, 'description')
        return term_description    

    
class CORUM_EnrichmentTerms(EnrichmentTerms):
    def __init__(self, terms=None, term_name=None, hierarchy_genes=None, min_comp_size=4
                 ):
        super().__init__(terms=terms, term_name = term_name, hierarchy_genes = hierarchy_genes, min_comp_size=min_comp_size)
        self.term_genes = self._get_term_genes(terms)
        self.term_description =  None
        
    def _get_term_genes(self, terms):
        term_genes_dict = {}
        for id, node in terms.get_nodes():
            genes = terms.get_node_attribute_value(node, 'subunits(Gene name)')
            if genes is None:
                continue
            genes = list(set(genes).intersection(set(self.hierarchy_genes)))
            if len(genes) < self.min_comp_size:
                continue
            term = node.get('n')
            term_genes_dict[term] = genes
        return term_genes_dict

class HPA_EnrichmentTerms(EnrichmentTerms):
    def __init__(self, terms=None, term_name=None, hierarchy_genes=None, min_comp_size=4
                 ):
        super().__init__(terms=terms, term_name = term_name, hierarchy_genes = hierarchy_genes, min_comp_size=min_comp_size)
        self.term_genes = self._get_term_genes(terms)
        self.term_description =  None
        
    def _get_term_genes(self, terms):
        term_genes_dict = {}
        for id, node in terms.get_nodes():
            node_name = node.get('n')
            if node_name not in self.hierarchy_genes:
                continue
            for a in ['Main location', 'Additional location']:
                annotations = terms.get_node_attribute_value(node, a)
                if annotations is None:
                    continue
                for c in annotations:
                    if c in term_genes_dict:
                        term_genes_dict[c].append(node_name)
                    else:
                        term_genes_dict[c] = [node_name]
        return term_genes_dict
 

    
class EncirhmentResult(object):
    """
    Base class for generating hierarchy
    that is output in CX format following
    CDAPS style
    """
    def __init__(self,
                 term=None,
                 pval=None,
                 jaccard_index = None,
                 overlap_genes = None):
        """
        Constructor
        """
        self.term = term
        self.description = None
        self.pval = pval
        self.jaccard_index = jaccard_index
        self.overlap_genes = overlap_genes
        self.adjusted_pval = None
        self.accepted = False
        
    def set_adjusted_pval(self, adjusted_pval):
        self.adjusted_pval = adjusted_pval
        
    def set_accepted(self, min_jaccard_index, max_fdr):
        if self.jaccard_index == 1:
            self.accepted = True
        elif (self.jaccard_index >= min_jaccard_index) & (self.adjusted_pval < max_fdr):
            self.accepted = True

    def set_description(self, description):
        self.description = description
        
class CellmapshierarchyevalRunner(object):
    """
    Class to run algorithm
    """
    def __init__(self, outdir=None, 
                 hierarchy_dir=None,
                 min_comp_size = 4,
                 max_fdr = 0.05,
                 min_jaccard_index = 0.1,
                 name=None,
                 organization_name=None,
                 project_name=None,
                 input_data_dict=None,
                 provenance_utils=ProvenanceUtil()):
        """
        Constructor

        :param outdir: Directory to create and put results in
        :type outdir: str
        :param hierarchy_dir: Output directory from cellmaps_generate_hierarchy 
        :param name:
        :param organization_name:
        :param project_name:
        :param provenance_utils:
        """
        logger.debug('In constructor')
        if outdir is None:
            raise CellmapshierarchyevalError('outdir is None')
        self._outdir = os.path.abspath(outdir)
        self._hierarchy_dir = hierarchy_dir
        self._min_comp_size = min_comp_size
        self._max_fdr = max_fdr
        self._min_jaccard_index = min_jaccard_index
        self._start_time = int(time.time())
        self._name = name
        self._project_name = project_name
        self._organization_name = organization_name
        self._input_data_dict = input_data_dict
        self._provenance_utils = provenance_utils

        
    def _get_hierarchy_file(self):
        """

        :return:
        """
        return os.path.join(self._hierarchy_dir,
                            constants.HIERARCHY_NETWORK_PREFIX)
    
    def _term_enrichment_hierarchy(self, hierarchy):
            
        hierarchy_cx = ndex2.create_nice_cx_from_file(hierarchy)
        
        #get genes in hierarchy
        hierarchy_genes = self._get_hierarchy_genes(hierarchy_cx)
        M = len(hierarchy_genes)
        
        term_name = 'CORUM'
        uuid = NDEX_UUID[term_name]
        terms_cx =  ndex2.create_nice_cx_from_server('http://www.ndexbio.org',uuid=uuid)
        terms = CORUM_EnrichmentTerms(terms_cx, term_name, hierarchy_genes, self._min_comp_size)
        enrichment_results = self._enrichment_test(hierarchy_cx, terms, M)
        self._add_results_to_hierarchy(hierarchy_cx, terms, enrichment_results)
        
        term_name = 'GO_CC'
        uuid = NDEX_UUID[term_name]
        terms_cx =  ndex2.create_nice_cx_from_server('http://www.ndexbio.org',uuid=uuid)
        terms = GO_EnrichmentTerms(terms_cx, term_name, hierarchy_genes, self._min_comp_size)
        enrichment_results = self._enrichment_test(hierarchy_cx, terms, M)
        self._add_results_to_hierarchy(hierarchy_cx, terms, enrichment_results)
        
        term_name = 'HPA'
        uuid = NDEX_UUID[term_name]
        terms_cx =  ndex2.create_nice_cx_from_server('http://www.ndexbio.org',uuid=uuid)
        terms = HPA_EnrichmentTerms(terms_cx, term_name, hierarchy_genes, self._min_comp_size)
        enrichment_results = self._enrichment_test(hierarchy_cx, terms, M)
        self._add_results_to_hierarchy(hierarchy_cx, terms, enrichment_results)
        
        return hierarchy_cx
            
    def _enrichment_test(self, hierarchy, terms, M):
        
        hierarchy_size = len(hierarchy.nodes)
        
        term_genes_dict = terms.term_genes
        term_size = len(term_genes_dict)
        term_names = list(term_genes_dict.keys())
        enrichment_results = np.empty((hierarchy_size, term_size), dtype=object)
        
        for hierarchy_index in np.arange(hierarchy_size):

            node = hierarchy.get_node(hierarchy_index)
            node_name = node['n']

            node_genes = set(hierarchy.get_node_attribute(node, 'CD_MemberList')['v'].split(' '))
            n = len(node_genes)

            for term_index in np.arange(term_size):
                term = term_names[term_index]
                value = term_genes_dict[term]

                term_genes = set(value)
                N = len(term_genes)
                overlap_genes = list(node_genes.intersection(term_genes))
                x = len(overlap_genes)
                pval = hypergeom.sf(x - 1, M, n, N)
                jaccard_index = len(overlap_genes) / len(node_genes.union(term_genes))
                result = EncirhmentResult(term, pval, jaccard_index, overlap_genes)

                if terms.term_description is not None:
                    result.set_description(terms.term_description[term])

                enrichment_results[hierarchy_index, term_index] = result

        pvals = np.array([[obj.pval for obj in row] for row in enrichment_results])
        fdr = multipletests(pvals.flatten(), method='fdr_bh')[1].reshape(pvals.shape)
        
        # filter results by fdr and JI, sort by max JI
        for hierarchy_index in np.arange(enrichment_results.shape[0]):

            #set adjusted p-values and if the enrichment is accepted
            for term_index in np.arange(enrichment_results.shape[1]):
                enrichment_results[hierarchy_index, term_index].set_adjusted_pval(fdr[hierarchy_index,term_index])
                enrichment_results[hierarchy_index, term_index].set_accepted(self._min_jaccard_index, self._max_fdr)
        
        return enrichment_results

    def _add_results_to_hierarchy(self, hierarchy, terms, enrichment_results):
        
        for hierarchy_index in np.arange(enrichment_results.shape[0]):
            node = hierarchy.get_node(hierarchy_index)
            sorted_results  = sorted(enrichment_results[hierarchy_index], key=lambda obj: obj.jaccard_index, reverse=True)
            sorted_results_threshold = [x for x in sorted_results if x.accepted]
            hierarchy.set_node_attribute(node, '{}_terms'.format(terms.term_name), '|'.join([x.term for x in sorted_results_threshold]))
            if terms.term_description is not None:
                hierarchy.set_node_attribute(node, '{}_descriptions'.format(terms.term_name), '|'.join([x.description for x in sorted_results_threshold]))
            hierarchy.set_node_attribute(node, '{}_FDRs'.format(terms.term_name), '|'.join([str(x.adjusted_pval) for x in sorted_results_threshold]))
            hierarchy.set_node_attribute(node, '{}_jaccard_indexes'.format(terms.term_name), '|'.join([str(x.jaccard_index) for x in sorted_results_threshold]))
            hierarchy.set_node_attribute(node, '{}_overlap_genes'.format(terms.term_name), '|'.join([','.join(x.overlap_genes) for x in sorted_results_threshold]))

    def _get_hierarchy_genes(self, hierarchy_cx):
        hierarchy_genes = set()
        for node in hierarchy_cx.nodes:
            node_genes = hierarchy_cx.get_node_attribute(node, 'CD_MemberList')['v'].split(' ')
            hierarchy_genes.update(node_genes)
        return list(hierarchy_genes)
    
    def _create_rocrate(self):
        """
        Creates rocrate for output directory

        :raises CellMapsProvenanceError: If there is an error
        """
        logger.debug('Registering rocrate with FAIRSCAPE')
        name, proj_name, org_name = self._provenance_utils.get_name_project_org_of_rocrate(self._hierarchy_dir)

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
                                                                    file_format='.py',
                                                                    url=cellmaps_hierarchyeval.__repo_url__)

    def _register_computation(self, generated_dataset_ids=[]):
        """
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

        
        
    def _write_and_register_annotated_hierarchy_as_nodelist(self, hierarchy, dest_path=None):
        """
        Writes out **hierarchy** passed in as node list file

        :param hierarchy:
        :type hierarchy: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork`
        :return: (dataset id, path to output file)
        :rtype: tuple
        """
        logger.debug('Writing hierarchy nodelist')
        dest_path = os.path.join(self._outdir, constants.HIERARCHY_NODES_PREFIX)

        # write node list to filesystem
        with open(dest_path, 'w') as f:
            
            # write headers
            f.write('Name' + '\t')
            for a in hierarchy.get_node_attributes(0):
                f.write(a['n'] + '\t')
            f.write('\n')
            
            #write node attributes
            for node_id, node_obj in hierarchy.get_nodes():
                f.write(node_obj['n'] + '\t')
                for a in hierarchy.get_node_attributes(node_obj):
                    f.write(a['v'] + '\t')
                f.write('\n')

        # register node list file with fairscape
        data_dict = {'name': os.path.basename(dest_path) + ' PPI edgelist file',
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

        :param network:
        :return:
        """
        logger.debug('Writing hierarchy')
        hierarchy_out_file = self.get_annotated_hierarchy_dest_file(hierarchy) + constants.CX_SUFFIX
        with open(hierarchy_out_file, 'w') as f:
            json.dump(hierarchy.to_cx(), f)
            # register ppi network file with fairscape
            data_dict = {'name': os.path.basename(hierarchy_out_file) + ' Hierarchy network file',
                         'description': 'Hierarchy network file',
                         'data-format': 'CX',
                         'author': cellmaps_hierarchyeval.__name__,
                         'version': cellmaps_hierarchyeval.__version__,
                         'date-published': date.today().strftime('%m-%d-%Y')}
            dataset_id = self._provenance_utils.register_dataset(self._outdir,
                                                                 source_file=hierarchy_out_file,
                                                                 data_dict=data_dict)
        return dataset_id
        

    def get_annotated_hierarchy_dest_file(self, hierarchy):
        """
        Creates file path prefix for hierarchy

        Example path: ``/tmp/foo/hierarchy``

        :param hierarchy: Hierarchy Network
        :type hierarchy: :py:class:`ndex2.nice_cx_network.NiceCXNetwork`
        :return: Prefix path on filesystem to write Hierarchy Network
        :rtype: str
        """
        return os.path.join(self._outdir, constants.HIERARCHY_NETWORK_PREFIX)
    
    def get_annotated_hierarchy_as_nodelist_dest_file(self, hierarchy):
        """
        Creates file path prefix for hierarchy

        Example path: ``/tmp/foo/hierarchy``

        :param hierarchy: Hierarchy Network
        :type hierarchy: :py:class:`ndex2.nice_cx_network.NiceCXNetwork`
        :return: Prefix path on filesystem to write Hierarchy Network
        :rtype: str
        """
        return os.path.join(self._outdir, constants.HIERARCHY_NODES_PREFIX)
    
    def get_hierarchy_input_file(self):
        """
        Creates file path prefix for hierarchy

        Example path: ``/tmp/foo/hierarchy``

        :param hierarchy: Hierarchy Network
        :type hierarchy: :py:class:`ndex2.nice_cx_network.NiceCXNetwork`
        :return: Prefix path on filesystem to write Hierarchy Network
        :rtype: str
        """
        return os.path.join(self._hierarchy_dir, constants.HIERARCHY_NETWORK_PREFIX)
    
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

            logutils.setup_filelogger(outdir=self._outdir,
                                      handlerprefix='cellmaps_image_embedding')
            logutils.write_task_start_json(outdir=self._outdir,
                                           start_time=self._start_time,
                                           data={'commandlineargs': self._input_data_dict},
                                           version=cellmaps_hierarchyeval.__version__)
            self._create_rocrate()
            self._register_software()
            
            generated_dataset_ids = []

            # annotate hierarchy
            hierarchy = self.get_hierarchy_input_file() + constants.CX_SUFFIX
            hierarchy = self._term_enrichment_hierarchy(hierarchy)

            # write out annotated hierarchy
            dataset_id = self._write_and_register_annotated_hierarchy(hierarchy)
            generated_dataset_ids.append(dataset_id)
            
            #write out nodes file
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
