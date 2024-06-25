
import os
import logging
import math
import pandas as pd
from scipy.stats import ranksums
from scipy import stats
import cellmaps_utils.music_utils as music_utils


logger = logging.getLogger(__name__)


class PerturbSeqAnalysis(object):
    """
    Contains utilities to compare Perturbation data
    against hierarchy passed in via constructor
    """
    def __init__(self, hierarchy, hierarchy_parent):
        """
        Constructor

        :param hierarchy:
        :type hierarchy: :py:class:`~ndex2.cx2.CX2Network`
        :param hierarchy_parent:
        :type hierarchy_parent: :py:class:`~ndex2.cx2.CX2Network`
        """
        self._hierarchy = hierarchy
        self._hierarchy_parent = hierarchy_parent

    def get_heatmap_for_given_hierarchy_system(self, hier_system_node_id,
                                               perturbseq_df, num_perturb_seq=25):
        """
        Given an id for a system in hierarchy **hier_system_node_id** and
        Perturb-seq data **perturbseq_df** create a heatmap of X most variable
        Perturb-seq proteins.

        This is done by filtering **perturbseq_df** for rows that match genes in
        given system and then keeping **num_perturb_seq** most variable columns

        :param hier_system_node_id: node id system to analyze
        :type hier_system_node_id: int
        :param perturbseq_df:
        :type perturbseq_df: :py:class:`pandas.DataFrame`
        :return: heat map table
        :rtype: :py:class:`pandas.DataFrame`
        """
        # TODO: load hierarchy and get genes for that system. put into variable assembly_genes
        assembly_genes = []
        cluster_genes_in_perturb = [x for x in assembly_genes if x in perturbseq_df.index.values]

        # from notebook but changed to match these variables
        variance_per_column = perturbseq_df.var()
        most_variable = variance_per_column.sort_values(ascending=False).head(num_perturb_seq).index.values
        data = perturbseq_df.loc[cluster_genes_in_perturb, most_variable]
        data = data.apply(stats.zscore, axis=1)
        return data

    def get_root_node_pair_similarities(self):
        """
        Gets root node pair similarities
        # Todo: explain this better

        :return:
        :rtype: :py:class:`pandas.DataFrame`
        """
        # get pairs only in root node
        for nodeid, node in self._hierarchy.get_nodes().items():
            if node['v']['HCX::isRoot'] == True:
                root_node = nodeid
                genes = node['v']['CD_MemberList'].split(' ')
                logger.debug('Number of genes: ' + str(len(genes)))

        communities_connected_to_root = []
        for edgeid, edge in self._hierarchy_parent.get_edges().items():
            if edge['s'] == root_node:
                communities_connected_to_root.append(edge['t'])
        # print(len(communities_connected_to_root))

        root_pairs = pd.DataFrame(1, index=genes, columns=genes)
        for community in communities_connected_to_root:
            community_node = self._hierarchy.get_node(community)
            community_genes = community_node['v']['CD_MemberList'].split(' ')
            root_pairs.loc[community_genes, community_genes] = 0
        return root_pairs

    def get_root_overlapping_pair_simliarities(self, root_pairs, perturbseq_df):
        """
        Get similarity scores from **perturbseq_df** that match genes attached to the root
        node of the hierarchy

        :param root_pairs:
        :type root_pairs: # TODO: what is this
        :param perturbseq_df:
        :type perturbseq_df: :py:class:`pandas.DataFrame`
        :return: (:py:func:`numpy.array`, None) # TODO: replace None with actual type and explanation
        :rtype: tuple
        """
        overlap_genes = list(set(root_pairs.index.values).intersection(set(perturbseq_df.index.values)))
        overlap_functional_data = perturbseq_df.loc[overlap_genes]
        functional_data_similarity = music_utils.cosine_similarity_scaled(overlap_functional_data)
        overlap_root_pairs = root_pairs.loc[overlap_genes, overlap_genes]

        return functional_data_similarity, overlap_root_pairs

    def get_data_simliarity_root(self, functional_data_similarity, overlap_root_pairs):
        """
        # TODO: explain what this is doing

        :param functional_data_similarity:
        :type functional_data_similarity: # TODO: what is this
        :param overlap_root_pairs:
        :type overlap_root_pairs: # TODO: what is this
        :return:
        :rtype: list
        """
        root_mask = overlap_root_pairs > 0
        functional_data_similarity_root = [x for x in
                                           music_utils.upper_tri_values(functional_data_similarity[root_mask]) if
                                           not math.isnan(x)]
        return functional_data_similarity_root

    def get_cluster_similarity(self, functional_data_similarity, hier_system_node_id):
        """
        # TODO: explain what this is doing

        :param functional_data_similarity:
        :type functional_data_similarity: # TODO: what is this
        :param hier_system_node_id:
        :type hier_system_node_id: int
        :return:
        :rtype: :py:func:`numpy.array`
        """
        cluster_genes = [] # TODO: set cluster_genes to list of genes in assembly with id hier_system_node_id
        cluster_genes_in_functional_data = [x for x in cluster_genes if x in functional_data_similarity.index.values]
        return music_utils.upper_tri_values(functional_data_similarity.loc[cluster_genes_in_functional_data,
                                                                           cluster_genes_in_functional_data])
