import logging
import math
import pandas as pd
from scipy import stats
import cellmaps_utils.music_utils as music_utils
from ndex2 import constants
from scipy.stats import ranksums

from cellmaps_hierarchyeval.exceptions import CellmapshierarchyevalError

logger = logging.getLogger(__name__)


class PerturbSeqAnalysis(object):
    """
    Contains utilities to compare Perturbation data
    against hierarchy passed in via constructor
    """

    def __init__(self, hierarchy, hierarchy_parent=None):
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
        :param num_perturb_seq:
        :type num_perturb_seq: int
        :return: heat map table
        :rtype: :py:class:`pandas.DataFrame`
        """
        node_values = self._hierarchy.get_node(hier_system_node_id)
        assembly_genes = node_values[constants.ASPECT_VALUES]['CD_MemberList'].split(' ')
        cluster_genes_in_perturb = [x for x in assembly_genes if x in perturbseq_df.index.values]

        # from notebook but changed to match these variables
        variance_per_column = perturbseq_df.var()
        most_variable = variance_per_column.sort_values(ascending=False).head(num_perturb_seq).index.values
        data = perturbseq_df.loc[cluster_genes_in_perturb, most_variable]
        data = data.apply(stats.zscore, axis=1)
        return data

    def get_root_gene_pair_similarities(self):
        """
        Calculates similarity scores between gene pairs in the root node of a hierarchy. Genes in the same community
        linked to the root node are marked with a similarity of 0, indicating they are directly related,
        while all other pairs are set to 1, suggesting no direct relation.

        :return: A DataFrame with genes as both rows and columns, populated with similarity scores.
        :rtype: :py:class:`pandas.DataFrame`
        """
        root_node = None
        genes = []
        # get pairs only in root node
        for nodeid, node in self._hierarchy.get_nodes().items():
            if node[constants.ASPECT_VALUES]['HCX::isRoot']:
                root_node = nodeid
                genes = node[constants.ASPECT_VALUES]['CD_MemberList'].split(' ')

        if root_node is None:
            raise CellmapshierarchyevalError('No root node detected!')

        # Get nodes that are directly connected to root node
        communities_connected_to_root = []
        for edgeid, edge in self._hierarchy.get_edges().items():
            if edge[constants.EDGE_SOURCE] == root_node:
                communities_connected_to_root.append(edge[constants.EDGE_TARGET])

        # Assign similarity scores
        root_pairs = pd.DataFrame(1, index=genes, columns=genes)
        for community in communities_connected_to_root:
            community_node = self._hierarchy.get_node(community)
            community_genes = community_node[constants.ASPECT_VALUES]['CD_MemberList'].split(' ')
            root_pairs.loc[community_genes, community_genes] = 0
        return root_pairs

    @staticmethod
    def get_root_overlapping_pair_similarities(root_pairs, perturbseq_df):
        """
        Get similarity scores from **perturbseq_df** that match genes attached to the root
        node of the hierarchy

        :param root_pairs: A DataFrame representing similarity scores between all genes in the root node,
                            where genes within the same community connected to the root have a score of 0,
                            indicating direct relation, and all other pairs have a score of 1,
                            indicating no direct relation.
        :type root_pairs: :py:class:`pandas.DataFrame`
        :param perturbseq_df:
        :type perturbseq_df: :py:class:`pandas.DataFrame`
        :return: A tuple containing:
             - A DataFrame of scaled cosine similarity scores for overlapping genes in communities direct to root
                and Perturb-seq data.
             - A DataFrame of root-associated similarity scores, filtered to only include overlapping genes.
        :rtype: tuple
        """
        overlap_genes = list(set(root_pairs.index.values).intersection(set(perturbseq_df.index.values)))
        overlap_functional_data = perturbseq_df.loc[overlap_genes]
        functional_data_similarity = music_utils.cosine_similarity_scaled(overlap_functional_data)
        overlap_root_pairs = root_pairs.loc[overlap_genes, overlap_genes]

        return functional_data_similarity, overlap_root_pairs

    @staticmethod
    def get_root_functional_data_similarity(functional_data_similarity, overlap_root_pairs):
        """
        Extracts and returns a list of functional similarity scores for gene pairs that are not in the same community,
            based on a filtered upper triangle extraction of the similarity matrix (ensures that only unique,
            non-redundant gene pair comparisons are considered).

        :param functional_data_similarity: A DataFrame of scaled cosine similarity scores for overlapping genes in
                                            communities direct to root and Perturb-seq data.
        :type functional_data_similarity: :py:class:`pandas.DataFrame`
        :param overlap_root_pairs: A DataFrame of root-associated similarity scores, filtered to only include
                                    overlapping genes. A score of 0 indicates a direct relation (same community)
                                    and scores greater than 0 indicate no direct relation
        :type overlap_root_pairs: :py:class:`pandas.DataFrame`
        :return: A list of non-NaN similarity scores for gene pairs that are not directly related.
        :rtype: list
        """
        root_mask = overlap_root_pairs > 0
        root_functional_data_similarity = [x for x in
                                           music_utils.upper_tri_values(functional_data_similarity[root_mask]) if
                                           not math.isnan(x)]
        return root_functional_data_similarity

    def get_cluster_similarity(self, functional_data_similarity, hier_system_node_id):
        """
        Retrieves the upper triangle similarity scores for genes within a specific cluster of a hierarchy.
        The scores are extracted from a DataFrame that contains scaled cosine similarity scores for genes that overlap
        between communities direct to root and Perturb-seq data.

        :param functional_data_similarity: A DataFrame of scaled cosine similarity scores for overlapping genes in
                                            communities direct to root and Perturb-seq data.
        :type functional_data_similarity: :py:class:`pandas.DataFrame`
        :param hier_system_node_id: The identifier for a specific node within a hierarchy.
        :type hier_system_node_id: int
        :return: An array of similarity scores from the upper triangle portion of the matrix for the specified cluster.
        :rtype: :py:func:`numpy.array`
        """
        # Retrieve node information and extract gene members from the specified cluster.
        node_values = self._hierarchy.get_node(hier_system_node_id)
        cluster_genes = node_values[constants.ASPECT_VALUES]['CD_MemberList'].split(' ')

        # Filter genes that are present in the functional data similarity matrix.
        cluster_genes_in_functional_data = [x for x in cluster_genes if x in functional_data_similarity.index.values]

        # Extract the relevant portion of the similarity matrix and return upper triangle values.
        cluster_functional_data_similarity = music_utils.upper_tri_values(functional_data_similarity.loc[
            cluster_genes_in_functional_data, cluster_genes_in_functional_data])

        return cluster_functional_data_similarity

    @staticmethod
    def compare_cluster_root_similarities(cluster_functional_data_similarity, root_functional_data_similarity):
        """
        Performs a rank-sum test to compare the distribution of functional data similarity scores
        between a specific cluster and gene pairs in root. This test helps determine if the similarity scores
        in the cluster are statistically significantly greater than those in the root.

        :param cluster_functional_data_similarity: An array of similarity scores within a specific cluster.
        :type cluster_functional_data_similarity: numpy.array
        :param root_functional_data_similarity: A list of non-NaN similarity scores for gene pairs not directly related
                                                in the root.
        :type root_functional_data_similarity: list
        :return: A tuple containing the test statistic and the p-value of the rank-sum test.
        :rtype: (float, float)
        """
        statistic, p_value = ranksums(cluster_functional_data_similarity, root_functional_data_similarity,
                                      alternative='greater')

        return statistic, p_value
