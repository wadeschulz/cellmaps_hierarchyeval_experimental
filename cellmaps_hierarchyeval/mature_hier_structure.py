
import math

import os
import argparse
import networkx as nx
#function for make ontologies
from cellmaps_hierarchyeval.make_hidefOnt_utils import *

def to_pandas_dataframe(G):
    e = G.edges(data=True)
    df = pd.DataFrame()
    df['source'] = [x[0] for x in e]
    df['target'] = [x[1] for x in e]
    df['type'] = [x[2]['type'] for x in e]
    return df

def get_termStats(G, hiergeneset):
    clusters = list(set(list(G.nodes())) - hiergeneset)
    tsize_list = []
    cgene_list = []
    descendent_list = []
    for c in clusters:
        infoset = nx.descendants(G, c)
        cgeneset = infoset.intersection(hiergeneset)
        tsize_list.append(len(cgeneset))
        cgene_list.append(list(cgeneset))
        descendent_list.append(list(infoset - cgeneset))
    df = pd.DataFrame(index=clusters)
    df['tsize'] = tsize_list
    df['genes'] = cgene_list
    df['descendent'] = descendent_list
    return df

def jaccard(A, B):
    if type(A) != set:
        A = set(A)
    if type(B) != set:
        B = set(B)
    return len(A.intersection(B)) / len(A.union(B))

def clean_shortcut(G):
    edge_df = to_pandas_dataframe(G)
    edge_df.columns = ['parent', 'child', 'type']
    for idx, row in edge_df.iterrows():
        if len(list(nx.all_simple_paths(G, row['parent'], row['child']))) > 1:
            G.remove_edge(row['parent'], row['child'])
            print('shortcut edges is removed between {} and {}'.format(row['parent'], row['child']))
    return

def reorganize(G, hiergeneset, ci_thre): # Add an edge if the nodes have containment index >=threshold
    iterate = True
    n_iter = 1
    while iterate:
        clear = True
        print('... starting iteration {}'.format(n_iter))
        ts_df = get_termStats(G, hiergeneset) # get the termStats from the networkx 
        ts_df.sort_values('tsize', ascending=False, inplace=True) 
        for comp, row in ts_df.iterrows():
            tmp = ts_df[ts_df['tsize'] < row['tsize']] # get all components smaller than this components 
            if tmp.shape[0] == 0:
                continue
            comp_geneset = set(row['genes']) # get the set of genes 
            descendent = row['descendent'] # get the list of descendent nodes 
            for tmp_comp, tmp_row in tmp.iterrows():
                if tmp_comp in descendent: # skip if already in descendent
                    continue
                tmp_comp_geneset = set(tmp_row['genes'])
                # Check if satisfy ci_thre
                if len(comp_geneset.intersection(tmp_comp_geneset))/tmp_row['tsize'] >= ci_thre: #intersection of two components divided by the term size of the smaller component
                    # Check if child having higher weight than parent
                    # if cluster_weight[comp] < cluster_weight[tmp_comp]: ## do not have weight in hidef 
                    print('{} is contained in {} with a CI bigger than threshold, add edge between'.format(tmp_comp, comp))
                    G.add_edge(comp, tmp_comp, type='default')
                    clear = False
                    descendent += tmp_row['descendent']
        # Further clean up using networkx to remove shortcut edges 
        clean_shortcut(G)
        # Update variables
        n_iter += 1
        if clear:
            iterate = False
    if n_iter == 2:
        modified = False
    else:
        modified = True
    return modified
    
def merge_parent_child(G, hiergeneset, ji_thre):
    # Delete child term if highly similar with parent term
    # One parent-child relationship at a time to avoid complicacies involved in potential long tail
    print('... start removing highly similar parent-child relationship')
    similar = True
    merged = False
    while similar:
        clear = True
        edge_df = to_pandas_dataframe(G)
        ts_df = get_termStats(G, hiergeneset)
        default_edge = edge_df[edge_df['type'] == 'default'] # edges
        for idx, row in default_edge.iterrows():
            if jaccard(ts_df.loc[row['source']]['genes'], ts_df.loc[row['target']]['genes']) >= ji_thre:
                print('# Cluster pair {}->{} failed Jaccard, removing cluster {}'.format(row['source'], row['target'], 
                                                                                         row['target']))
                clear = False
                merged = True
                parents = edge_df[edge_df['target'] == row['target']]['source'].values
                children = edge_df[edge_df['source'] == row['target']]['target'].values
                # Remove all parent->node edges
                for pnode in parents:
                    G.remove_edge(pnode, row['target'])
                for child_node in children:
                    etype = G[row['target']][child_node]['type']
                    # Remove all node->child edges
                    G.remove_edge(row['target'], child_node)
                    # Add all parent->child edges
                    for pnode in parents:
                        G.add_edge(pnode, child_node, type=etype)
                # Remove target node
                G.remove_node(row['target'])
                break
        if clear:
            similar = False
    # Clean up shortcuts introduced during node deleteing process
    clean_shortcut(G)
    return merged

def collapse_redundant(G, hiergeneset, min_diff):
    # Delete child term if highly similar with parent term
    # One parent-child relationship at a time to avoid complicacies involved in potential long tail
    print('... start removing highly redundant systems')
    while True:
        edge_df = to_pandas_dataframe(G)
        ts_df = get_termStats(G, hiergeneset)
        default_edge = edge_df[edge_df['type'] == 'default']
        to_collapse = []
        for idx, row in default_edge.iterrows():
            parentSys, childSys, _ = row.values
            if ts_df.loc[parentSys]['tsize'] - ts_df.loc[childSys]['tsize'] < min_diff:
                to_collapse.append([parentSys, childSys])
        if len(to_collapse) == 0:
            print('nothing to collapse')
            return
        to_collapse = pd.DataFrame(to_collapse, columns=['parent', 'child'])
        # print(to_collapse)
        # cidx = to_collapse['weight'].idxmin()
        deleteSys = to_collapse.loc['child']
        print('# Cluster pair {}->{} highly redundant, removing cluster {}'.format(to_collapse.loc['parent'], 
                                                                                   to_collapse.loc['child'], 
                                                                                   deleteSys))
        parents = edge_df[edge_df['target'] == deleteSys]['source'].values
        children = edge_df[edge_df['source'] == deleteSys]['target'].values
        # Remove all parent->node edges
        for pnode in parents:
            G.remove_edge(pnode, deleteSys)
        for child_node in children:
            etype = G[deleteSys][child_node]['type']
            # Remove all node->child edges
            G.remove_edge(deleteSys, child_node)
            # Add all parent->child edges
            for pnode in parents:
                G.add_edge(pnode, child_node, type=etype)
        # Remove target node
        G.remove_node(deleteSys)

parser = argparse.ArgumentParser()
parser.add_argument('--outprefix', help='output_dir/file_prefix for the output file')
parser.add_argument('--ci_thre', type=float, default=0.75, help='Containment index threshold')
parser.add_argument('--ji_thre', type=float, default=0.9, 
                    help='Jaccard index threshold for merging similar clusters')
parser.add_argument('--minSystemSize', type=int, default=4, 
                    help='Minimum number of proteins requiring each system to have.')
parser.add_argument('--path_to_alignOntology', default = '/cellar/users/mhu/MuSIC/ddot/ddot/alignOntology', help='Full path to alignOntology.')
parser.add_argument('--min_diff', type=int, default=1, help='Minimum difference in number of proteins for every parent-child pair.')
args = parser.parse_args()

outprefix = args.outprefix
minSystemSize = args.minSystemSize




ci_thre = args.ci_thre
ji_thre = args.ji_thre
print('Containment index threshold: {}'.format(ci_thre))
print('Jaccard index threshold: {}'.format(ji_thre))

f = outprefix
ont = create_ontology(outprefix, minSystemSize)
hiergeneset = set(ont[ont['type'] == 'gene']['child'].values)

G = nx.from_pandas_edgelist(ont, source='parent', target='child', edge_attr='type', create_using=nx.DiGraph())

if not nx.is_directed_acyclic_graph(G):
    raise ValueError('Input hierarchy is not DAG!')

while True:
    modified = reorganize(G, hiergeneset, ci_thre)
    merged = merge_parent_child(G, hiergeneset, ji_thre)
    if not modified and not merged:
        break

collapse_redundant(G, hiergeneset, args.min_diff)
# Output as ddot edge file
clean_shortcut(G)
edge_df = to_pandas_dataframe(G)
edge_df.to_csv('{}_pruned.ont'.format(outprefix), header=False, index=False, sep='\t')
run_termStats = '{}/ontologyTermStats {} genes > {}'.format(args.path_to_alignOntology, 
                                                  '{}_pruned.ont'.format(outprefix), 
                                                  '{}_pruned.nodes'.format(outprefix))
os.system(run_termStats)

## step to clean the termStats file and make it the same format as hidef output 
nodes = pd.read_csv(outprefix+'_pruned.nodes', sep = '\t', header = None)
nodes.columns = ['terms', 'tsize', 'genes']
nodes = nodes.set_index('terms')
cleaned= []
for i, rows in nodes.iterrows():
    genes = sorted(rows['genes'].split(',')[:-1])
    cleaned.append(' '.join(genes))

nodes['genes'] = cleaned
nodes['logsize'] = [math.log2(x) for x in nodes['size']] ## add logsize to the nodes file (to read in cytoscape) 
nodes.to_csv(outprefix+'_pruned.nodes', header=False, sep='\t')


edges = edge_df[edge_df['type']=='default']# create the hidef format edges
edges.to_csv(outprefix+'_pruned.edges',sep = '\t', header=None,index=None) 

print(f'Number of edges is {len(edges)}, number of nodes are {len(nodes)}')
print('=== finished mature_hier_structure.py ====')
