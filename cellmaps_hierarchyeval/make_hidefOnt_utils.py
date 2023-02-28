# Utils for making hidef ontology from nodes and edges file 
import pandas as pd
import numpy as np
from collections import defaultdict
import argparse

### One step making ontology from the hidef Hierarchy 
## make a long table of (parent, child, edgetype)

def create_ontology(path, minTermsize= 4):
    ## read node table 
    node_table = pd.read_csv(path+'.nodes', header =None, sep ='\t')
    # print(node_table.head())
    # print(len(node_table))
    node_table.columns = ['terms', 'tsize', 'genes', 'stability']
    # node_table.genes = [g.split(' ') for g in node_table.genes] # change the genes column to list of genes 
    # print(node_table.head())
    node_table_filtered = node_table[node_table['tsize']>=minTermsize]
    # print(node_table[node_table['Num_genes']<minTermsize]['Node'])
    print(len(node_table_filtered))
    node_list = set(node_table_filtered['terms'] ) ## get the list of nodes 
    # print(node_list)

    new_rows = []
    for i, row in node_table_filtered.iterrows():
        for c in row['genes']:
            new_rows.append([row['terms'], c]) # for each gene, give the corresponding cluster number 

    node_table_final = pd.DataFrame(new_rows, columns = ['parent', 'child'])
    node_table_final['type'] = 'gene'
    
    ### Read edge table 
    edge_table = pd.read_csv(path+'.edges', header=None, sep ='\t')
    edge_table.columns = ['parent', 'child', 'type']
    # print(len(edge_table))
    #filtered so that only nodes bigger >=4 will be in here 
    edge_table_filtered = edge_table[edge_table[['parent', 'child']].isin(list(node_list)).all(axis=1)]
    
    

    add_gene = []
    genecount = 0
    ## Find leaves and get their genes 
    leaves =  (set(edge_table_filtered['child'].unique()))-set(edge_table_filtered['parent'].unique())
    # print(leaves)
    for leaf in leaves:
        genes = get_genes(node_table_filtered, leaf)
        genecount +=len(genes)
        # print(leaf, len(genes))
        for gene in genes:
            add_gene.append([leaf, gene, 'gene'])
    # print(len(add_gene))        
    
    parent_to_child = edge_table_filtered.groupby('parent')['child'].apply(list)  # group the parents
    # print(parent_to_child)
    for parent, children in parent_to_child.iteritems(): 
        parent_genes = get_genes(node_table_filtered, parent)
        # print(parent_genes)
        child_genes = []
        for child in children:
            genes = get_genes(node_table_filtered, child)
            # print(f'length of child {child} is {len(genes)}')
            child_genes = set(child_genes).union(set(genes))
        # print(f'full length of child genes is {len(child_genes)}')
        only_parent= set(parent_genes)-set(child_genes)## genes only in parent did not pass to child
        genecount +=len(only_parent)
        # print(parent, len(only_parent))
        if len(only_parent) >=1:
            for gene in only_parent:
                add_gene.append([parent, gene, 'gene'])
    # print(len(add_gene))
    # print(genecount)
    add_rows = pd.DataFrame(add_gene, columns=['parent', 'child', 'type'])
    final_df = edge_table_filtered.append(add_rows)
    print(len(final_df))
    return final_df

#get the genes from the nodes file, split each gene by the space
def get_genes(nodes_df, term):
    genes = nodes_df.loc[nodes_df.terms == term]['genes']
    gene_list =[]
    for g in genes:
        gene_list.extend(g.split(' '))
    return gene_list


## Load Edge Files and create ontology
def read_edge_table(path):
    edge_table = pd.read_csv(path, header=None, sep ='\t')
    edge_table.columns = ['Parent', 'Child', 'EdgeType']
    return edge_table

## Load Node Files and make a long table of (parent, child, edgetype)
def read_node_table(path):
    node_table = pd.read_csv(path, header =None, sep ='\t')
    #print(node_table.head())
    cd_list =[]
    node_table.columns = ['Node', 'Num_genes', 'Genes', 'Persistence']
    for i, row in node_table.iterrows():
        cd = row['Genes']
        all_cds = cd.split(' ')
        cd_list.append(all_cds)
    ## each gene will be in one row
    node_table['Genes_list'] = cd_list

    new_rows = []
    for i, row in node_table.iterrows():
        for c in row['Genes_list']:
            new_rows.append([row['Node'], c]) # for each gene, give the corresponding cluster number 

    node_table_final = pd.DataFrame(new_rows, columns = ['Parent', 'Child'])
    node_table_final['EdgeType'] = 'gene'
    
    return node_table_final

## clean up the table, remove duplicated genes 
def clean_table(df):
    only_nodes = df[df['EdgeType'] == 'default']
    only_genes = df[df['EdgeType'] == 'gene']
    
    ## for nodes
    node_to_gene = only_genes.groupby('Parent')['Child'].apply(list) # each cluster, genes inside
    gene_to_node = only_genes.groupby('Child')['Parent'].apply(list) # each gene, clusters they are present

    genes = list(only_genes['Child'].value_counts().keys()) # list of all genes
    print(len(genes))
    
    ## for edges 
    parent_to_child = only_nodes.groupby('Parent')['Child'].apply(list) # each cluster, the child cluster(s) they connected to
    child_to_parent = only_nodes.groupby('Child')['Parent'].apply(list) # each cluster, the parent cluster(s) they connected to

    nodes = set(only_nodes['Parent'].value_counts().keys()).union(set(only_nodes['Child'].value_counts().keys())) # get the union between parent clusters and child clusters (all clusters without counting duplicates)
    nodes_priority = defaultdict(int)

    only_genes.loc[:,'Priority'] = only_genes['Parent'].map(lambda x: get_depth(x, 0, child_to_parent)) # get the depth of each child

    keep_rows = []
    for g in genes:
        keep_rows.append(only_genes[only_genes['Child'] == g].sort_values(by='Priority', ascending=False).iloc[0]) # sort genes based on depth (choose the top)

    kept_genes = pd.DataFrame(keep_rows)
    kept_genes = kept_genes[['Parent', 'Child', 'EdgeType']]
    final_df = only_nodes.append(kept_genes)
    
    return final_df

def get_depth(node, count, child_to_parent):
    if node in child_to_parent.keys():
        max_d = get_depth(child_to_parent[node][0], count+1, child_to_parent)
        #save_n = child_to_parent[node][0]
        for n in child_to_parent[node]:
            d = get_depth(n, count+1, child_to_parent)
            if d> max_d:
                max_d = d
                
        return max_d
    else:
        return count

