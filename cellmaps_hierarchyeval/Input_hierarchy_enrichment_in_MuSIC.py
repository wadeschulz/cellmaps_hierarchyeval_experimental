import pandas as pd
from tqdm import tqdm
import sys
import os
import argparse
from cellmaps_utils.music_tools import load_obj, upper_tri_values
from statsmodels.stats.multitest import multipletests
from scipy.stats import hypergeom
import json

## calculate jaccard score
def jaccard(setA, setB):
    return len(setA.intersection(setB)) / len(setA.union(setB))

## function to load hidef hierarchy nodes and edges 
        
def load_hidef(f, minTermSize, w_root = False):
    '''
    load the hidef nodes file, add column names, and prefiltering
    args:
    f -- full path of the hidef nodes file 
    minTermSize -- term size threshold
    w_root -- whether to keep the root term, default= False
    
    return:
    1. filtered nodes table by minTermsize, annotate columns and index is the cluster name 
    2. hierarchy genes -- all genes in the root node (all genes in the hierarchy)
    3. filtered edges table, removing edges with the node that is smaller than minTermsize
    '''
    df = pd.read_table(f, header=None) ## load the input nodes 

    if len(df) > 800:
        print('=== Hierarchy too big! ===')
        sys.exit()
    
    df.columns = ['term', 'tsize', 'genes', 'stability']

    root_size = df['tsize'].max()

    hierarchygenes = df[df['tsize'] == root_size]['genes'].values[0].split(' ')  ## select the root node and collect all genes there (all genes included in the map)
    print(f'number of hierarchy genes = {len(hierarchygenes)}')
    
    node_rm = None
    if len(df[df['tsize']<minTermSize]) >0:
        node_rm = df[df['tsize']<minTermSize]['term']#the node need to be removed
        print('{} is smaller than minimum term size, Removed'.format(node_rm))  


    if w_root:
        df = df[df['tsize'] >= minTermSize]
    else:
        df = df[(df['tsize'] >= minTermSize) & (df['tsize'] < root_size)]

    if df.shape[0] == 0:
        print('=== No system left after size filter ===')
        sys.exit()

    df.set_index('term', inplace=True)
    print(f'Number of hierarchy nodes: {len(df)}')
    
    ### Load edge file
    edges_f = f[:-5] + 'edges'
    edges = pd.read_table(edges_f, header=None)
    if node_rm:
        edges = edges[~edges[[0, 1]].isin(list(node_rm)).any(axis=1)]
    
    return df, hierarchygenes, edges 


## analyze MuSIC hierarchy genes enriched in communities detected in DenseNet, Node2Vec and Bioplex hierarchies 
## default fdr threshold is 0.01 and the JI threshold is 0 

def run_hypergeo_enrichr(ont_ts, hierarchygenes, ref_file, fdr_thre=0.01, ji_thre = 0, minCompSize=4):  
    '''
    ont_ts: the result hierarchy from the community detection 
    hierarchy genes: total number of genes in the root of the hierarchy 
    ref_file: the reference cellular component/protein complexes table 
    fdr_thre: the fdr cutoff to be collected in the enriched terms (default = 0.01)
    ji_thre: the threshold for Jaccard Idex (default = 0, will set the threshold in the organization step)
    minCompSize: the smallest component size to be considered for the enrichment (default = 4)
    '''
    
    M = len(hierarchygenes)
    
    track = 0
    ref_df = pd.DataFrame(index=ont_ts.index, columns=ref_file.index, dtype=float)
    ji_df = pd.DataFrame(index=ont_ts.index, columns=ref_file.index, dtype=float)
    genecount_df = pd.DataFrame(index=ont_ts.index, columns=ref_file.index, dtype=float)
    for ont_comp, ont_row in tqdm(ont_ts.iterrows(), total=ont_ts.shape[0]):
        track += 1
        ont_comp_genes = ont_row['genes'].split(' ')
        ont_comp_genes = [x for x in ont_comp_genes if len(x) > 1] #in case theres a weird comma
        n = ont_row['tsize']
        for comp, row in ref_file.iterrows():
            comp_genes = row['genes'].split(' ')
            comp_genes = [x for x in comp_genes if len(x) > 1] #in case theres a weird comma
            N = len(set(comp_genes).intersection(set(hierarchygenes)))
            x = len(set(ont_comp_genes).intersection(set(comp_genes)))
            ref_df.at[ont_comp, comp] = hypergeom.sf(x - 1, M, n, N) # calculare the hypergeometric distribution 
            ji_df.at[ont_comp, comp] = jaccard(set(ont_comp_genes), set(comp_genes)) # calculate the jaccard score
            genecount_df.at[ont_comp, comp] = x
    fdr = multipletests(ref_df.values.flatten(), method='fdr_bh')[1]
    fdr_df = pd.DataFrame(fdr.reshape(ref_df.shape), index=ref_df.index, columns=ref_df.columns, dtype=float)
    for comp, row in fdr_df.iterrows():
        ont_ts.at[comp, 'enriched'] = 0
        ont_ts.at[comp, 'terms'] = ""
        ont_ts.at[comp, 'gene_counts'] = ""
        ont_ts.at[comp, 'adjPvalue'] = ""
        ont_ts.at[comp, 'ji_indexes'] = ""

        tmp = row[row <= fdr_thre]
        
        gene_counts = []
        terms = []
        ji = []
        pvalues = []
        for term in tmp.index:
            if ji_df.at[comp, term] >= ji_thre:
                gene_counts.append(genecount_df.at[comp, term])
                terms.append(term)
                ji.append(ji_df.at[comp,term])
                pvalues.append(tmp[term])
        tmp_df = pd.DataFrame({'gene':gene_counts,'ji': ji, 'terms': terms, 'pvalue': pvalues})
        tmp_df.sort_values(ascending=False, by='ji', inplace=True)
        ont_ts.at[comp, 'enriched'] = len(terms)
        ont_ts.at[comp, 'terms'] = "; ".join(tmp_df['terms'])
        ont_ts.at[comp, 'gene_counts'] =  "; ".join(tmp_df['gene'].astype(str))
        ont_ts.at[comp, 'adjPvalue'] =  "; ".join(tmp_df['pvalue'].astype(str))
        ont_ts.at[comp, 'ji_indexes'] =  "; ".join(tmp_df['ji'].astype(str))

    for comp in ont_ts.index:
        enrich_terms = ji_df.columns[ji_df.loc[comp] == 1].tolist()
        for term in enrich_terms:
            ont_ts.at[comp, 'enriched'] = ont_ts.at[comp, 'enriched'] + 1
            ont_ts.at[comp, 'terms'] = term + "; " + ont_ts.at[comp, 'terms']
            ont_ts.at[comp, 'ji_indexes'] = str(ji_df.at[comp, term]) + "; " + ont_ts.at[comp, 'ji_indexes'] 
            ont_ts.at[comp, 'adjPvalue'] = str(fdr_df.at[comp, term]) + "; " + ont_ts.at[comp, 'adjPvalue']
            ont_ts.at[comp, 'gene_counts'] = str(genecount_df.at[comp, term]) + "; " + ont_ts.at[comp, 'gene_counts']
    print(f'Number of nodes: {len(ont_ts)}' )
    print('... finished running geneset enrichment')
    return ont_ts

###edge enrichment test 
## Load data from PPI edges
def load_data(filedir, rootgenes):
    df = load_obj(filedir)
    ref_genes = df.index.values
    intersect_genes = list(set(ref_genes).intersection(set(rootgenes)))
    print(len(intersect_genes))
    df = df.loc[intersect_genes, intersect_genes]
    return df, intersect_genes

def num_comb(x): #calculate all possible edges between the list of genes (x^2-x)/2 
    # minus x because we do not need the edges between the same genes 
    # divided by 2 because edges between A-B and B-A are the same edge, do not want to double count 
    return x*(x-1)/2

### function for get the edges that are in the enriched subsytems and the pvalues
def get_count_pval(ref_matrix, ref_genes, hier_matrix):
    M = len(upper_tri_values(ref_matrix))
    total = upper_tri_values(ref_matrix).sum()

    ref_x_termgenes = list(set(termgenes).intersection(set(ref_genes)))  ## filter for ref genes that are in the parent nodes
    if len(ref_x_termgenes) < 2: # if only one gene enriched in the nodes then p val is high, and not count as enriched node
        count_value = 0
        ref_pval = 1
        ref_unique_count = 0
    else:
        count_value = 0
        for i in range(len(ref_x_termgenes)-1):
            ga = ref_x_termgenes[i]
            for j in range(i+1, len(ref_x_termgenes)):
                gb = ref_x_termgenes[j]
                # check if edge is specific to this system...
                if children_edges.at[ga,gb] == 1: ## is there an edge in the children nodes? if yes, do not count
                    continue
                if ref_matrix.at[ga, gb] == 1: ## is the edge in the hcm system? if yes, add to the edge count
                    count_value += 1
                    hier_matrix.at[ga, gb] = 1
                    hier_matrix.at[gb, ga] = 1

        N = num_comb(len(ref_x_termgenes)) - sum(upper_tri_values(children_edges.loc[ref_x_termgenes,ref_x_termgenes])) ## number all possible edges of the overlapped genes in the term - number of edges that are in the child nodes
        ref_pval= hypergeom.sf(count_value -1, M, total, N) # hypermeometric stats
        return count_value, sum(upper_tri_values(hier_matrix)), ref_pval # return edge counts, unique edges and p value from hyermeometric 


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze each system in the given hierarchy.')
    parser.add_argument('--infname', help='Full path of the input nodes file.')
    parser.add_argument('--outprefix', help='Prefix of output file path. Format: path_to_outdir/file_identifier')
    parser.add_argument('--refhier', required=True, nargs='+', 
                        help='Full path of the reference hidef hierarchy nodes. Can be multiple references paths')
    parser.add_argument('--w_root', action='store_true', help='Do analysis for root term.')
    parser.add_argument('--minTermSize', default=4, type=int)
    parser.add_argument('--FDRthre', default=0.01, type=float)
    args = parser.parse_args()


    ### load the nodes data ###
    f = args.infname
    minTermSize = args.minTermSize
    outprefix = args.outprefix
    fdrthre = args.FDRthre
    refdir = args.refhier

    if not os.path.exists(f):
        print(f)
        raise ValueError('Input termStats file does not exist!')


    if os.path.getsize(f) == 0:
        print('=== No term in hierarchy! ===')
        sys.exit()

    df, hierarchygenes, edges = load_hidef(f, minTermSize, w_root =args.w_root)
    
    #### hypergeometric test on clusters in DenseNet, Node2Vec, and Bioplex Hierarchy 
    for ref in refdir:
        label = ref.rsplit('/',1)[-1].rsplit('_')[2].rsplit('.')[0] ##Get the type of reference hierarchy 
        print(f'type of ref is {label}')
        
        ref_file, _ , _ = load_hidef(ref, minTermSize, w_root =args.w_root) #load the input ref hidef hierarchy
        
        hypergeom_result = run_hypergeo_enrichr(df.copy(), hierarchygenes, ref_file, fdr_thre=fdrthre, minCompSize=minTermSize)

        df[f'hypergeom_{label}_enrich'] = hypergeom_result['enriched']
        df[f'hypergeom_{label}_terms'] = hypergeom_result['terms']
        df[f'{label}_gene_counts'] = hypergeom_result['gene_counts']
        df[f'hypergeom_{label}_adjPvalues'] = hypergeom_result['adjPvalue']
        df[f'{label}_ji_indexes'] = hypergeom_result['ji_indexes']
        
    #### load the u2os HPA locations file 
    with open("/cellar/users/mhu/MuSIC/Resources/u2os_loc.txt", "r") as f:
        geneloc =json.loads(f.read()) # load the location dictionary
    loc_df=pd.DataFrame(geneloc.keys(), columns = ['Locations']) #the keys are the locations
    loc_df['genes'] = [x for x in geneloc.values()] #values are genes
    loc_df['tsize'] = loc_df.genes.apply(lambda x: len(x))# get the number of genes in each location as the term size
    loc_df['genes'] = loc_df.genes.apply(lambda x: ' '.join(x)) # join the genes with space (for consistency)
    loc_df = loc_df[loc_df['tsize'] >=minTermSize] # remove locations with fewer than mintermsize genes
    loc_df.set_index('Locations', inplace=True) #set the location (term) as the index
    
    #### hypergeometric test on the HPA locations 
    hypergeom_hpa_enrich_result = run_hypergeo_enrichr(df.copy(), hierarchygenes, loc_df, fdr_thre=fdrthre, minCompSize=minTermSize)
    df['hypergeom_hpa_enrich'] = hypergeom_hpa_enrich_result['enriched']
    df['hypergeom_hpa_terms'] = hypergeom_hpa_enrich_result['terms']
    df['hpa_gene_counts'] = hypergeom_hpa_enrich_result['gene_counts']
    df['hypergeom_hpa_adjPvalues'] = hypergeom_hpa_enrich_result['adjPvalue']
    df['hpa_ji_indexes'] = hypergeom_hpa_enrich_result['ji_indexes']
    
    
    print('... finished running enrichment analysis')

    # save the dataframe and the genes that are intersect with the ones in the hierarchy root
    bioplex, bioplex_genes = load_data('/cellar/users/mhu/MuSIC/U2OS/resources/U2OS_bioplex_012623.M.pkl', hierarchygenes)

    print('... finished loading PPI data')


    ##run the edge enrichment analysis

    bioplex_count = []
    bioplex_unique_count = []
    bioplex_pvalue = []

    for idx, row in tqdm(df.iterrows(), total=df.shape[0]):
        termgenes = row['genes'].split(' ') #list all genes in the parent nodes 

        children = edges[edges[0] == idx] # children is all the edges the system have in the edges file 
        children_edges = pd.DataFrame(0, index=termgenes, columns=termgenes, dtype=int)
        for child in children[1].values: #children[0] is the parent children[1] are the child nodes
            if child in df.index.values: #b/c node table is filtered by the term size >= 4 but edges may still contain small components 
                child_genes = df.loc[child, 'genes'].split(' ') # find genes in each child node
                # print(child_genes)
                for i in range(len(child_genes)-1):
                    ga = child_genes[i]
                    if ga not in termgenes: #skip the genes that is not in parent
                        continue
                    for j in range(i+1, len(child_genes)):
                        gb = child_genes[j]
                        if gb not in termgenes:
                            continue
                        children_edges.at[ga, gb] = 1  ## get an edge between each pair of genes (that is also in the parents)
                        children_edges.at[gb, ga] = 1

        #start with smallest nodes to get children....
        hier_bioplex = pd.DataFrame(0, index=bioplex_genes, columns=bioplex_genes, dtype=int)

        bioplex_edges, bioplex_unique_edges, bioplex_pval = get_count_pval(bioplex, bioplex_genes, hier_bioplex)
        bioplex_count.append(bioplex_edges)
        bioplex_unique_count.append(bioplex_unique_edges)
        bioplex_pvalue.append(bioplex_pval)


    #Bioplex
    df['bioplex_count'] = bioplex_count
    df['bioplex_unique_count'] = bioplex_unique_count
    df['bioplex_pvalue'] = bioplex_pvalue
    df['bioplex_adjp'] = multipletests(bioplex_pvalue, method='fdr_bh')[1]


    print('... finished running edge enrichment analysis')

    
    root_file_label = 'noRoot'
    if args.w_root:
        root_file_label = 'wRoot'

    # from datetime import date
    # day = date.today().strftime("%m%d%Y")
    df.to_csv('{}.{}.enrich_input.csv'.format(outprefix, root_file_label), sep='\t')
    save_obj(df, '{}.{}.enrich_input.pkl'.format(outprefix, root_file_label))
    print('=== finished analyze_hidef_output ===')
