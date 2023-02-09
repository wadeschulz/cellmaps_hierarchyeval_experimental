import pandas as pd
import numpy as np
import os
from cellmaps_utils.music_utils import load_obj
from glob import glob
from tqdm import tqdm
import re

##Functions to use for enrichment analysis 

#return the hierarchy depth
def get_depth(node, count, child_to_parent):
    if node in child_to_parent.keys():
        max_d = get_depth(child_to_parent[node][0], count+1, child_to_parent)
        for n in child_to_parent[node]:
            d = get_depth(n, count+1, child_to_parent) ##some nodes connect to the same level 
            # print(d)
            if d< max_d: ##only consider the shortest route
                max_d = d                
        return max_d
    else:
        return count

##########the following functions are no longer used because it is depend on the cluster name that hidef spits out######
## return maximum depth    
# def get_maxdepth(hier):
#     cluster_num = [int(re.findall("([0-9]+)", x)[0]) for x in hier.index]
#     depth = max(cluster_num )
#     return depth

# #return the hierarchy median depth
# def get_mediandepth(edges_df):
#     edges_df.columns = ['parent', 'child', 'type']
#     leaf_nodes = list((set(edges_df['child'].unique()))-set(edges_df['parent'].unique()))
#     depth = [int(re.findall("([0-9]+)", x)[0]) for x in leaf_nodes]
#     med_depth = statistics.median(depth)
#     return med_depth
    
#return the hierarchy breadth
def get_maxbreadth(hier):
    nodes_num = [int(re.findall("([0-9]+)", x)[1]) for x in hier.index]
    breadth = max(nodes_num)
    return breadth

# return the number of genes that are assigned to clusters that are below level one 
def get_numgenes_belowone(hier):
    idx = [x.startswith('Cluster1') == False | x.startswith('Cluster0') == False for x in hier.index] #get rows that is below 1 
    belowone = hier[idx]
    genes = []
    for idx, row in belowone.iterrows():
        genes.append(row['genes'].split(' '))
    all_genes = sum(genes, [])
    return len(set(all_genes)) ## number of unique genes in below level 1

# calculate the number of child a cluster has 
def num_child(edges, hier, cluster):
    edges.columns= ['parent', 'child', 'type']
    parent_to_child = edges.groupby('parent')['child'].apply(list)
    idx = [x.startswith(cluster) for x in hier.index]
    levelone = hier[idx].index.values
    num_child =[]
    for term in levelone: 
        if term in parent_to_child.keys():
            num_child.append(len(parent_to_child[term]))
        else:
            num_child.append(0)
    return np.mean(num_child)

# find number of enriched ref components (large vs small)
def enrich_small_large(enrich_type, df, ref_ts, known_comp = None, large = 80, fdr = 0.01, ji = 0.2):
    term_term_mapping = []
    for index,row in df.iterrows():
        pvalues = str(row[f'hypergeom_{enrich_type}_adjPvalues']).split('; ')
        jaccard_indexes = str(row[f'{enrich_type}_ji_indexes']).split('; ')
        hypergeom_terms = str(row[f'hypergeom_{enrich_type}_terms']).split('; ')
        for i in range(len(pvalues) - 1):
            term_term_mapping.append((index, pvalues[i], hypergeom_terms[i], jaccard_indexes[i]))

    term_term_mapping = pd.DataFrame(term_term_mapping, columns=['term','pval', 'gs_term', 'ji'])
    term_term_mapping = term_term_mapping.merge(ref_ts['tsize'], left_on='gs_term', right_on=ref_ts.index, how='left')
    num_enriched = 0
    num_sig = 0
    num_enriched_small = 0
    num_enriched_large = 0
    num_sig_small = 0
    num_sig_large = 0
    taken_gs_terms = []
    taken_terms = []

    for index,row in term_term_mapping.sort_values(by='ji', ascending=False).iterrows():
        if (row['gs_term'] not in taken_gs_terms) & (row['term'] not in taken_terms):
            taken_gs_terms.append(row['gs_term'])
            taken_terms.append(row['term'])
            num_enriched += 1
            if (float(row['pval']) <= fdr):#set p < fdr thre as the enriched term
                if row['tsize'] >= large:
                    num_enriched_large += 1
                else:
                    num_enriched_small += 1
                if (float(row['pval']) <= fdr) & (float(row['ji']) >= ji):
                    num_sig +=1
                    if row['tsize'] >= 80:
                        num_sig_large += 1
                    else:
                        num_sig_small += 1
    return num_enriched, num_sig, num_enriched_small, num_sig_small, num_enriched_large, num_sig_large

## Enrichment for known (Gold standard GO terms) -----As for Dec 2022, these are manually put in 

'''
Proteasome =  'GO:0000502'   (size 45)
Ribosome = 'GO:0005840'    (size 95)
Nuclear Speck = 'GO:0016607' (size 114)
mitochondrion = ‘GO:0005739’ (size 591)
nucleus = ‘GO:0005634’ (size 2555)
Nuclear pore = 'GO:0005643' (size 45)
cytoplasm = 'GO:0005737' (size 3650)
ER = ‘GO:0005783’  (size 522)
cytosol = ‘GO:0005829’ (size 1724)
golgi_apparatus = ‘GO:0005794’ (size 413)
lysosome = 'GO:0005764' (size 241)
membrane = 'GO:0016020' (size 2270)

'''

known_comp = ['GO:0000502', 'GO:0005840', 'GO:0016607', 'GO:0005739','GO:0005634','GO:0005643', 'GO:0005737', 'GO:0005783', 'GO:0005829','GO:0005794', 'GO:0005764','GO:0016020']
def enrich_GS(enrich_type, df, known_comp, fdr = 0.01, ji = 0.2):
        term_term_mapping = []
        for index,row in df.iterrows():
            pvalues = str(row[f'hypergeom_{enrich_type}_adjPvalues']).split('; ')
            jaccard_indexes = str(row[f'{enrich_type}_ji_indexes']).split('; ')
            hypergeom_terms = str(row[f'hypergeom_{enrich_type}_terms']).split('; ')
            for i in range(len(pvalues) - 1):
                term_term_mapping.append((index, pvalues[i], hypergeom_terms[i], jaccard_indexes[i]))

        term_term_mapping = pd.DataFrame(term_term_mapping, columns=['term','pval', 'gs_term', 'ji'])
        num_enriched = 0
        num_sig = 0
        taken_gs_terms = []
        taken_terms = []
        enriched_GO = dict.fromkeys(known_comp, 0)
        sig_enriched_GO = dict.fromkeys(known_comp, 0) 

        for index,row in term_term_mapping.sort_values(by='ji', ascending=False).iterrows():
            if (row['gs_term'] not in taken_gs_terms) & (row['term'] not in taken_terms):
                taken_gs_terms.append(row['gs_term'])
                taken_terms.append(row['term'])
                if (row['gs_term'] in known_comp)&(float(row['pval']) <= fdr): #set p < fdr thre as the enriched term
                    # enriched_GO[row['gs_term']]+=1
                    enriched_GO[row['gs_term']]=[row['pval'],row['ji']]
                    num_enriched +=1
                    if (float(row['pval']) <= fdr) & (float(row['ji']) >= ji):
                        sig_enriched_GO[row['gs_term']]+=1
                        num_sig +=1
        return num_enriched, num_sig, sig_enriched_GO

## number of enriched components regardless of size 
def enrich_eval(enrich_type, df, fdr = 0.01, ji = 0.2):
        term_term_mapping = []
        for index,row in df.iterrows():
            pvalues = str(row[f'hypergeom_{enrich_type}_adjPvalues']).split('; ')
            jaccard_indexes = str(row[f'{enrich_type}_ji_indexes']).split('; ')
            hypergeom_terms = str(row[f'hypergeom_{enrich_type}_terms']).split('; ')
            for i in range(len(pvalues) - 1):
                term_term_mapping.append((index, pvalues[i], hypergeom_terms[i], jaccard_indexes[i]))

        term_term_mapping = pd.DataFrame(term_term_mapping, columns=['term','pval', 'gs_term', 'ji'])
        num_enriched = 0
        num_sig = 0
        taken_gs_terms = []
        taken_terms = []

        for index,row in term_term_mapping.sort_values(by='ji', ascending=False).iterrows():
            if (row['gs_term'] not in taken_gs_terms) & (row['term'] not in taken_terms):
                taken_gs_terms.append(row['gs_term'])
                taken_terms.append(row['term'])
                if (float(row['pval']) <= fdr): #set p < fdr thre as the enriched term
                    num_enriched += 1
              
                    if (float(row['pval']) <= fdr) & (float(row['ji']) >= ji):
                        num_sig +=1
        return num_enriched, num_sig

## number of enriched edges for nodes smaller than a particular size 
def enrich_edge(enrich_type, df, large = 80, fdr = 0.01, ji = 0.2):
    enriched_comp = df[(df[f'{enrich_type}_count'] > 0) & (df[f'{enrich_type}_adjp'] <= fdr) & (df['tsize'] < large)] # significant components 
    num_enriched = len(df[(df[f'{enrich_type}_count'] > 0) & (df['tsize'] < large)])
    num_sig = len(enriched_comp) # number of the significant hidef components that is enriched 
    num_enrich_edges = sum(enriched_comp[f'{enrich_type}_count']) # number of edge counts that is enriched in the hidef hierarchy 
    num_unique_edges = sum(enriched_comp[f'{enrich_type}_unique_count'])
    # hcm_stats.at[idx, 'ratio_enriched'] = len(enriched_comp) / hs_df.shape[0] # number of component enriched in hcm vs total hidef hierarch
    return (num_enriched, num_sig, num_enrich_edges, num_unique_edges)


#### Function for picking one best enrichment term and annotate the node with it

def annotation_enrich(enrich_type, df, fdr = 0.01, ji = 0.2):
    '''
    used for picking one best enriched term and annotate the system 
    enrich_type: type of databased use for annotation, options: 'cc', 'corum'
    df: the hidef hierarchy nodes file after enrichment analysis
    fdr: cutoff for calling enriched 
    ji: JI for calling significant
    '''
    term_term_mapping = []
    for index,row in df.iterrows():
        pvalues = str(row[f'hypergeom_{enrich_type}_adjPvalues']).split('; ')
        jaccard_indexes = str(row[f'{enrich_type}_ji_indexes']).split('; ')
        hypergeom_terms = str(row[f'hypergeom_{enrich_type}_terms']).split('; ')
        if pvalues[0]:
            if (float(pvalues[0])<= fdr) & (float(jaccard_indexes[0])>=ji):
                term_term_mapping.append((index, hypergeom_terms[0], pvalues[0], jaccard_indexes[0]))
            else:
                term_term_mapping.append((index, '', '', ''))
        else: 
            term_term_mapping.append((index, '', '', ''))
    term_term_mapping = pd.DataFrame(term_term_mapping, columns=['term',f'{enrich_type}_term',f'{enrich_type}_pval', f'{enrich_type}_ji'])
    return term_term_mapping



####################function to analyze # of enrichments per hierarchy ##############################
def analyze_enrichment(prefix,  workdir,  analyze_day, refdir = '/cellar/users/mhu/MuSIC/U2OS/resources/',  minTermSize = 4,fdr= 0.01, ji= 0.2, large = 80):
    
    
    '''
    prefix: the prefix used for the output file 
    workdir: the directory where the inputs are and output will be stored 
    analyze_day: the date for running the gene set and edge enrichment analysis 
    refdir: the directory where the reference GO file is stored
    minTermSize: the minimum size of the node to consider(default = 4)
    fdr: the fdr cutoff to consider significantly enriched (default = 0.01)
    ji: the Jaccard Index cutoff to consider significantly enriched terms (only for GO and corum terms, default = 0.2)
    large: the enriched terms that are bigger than this value will be considered as large (default = 80) (used for separate large and small GO, and collect enriched edges in smaller nodes_
    
    '''
    cc_ts = pd.read_table(f'{refdir}cc_noHPA.5183.dropDupTerm.termStats', header = None, index_col = 0)
    cc_ts.columns = ['tsize', 'genes']
 
    cc_ts = cc_ts[cc_ts['tsize']>=minTermSize]

    common_df = []
    new_df = []
    fnames = glob('{}*.nodes'.format(workdir))
    index = [x.split('/')[-1][:-6] for x in fnames] 

    for filename in tqdm(index): 
        nodes_df = pd.read_csv(workdir + filename +'.nodes', sep = '\t', header = None) ## open nodes file
        edges_df = pd.read_csv(workdir + filename +'.edges', sep = '\t', header = None) ## open edges 
        num_nodes = len(nodes_df) # total number of nodes
        num_edges = len(edges_df) # total number of edges between clusters 
        edges_df.columns= ['parent', 'child', 'type']
        child_to_parent = edges_df.groupby('child')['parent'].apply(list) #list of child parent 
        leaves = (set(edges_df['child'].unique()))-set(edges_df['parent'].unique())

        # if 'RF' in filename:
        #     method = 'RF'
        # else: 
        #     method = label +re.search('5183_(.+?)_percthre', filename).group(1)
        # # print (method)
        # cutoffthre= re.search('percthre_(.+?).chi', filename).group(1) ##get the cutoff threshold from the name
        # maxres= filename.split('_')[-1] ## get the maxresolution in the name
        # stability= re.search('chi_(.+?).maxres', filename).group(1) #presistent threshold 
        # # print(filename, num_edges, num_nodes, cutoffthre, maxres, stability)

        if not os.path.exists('{}/{}.noRoot{}.pkl'.format(workdir, filename, analyze_day)):
            continue
        hs_df = load_obj('{}/{}.noRoot{}.pkl'.format(workdir, filename, analyze_day))
        # hs_df = hs_df[hs_df['tsize'] <= max_comp_size] # maximum component size 
    #     # max_depth = get_maxdepth(hs_df) ##calculate the depth of each hierarchy 
        depth = [get_depth(x, 0, child_to_parent) for x in leaves]
        med_depth =np.median(depth) ##calculate the median depth
           ##size of leaves
        leaves = leaves.intersection(hs_df.index.values)
        leaf_size= [hs_df.loc[x, 'tsize'] for x in leaves]
        avg_leaf_size = np.mean(leaf_size)
        # avg_leaf_size = np.mean([y / x for x, y in zip(depth, leaf_size)])
        max_breadth = get_maxbreadth(hs_df) ##calculate the maximun breadth of each hierarchy (do not want too deep or too wide hierarchy)
        num_below_level1 = get_numgenes_belowone(hs_df)
        num_child_level1 = num_child(edges_df, hs_df, 'Cluster1')
        ##first column is the file name, including the cutoffs, stability, maxres, algorithm, and coembedding method etc.
        common_df.append([filename, num_nodes, num_edges, med_depth, max_breadth, avg_leaf_size, num_below_level1, num_child_level1]) ## collect the informations in common
        # common_df.append([method, cutoffthre, maxres, stability, num_nodes, num_edges, med_depth, max_breadth, avg_leaf_size, num_below_level1, num_child_level1]) ## collect the informations in common

        temp = {}

        for ref in ['cc', 'corum', 'hcm', 'achilles_2std', 'bioplex', 'pcnet', 'opencell']:
            if ref =='cc': ## collect num of enriched cc, and corum system, significant and the ratio
                num_enriched, num_sig, num_enriched_small, num_sig_small, num_enriched_large, num_sig_large= enrich_small_large(ref, hs_df, cc_ts, large, fdr, ji)          
                frac_enrich =num_enriched/num_nodes ## ratio of enriched systems
                frac_sig = num_sig/num_nodes ## ratio of sig enriched
                num_large_cc = len(cc_ts[cc_ts.tsize>=large])
                num_small_cc = len(cc_ts[cc_ts.tsize<large])
                num_enriched_known, n_sig_known, sig_enriched_terms=enrich_GS(ref, hs_df,known_comp, fdr = 0.01, ji = 0.2)
                frac_sig_known = n_sig_known/len(known_comp)
                # frac_sig_small = num_sig_small/num_small_cc
                # frac_sig_large = num_sig_large/num_large_cc
                # small_efficiency = num_sig_small/num_nodes

                temp = {**temp, **{#f"Number of components enriched in {ref.upper()}(FDR <= {fdr})":num_enriched, 
                                #f"Frac of Enriched Systems ({ref.upper()})":frac_enrich, 
                                f"Number of Significant Systems ({ref.upper()}) with FDR <= {fdr}, JI >= {ji}": num_sig, 
                                f"Frac of Significant Systems ({ref.upper()})":frac_sig,
                                f"Number of Significant small Systems ({ref.upper()}) with FDR <= {fdr}, JI >= {ji}": num_sig_small, 
                                # f"Frac of Significant small Systems in {ref.upper()} small":frac_sig_small,
                                # f"Significant Enriched small {ref.upper()} Systems Efficiency":small_efficiency,
                                f"Number of Significant large Systems ({ref.upper()}) with FDR <= {fdr}, JI >= {ji}": num_sig_large, 
                                # f"Frac of Significant large Systems in {ref.upper()} large":frac_sig_large
                                f"Number of enriched well known Systems in {ref.upper()} (FDR <= {fdr})": num_enriched_known,
                                f"Number of Significant well known Systems ({ref.upper()}) with FDR <= {fdr}, JI >= {ji}": n_sig_known,
                                f"Significant well known Systems ({ref.upper()}) with FDR <= {fdr}, JI >= {ji}": sig_enriched_terms,
                                f"Frac of Significant large Systems in {ref.upper()} large":frac_sig_known
                          }}
            elif ref == 'corum':
                num_enriched, num_sig = enrich_eval(ref, hs_df, fdr, ji)          
                frac_enrich =num_enriched/num_nodes ## ratio of enriched systems
                frac_sig = num_sig/num_nodes ## ratio of sig enriched

                temp = {**temp, **{f"Number of components enriched in {ref.upper()} (FDR <= {fdr})":num_enriched, 
                                        f"Ratio of Enriched Systems ({ref.upper()})":frac_enrich, 
                                        f"Number of Significant Systems ({ref.upper()}) with FDR <= {fdr}, JI >= {ji}": num_sig, 
                                        f"Ratio of Significant Systems ({ref.upper()})":frac_sig}}
            else:
                num_enriched, num_sig, num_enrich_edges, num_unique_edges = enrich_edge(ref, hs_df,large ,fdr, ji)
                frac_sig = num_sig/num_nodes 
                if ref == 'achilles_2std':
                    temp = {**temp, **{#"Number of components enriched in DepMap":num_enriched, 
                                       # "Frac of Enriched Systems (DepMap)":frac_enrich, 
                                        f"Number of Significant Systems (DepMap) with FDR <= {fdr}, size < {large}": num_sig, 
                                        "Frac of Significant Systems (DepMap)":frac_sig,
                                        # f"Number of edges in DepMap enriched systems (size < {large})":num_enrich_edges,
                                        f"Number of unique edges in DepMap enriched systems (size < {large})": num_unique_edges}}
                else:
                    temp = {**temp, **{#f"Number of components enriched in {ref.upper()}":num_enriched, 
                                       # f"Frac of Enriched Systems ({ref.upper()})":frac_enrich, 
                                        f"Number of Significant Systems ({ref.upper()}) with FDR <= {fdr}": num_sig, 
                                        f"Frac of Significant Systems ({ref.upper()})":frac_sig,
                                        # f"Number of edges in {ref.upper()} enriched systems (size < {large})": num_enrich_edges,
                                        f"Number of unique edges in {ref.upper()} enriched systems (size < {large})": num_unique_edges}}




        new_df.append(temp)

        # print(f'------ Analyze {filename} enrichment DONE ------')
    print('DONE')
    common_df = pd.DataFrame(common_df, columns = ["Hierarchy name", "Number of Nodes", "Number of Edges", 
                                                 "MedianDepth", "MaxBreadth","Avg Leaf size", "Number of genes assigned to nodes with depth >= 2", 'Avg number of child per level one node'])
    # common_df = pd.DataFrame(common_df, columns = ["Method", "Percent Edge Cutoff", "MaxRes", "Presistent threshold", "Number of Nodes", "Number of Edges", 
    #                                              "MedianDepth", "MaxBreadth","Avg Leaf size", "Number of genes below level one of the hierarchy", 'Avg number of child per level one node'])
    new_df= pd.concat([common_df, pd.DataFrame(new_df)], axis = 1)
    sort_new_df = new_df.sort_values(by= ['Hierarchy name'])
    # sort_new_df = new_df.sort_values(by=["Method", 'Percent Edge Cutoff','MaxRes', 'Presistent threshold'])
    sort_new_df.to_csv(workdir + prefix+ 'hidef_enrichment_analysis.csv')
    return sort_new_df 


