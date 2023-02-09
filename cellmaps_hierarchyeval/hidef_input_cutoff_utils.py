import pandas as pd
import numpy as np
import pickle
from cellmaps_utils.music_tools import euclidean_similarity, cosine_similarity_scaled, load_obj
import matplotlib.pyplot as plt

import math

##Save similarity matrix 

def make_sim_matrix(X,savedir, prefix, sim_type = 'euc_sim'):
    '''
    X: the input embedding 
    savedir: path to save the file 
    prefix: the prefix to name the similarity matrix 
    sim_type: the similarity calculation between pairs (euc_sim, or cos_sim), default: 'euc_sim' 
    '''
    if sim_type == 'euc_sim':
    
        with open('{}/{}.{}.scaled.pkl'.format(savedir, prefix, sim_type), 'wb') as f:
            pickle.dump(euclidean_similarity(X), f)
    
    if sim_type == 'cos_sim':
        with open('{}/{}.{}.scaled.pkl'.format(savedir, prefix, sim_type), 'wb') as f:
            pickle.dump(cosine_similarity_scaled(X), f)
    

## convert similarity matrix to pairwise distances
def mat_to_pairs(sim_mat_dir, sort = True):
    '''
    function to convert similarity matrix to pairwise distances/similarities
    
    sim_mat_dir: directory of the input sim_mat
    sort: sort the pair-wise similarities by descending weight, type= boolean, default is True
    
    '''
    sim_mat = load_obj(sim_mat_dir)
    keep = np.triu(np.ones(sim_mat.shape)).astype(bool) ##keep the upper triangle  
    #keep
    sim_mat = sim_mat.where(keep)
    pair_mat = sim_mat.stack().reset_index().rename(columns={'level_0': 'GeneA', 'level_1': 'GeneB', 0: 'Weight'})
    pair_mat = pair_mat[pair_mat['GeneA'] != pair_mat['GeneB']]
    if sort: 
        pair_mat = pair_mat.sort_values('Weight', ascending=False) 
    return pair_mat

######Check how number of proteins changes and the percentage of edges changes at different weight cutoff#####

#function to do titration at different cutoff (record the number of proteins and percentage of edges) 
def tit_cutoff(sorted_df, cutoff):
    sorted_df.columns = ['geneA', 'geneB', 'weight']
    total_edges = len(sorted_df)
    num_genes = []
    perc_edges = []
    for i in cutoff:
        df_cutoff = sorted_df[sorted_df['weight'] > i][['geneA', 'geneB', 'weight']]
        gene_list = set(df_cutoff['geneA']).union(set(df_cutoff['geneB']))
        num_genes.append(len(gene_list))
        num_edges_keep = len(df_cutoff)
        perc_edges.append((num_edges_keep/total_edges)*100)
    df = pd.DataFrame({'Cutoffs': cutoff, 'Number of proteins': num_genes, 'Percentage of edges retained (%)': perc_edges})
    return df

# function to plot the graph to visualize the changes 
from datetime import date
def plot_gene_edge_loss(df, batch_type, savedir):
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax2 = ax1.twinx()
    ax1.plot('Cutoffs', 'Number of proteins', data = df , marker='o',  color = 'navy', linewidth=3, markersize=10)
    ax2.plot('Cutoffs', 'Percentage of edges retained (%)', data = df , marker='o',  color="coral", linewidth=3, markersize=10)

    ax1.set_xlabel('Weight Cutoffs', fontsize=18)
    ax1.set_ylabel('Number of proteins remaining', color='navy', fontsize=18)
    ax2.set_ylabel('Percentage of edges remaining', color='coral', fontsize=18)
    # ax2.hlines(y=[0.5, 20], xmin=0, xmax=max(cutoff), colors='purple', linestyles='--', lw=2)
    plt.title(f'{batch_type} result', fontsize = 18)
    plt.savefig('{}/Figures/{}_dist_cutoff_vs_numgenes_perc_edges_{}.pdf'.format(savedir, batch_type, date.today().strftime("%m%d%Y")), format='pdf', transparent=True)
    plt.show()
   
######Check how number of proteins changes and the min weight changes at different percentage of edges cutoff#####

## function to do titration at different percent cutoff 
def perc_cut_tit(sorted_df, perc_cutoff):
    sorted_df.columns = ['geneA', 'geneB', 'weight']
    total_edges = len(sorted_df)
    num_genes = []
    dist_max = []
    for i in perc_cutoff:
        keep =math.ceil((i/100)*total_edges)
        # print(keep) 
        df_cutoff = sorted_df.iloc[0:keep]
        gene_list = set(df_cutoff['geneA']).union(set(df_cutoff['geneB']))
        num_genes.append(len(gene_list))
        dist_max.append(min(df_cutoff['weight']))
        # perc_edges.append((num_edges_keep/total_edges)*100)
    df = pd.DataFrame({'Cutoffs': perc_cutoff, 'Number of proteins': num_genes, 'Min weight retained': dist_max})
    return df

## function to plot 
def plot_gene_weight_loss(df, batch_type, savedir):
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax2 = ax1.twinx()
    ax1.plot('Cutoffs', 'Number of proteins', data = df, marker='o',  color = 'navy', linewidth=3, markersize=10)
    ax2.plot('Cutoffs', 'Min weight retained', data = df , marker='o',  color="coral", linewidth=3, markersize=10)

    ax1.set_xlabel('Percent of edges Cutoffs', fontsize=18)
    ax1.set_ylabel('Number of proteins remaining', color='navy', fontsize=18)
    ax2.set_ylabel('Minimum weight remaining', color='coral', fontsize=18)
    plt.title(f'{batch_type} result', fontsize = 18)
    plt.savefig('{}/Figures/{}_sim_edge_perccutoff_vs_numgenes_weight_{}.pdf'.format(savedir, batch_type, date.today().strftime("%m%d%Y")), format='pdf', transparent=True)
    plt.show()
