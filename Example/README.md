# What do each example file mean and How to use them 

This folder provide examples of parameter sweeping with HiDeF 

The input to HiDeF is U2OS_5183_MUSE_batch_hard_percthre_0.2.tsv, and sweep the stability (chi) in the range of 5, 10, 15, maxresolution(maxres) is either 20 or 40. 

All the nodes, edges, and gml files are output from HiDeF. The file direction of nodes and edges files are used as input for the function [U2OS5183_analyze_hidef_enrichment.py](../U2OS5183_analyze_hidef_enrichment.py)

The output of this function are the '{file prefix}.noRoot.{date}.pkl' or '{file prefix}.noRoot.{date}.csv', which is the input to the enrichment analysis ([hidef_enrichment_analysis_utils.py](../hidef_enrichment_analysis_utils.py))

NOTE: if use the .pkl file, the pandas version I use is pandas==1.3.5, python3.7

The output of the final enrichment analysis is [U2OS_5183_MUSE_256_w_PCNet_Cutoff_hidef_enrichment_analysis.csv](./U2OS_5183_MUSE_256_w_PCNet_Cutoff_hidef_enrichment_analysis.csv)
