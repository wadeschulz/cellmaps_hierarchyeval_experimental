# Step by Step Guide for HiDeF hierarchy evaluation

Evaluate hierarchy from HiDeF 

Check the [Notebook](./HiDeF_Hierarchy_eval_pipeline.ipynb) for the cleaned pipeline to run the evaluation analysis 

### Step 0:

Use [hidef_finder.py](https://github.com/fanzheng10/HiDeF/blob/master/hidef/hidef_finder.py) from HiDeF repo to generate the hierarchical structures from protein interaction networks. 

```
python -u <Function PATH>/hidef_finder.py --g < PATH one input network or list of input networks > --k 5 --maxres 40 --alg leiden --o <OUTPATH>
```
where k is the stability, maxres is the maximum resolution, alg is the community detection algorithm. Other parameters in HiDeF can also be tuned check out the paper for details: Zheng, F., Zhang, S., Churas, C. et al., HiDeF: identifying persistent structures in multiscale ‘omics data. Genome Biol 22, 21 (2021).

### Step 1:

**Run analyze hidef enrichment**

```
python -u ./U2OS5183_analyze_hidef_enrichment.py $PARAM
```

PARAM: 

--infname: input path and the prefix of HiDeF nodes and edges 

--outprefix: output path and prefix 

--w_root: default noRoot (do remove the root from the enrichment analysis)

--minTermSize: minimum size of the term to consider when running enrichment (default = 4)

--FDRthre: cutoff of FDR, default is 0.01


### Step 2:

**Organize hidef enrichment into a single table**

Use function 'analyze_enrichment' in the utils file [hidef_enrichment_analysis_utils.py](hidef_enrichment_analysis_utils.py)

Check the [Notebook](./HiDeF_Hierarchy_eval_pipeline.ipynb) for the cleaned pipeline to run the evaluation analysis 

