# Evaluate_Hierarchy

Evaluate hierarchy from HiDeF 

### Step 1:

**Run analyze hidef enrichment**

--infname: input path and the prefix of HiDeF nodes and edges 

--outprefix: output path and prefix 

--w_root: default noRoot (do remove the root from the enrichment analysis

--minTermSize: minimum size of the term to consider when running enrichment (default = 4)

--FDRthre: cutoff of FDR, default is 0.01


### Step 2:

**Organize hidef enrichment into a single table**

Use function 'analyze_enrichment' in the utils file [hidef_enrichment_analysis_utils.py](hidef_enrichment_analysis_utils.py)

Check the [Notebook](./HiDeF_Hierarchy_eval_pipeline.ipynb) for the cleaned pipeline to run the evaluation analysis 

