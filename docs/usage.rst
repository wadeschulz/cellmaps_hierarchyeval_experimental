=====
Usage
=====

The **cellmaps_hierarchyeval** tool is designed for conducting enrichment analyses on hierarchical networks.
Its objective is to show the biological significance of various terms and components (such as genes and proteins)
within these networks. Terms are used from Gene Ontology (GO), CORUM, and the Human Protein Atlas (HPA).

In a project
--------------

To use cellmaps_hierarchyeval in a project::

    import cellmaps_hierarchyeval

On the command line
---------------------

For information invoke :code:`cellmaps_hierarchyevalcmd.py -h`

**Usage**

.. code-block::

    cellmaps_hierarchyevalcmd.py [outdir] [--hierarchy_dir HIERARCHYDIR] [OPTIONS]

**Arguments**

- ``outdir``
    The directory where the output, including enriched hierarchies, will be written.

- ``--hierarchy_dir HIERARCHYDIR``
    Directory containing the generated hierarchy from cellmaps_generate_hierarchy.

- ``--max_fdr``
    Maximum false discovery rate for enrichment analysis. Default is 0.05.

- ``--min_jaccard_index``
    Minimum Jaccard index for considering an enrichment result significant. Default is 0.1.

- ``--min_comp_size``
    Minimum size of a term to be considered for enrichment. Default is 4.

- ``--corum``
    UUID for CORUM network. Default is provided.

- ``--go_cc``
    UUID for GO-CC network. Default is provided.

- ``--hpa``
    UUID for HPA network. Default is provided.

- ``--ndex_server``
    NDEx server to use. Default is http://www.ndexbio.org.

- ``--skip_logging``
    If set, disables the creation of log files.

- ``--skip_term_enrichment``
    If set, SKIP enrichment against networks set via --corum, --go_cc, --hpa

- ``--ollama``
    Path to ollama command line binary or REST service. If value starts with http it is assumed to be a REST url and
    all prompts will be passed to service. For REST url the suffix api/generate must be appended.
    Example: http://foo/api/generate. NOTE: ollama integration with this tool is EXPERIMENTAL and interface may be
    changed or removed in the future.

- ``--ollama_user``
    Username to pass as basic auth to ollama REST service

- ``--ollama_password``
    Password to pass via basic autho to ollama REST service

- ``--ollama_prompts``
    Comma delimited value of format <MODEL NAME> or <MODEL NAME>,<PROMPT> where <PROMPT> can be path to prompt file or
    prompt to run. For insertion of gene set please include {GENE_SET} in prompt and tell LLM to put Process: <name> on
    first line with name assigned to assembly and Confidence Score: <score> on 2nd line with confidence in the name
    given. If just <MODEL NAME> is set, then default prompt is used with model specified. NOTE: if <MODEL NAME> is set
    to FAKE then a completely fake agent will be used. Also note: ollama integration with this tool is EXPERIMENTAL and
    interface may be changed or removed in the future.

- ``--provenance``
    Path to file containing provenance information about input files in JSON format. This is required if inputdir
    does not contain ro-crate-metadata.json file.

**Optional**

Logging and verbosity options.

Via Docker
---------------

**Example usage**

**TODO:** Add information about example usage


.. code-block::

   Coming soon...

