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

**Optional**

Logging and verbosity options.

Via Docker
---------------

**Example usage**

**TODO:** Add information about example usage


.. code-block::

   Coming soon...

