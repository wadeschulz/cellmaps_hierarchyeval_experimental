======================
cellmaps_hierarchyeval
======================


.. image:: https://img.shields.io/pypi/v/cellmaps_hierarchyeval.svg
        :target: https://pypi.python.org/pypi/cellmaps_hierarchyeval

.. image:: https://img.shields.io/travis/idekerlab/cellmaps_hierarchyeval.svg
        :target: https://travis-ci.com/idekerlab/cellmaps_hierarchyeval

.. image:: https://readthedocs.org/projects/cellmaps-hierarchyeval/badge/?version=latest
        :target: https://cellmaps-hierarchyeval.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




evaluate hidef hierarchy 


* Free software: MIT license
* Documentation: https://cellmaps-hierarchyeval.readthedocs.io.



Dependencies
------------

* TODO add

Compatibility
-------------

* Python 3.3+

Installation
------------

.. code-block::

   git clone https://github.com/idekerlab/cellmaps_hierarchyeval
   cd cellmaps_hierarchyeval
   make dist
   pip install dist/cellmaps_hierarchyevalcmd*whl


Run **make** command with no arguments to see other build/deploy options including creation of Docker image 

.. code-block::

   make

Output:

.. code-block::

   clean                remove all build, test, coverage and Python artifacts
   clean-build          remove build artifacts
   clean-pyc            remove Python file artifacts
   clean-test           remove test and coverage artifacts
   lint                 check style with flake8
   test                 run tests quickly with the default Python
   test-all             run tests on every Python version with tox
   coverage             check code coverage quickly with the default Python
   docs                 generate Sphinx HTML documentation, including API docs
   servedocs            compile the docs watching for changes
   testrelease          package and upload a TEST release
   release              package and upload a release
   dist                 builds source and wheel package
   install              install the package to the active Python's site-packages
   dockerbuild          build docker image and store in local repository
   dockerpush           push image to dockerhub




Needed files
------------

**TODO:** Add description of needed files


Usage
-----

For information invoke :code:`cellmaps_hierarchyevalcmd.py -h`

**Example usage**

**TODO:** Add information about example usage

.. code-block::

   cellmaps_hierarchyevalcmd.py # TODO Add other needed arguments here


Via Docker
~~~~~~~~~~~~~~~~~~~~~~

**Example usage**

**TODO:** Add information about example usage


.. code-block::

   docker run -v `pwd`:`pwd` -w `pwd` idekerlab/cellmaps_hierarchyeval:0.1.0 cellmaps_hierarchyevalcmd.py # TODO Add other needed arguments here



Step by Step Guide for HiDeF hierarchy evaluation
--------------------------------------------------
Evaluate hierarchy from HiDeF 

Check the [Notebook](./HiDeF_Hierarchy_eval_pipeline.ipynb) for the cleaned pipeline to run the evaluation analysis 

Step 0:
^^^^^^^^^^^^^^

Use [hidef_finder.py](https://github.com/fanzheng10/HiDeF/blob/master/hidef/hidef_finder.py) from HiDeF repo to generate the hierarchical structures from protein interaction networks. 

.. code-block::

   python -u <Function PATH>/hidef_finder.py --g < PATH one input network or list of input networks > --k 5 --maxres 40 --alg leiden --o <OUTPATH>


where k is the stability, maxres is the maximum resolution, alg is the community detection algorithm. Other parameters in HiDeF can also be tuned check out the paper for details: Zheng, F., Zhang, S., Churas, C. et al., HiDeF: identifying persistent structures in multiscale ‘omics data. Genome Biol 22, 21 (2021).

Step 1
^^^^^^^

**Run analyze hidef enrichment**

.. code-block::

   python -u ./U2OS5183_analyze_hidef_enrichment.py $PARAM

PARAM: 

--infname: input path and the prefix of HiDeF nodes and edges 

--outprefix: output path and prefix 

--w_root: default noRoot (do remove the root from the enrichment analysis)

--minTermSize: minimum size of the term to consider when running enrichment (default = 4)

--FDRthre: cutoff of FDR, default is 0.01


Step 2:
^^^^^^^^

**Organize hidef enrichment into a single table**

Use function 'analyze_enrichment' in the utils file [hidef_enrichment_analysis_utils.py](hidef_enrichment_analysis_utils.py)

Check the [Notebook](./HiDeF_Hierarchy_eval_pipeline.ipynb) for the cleaned pipeline to run the evaluation analysis 



Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _NDEx: http://www.ndexbio.org
