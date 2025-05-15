===============================
Cell Maps Hierarchy Evaluation
===============================
The Cell Maps Hierarchy Evaluation is part of the Cell Mapping Toolkit

.. image:: https://img.shields.io/pypi/v/cellmaps_hierarchyeval.svg
        :target: https://pypi.python.org/pypi/cellmaps_hierarchyeval

.. image:: https://app.travis-ci.com/idekerlab/cellmaps_hierarchyeval.svg?branch=main
        :target: https://app.travis-ci.com/idekerlab/cellmaps_hierarchyeval

.. image:: https://readthedocs.org/projects/cellmaps-hierarchyeval/badge/?version=latest
        :target: https://cellmaps-hierarchyeval.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




evaluate hidef hierarchy


* Free software: MIT license
* Documentation: https://cellmaps-hierarchyeval.readthedocs.io.

Dependencies
------------

* `cellmaps_utils <https://pypi.org/project/cellmaps-utils>`__
* `statsmodels <https://pypi.org/project/statsmodels>`__
* `tqdm <https://pypi.org/project/tqdm>`__
* `scipy <https://pypi.org/project/scipy>`__
* `numpy <https://pypi.org/project/numpy>`__
* `ndex2 <https://pypi.org/project/ndex2>`__

Compatibility
-------------

* Python 3.8+

Installation
------------

.. code-block::

   git clone https://github.com/idekerlab/cellmaps_hierarchyeval
   cd cellmaps_hierarchyeval
   pip install -r requirements_dev.txt
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

Before running tests, please install: ``pip install -r requirements_dev.txt``.

For developers
-------------------------------------------

To deploy development versions of this package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below are steps to make changes to this code base, deploy, and then run
against those changes.

#. Make changes

Modify code in this repo as desired

#. Build and deploy

.. code-block::

    # From base directory of this repo cellmaps_hierarchyeval
    pip uninstall cellmaps_hierarchyeval -y ; make clean dist; pip install dist/cellmaps_hierarchy*whl



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

   Coming soon...

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _NDEx: http://www.ndexbio.org
