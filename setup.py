#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""
import os
import re
from setuptools import setup, find_packages


with open(os.path.join('cellmaps_hierarchyeval', '__init__.py')) as ver_file:
    for line in ver_file:
        if line.startswith('__version__'):
            version=re.sub("'", "", line[line.index("'"):])

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['pandas',
                'networkx',
                'numpy',
                'matplotlib',
                'seaborn',
                'tqdm',
                'statsmodels',
                'cellmaps_utils']

setup_requirements = [ ]

test_requirements = [ ]

setup(
    author="Clara Hu",
    author_email='mhu@health.ucsd.edu',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="evaluate hidef hierarchy ",
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    long_description_content_type = 'text/x-rst',
    include_package_data=True,
    keywords='cellmaps_hierarchyeval',
    name='cellmaps_hierarchyeval',
    packages=find_packages(include=['cellmaps_hierarchyeval']),
    package_dir={'cellmaps_hierarchyeval': 'cellmaps_hierarchyeval'},
    scripts=['cellmaps_hierarchyeval/cellmaps_hierarchyevalcmd.py',
             'cellmaps_hierarchyeval/mature_hier_structure.py',
             'cellmaps_hierarchyeval/U2OS5183_analyze_hidef_enrichment.py'],
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/idekerlab/cellmaps_hierarchyeval',
    version=version,
    zip_safe=False)
