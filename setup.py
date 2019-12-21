import os
import codecs
from setuptools import setup

here = os.path.abspath(os.path.dirname(__file__))
with codecs.open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

msprime_ver = "msprime>=0.7.1"

setup(
    name='stdpopsim',
    description='A library of population genetic simulation models',
    long_description=long_description,
    author='PopSim Consortium',
    # TODO probably should have a different email address?
    author_email='popgen_benchmark@lists.uoregon.edu',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    keywords='simulations, recombination map, models',
    packages=['stdpopsim'],
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'stdpopsim=stdpopsim.cli:stdpopsim_main',
        ]
    },
    # NOTE: make sure this is the 'attrs' package, not 'attr'!
    install_requires=[msprime_ver, "attrs", "appdirs", "humanize"],
    url='https://github.com/popsim-consortium/stdpopsim',
    project_urls={
        'Bug Reports': 'https://github.com/popsim-consortium/stdpopsim/issues',
        'Source': 'https://github.com/popsim-consortium/stdpopsim',
    },
    license="GNU GPLv3+",
    platforms=["POSIX", "MacOS X", "Windows"],
    setup_requires=['setuptools_scm'],
    use_scm_version={"write_to": "stdpopsim/_version.py"},
)
