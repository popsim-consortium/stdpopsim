#!/usr/bin/env python3
from setuptools import setup

setup(
    # Set the name so that github correctly tracks reverse dependencies.
    # https://github.com/popsim-consortium/stdpopsim/network/dependents
    name="stdpopsim",
    use_scm_version={"write_to": "stdpopsim/_version.py"},
)
