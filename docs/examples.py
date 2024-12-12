"""
The examples used in the tutorial section.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(".."))

import stdpopsim  # noqa: E402


def generic_models_example():
    species = stdpopsim.get_species("HomSap")
    contig = species.get_contig("chr22", length_multiplier=0.1)
    model = stdpopsim.PiecewiseConstantSize(species.population_size)
    samples = {"pop_0": 5}
    engine = stdpopsim.get_default_engine()
    ts = engine.simulate(model, contig, samples)
    print("num_trees =", ts.num_trees)
    print("num_sites =", ts.num_sites)


generic_models_example()
