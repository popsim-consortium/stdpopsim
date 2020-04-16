"""
The examples used in the tutorial section.
"""
import os
import sys
sys.path.insert(0, os.path.abspath('..'))

import stdpopsim


def generic_models_example():
    species = stdpopsim.get_species("HomSap")
    contig = species.get_contig("chr22", length_multiplier=0.1)
    model = stdpopsim.PiecewiseConstantSize(species, species.population_size)
    samples = model.get_samples(10)
    engine = stdpopsim.get_default_engine()
    ts = engine.simulate(model, contig, samples)



generic_models_example()
