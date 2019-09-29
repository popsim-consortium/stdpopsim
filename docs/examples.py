"""
The examples used in the tutorial section.
"""
import os
import sys
sys.path.insert(0, os.path.abspath('..'))

import stdpopsim


def generic_models_example():
    species = stdpopsim.get_species("homsap")
    contig = species.get_contig("chr22", length_multiplier=0.1)
    model = stdpopsim.PiecewiseConstantSize(species.population_size)
    samples = model.get_samples(10)
    ts = model.run(contig, samples)



generic_models_example()
