"""
Human models, recombination and mutation rates.
"""

from . import models


# TODO this infrastructure should live somewhere else in the package
# hierarchy so it can be reused across the different species.

class Chromosome(object):

    def __init__(self, length, mean_recombination_rate, mean_mutation_rate):
        self.length = length
        self.mean_recombination_rate = mean_recombination_rate
        self.mean_mutation_rate = mean_mutation_rate

    # TODO add methods to return recombination maps

    # Add methods to print this out. __str__ should give a nice summary.


# Define the chromosomes.

chr22 = Chromosome(
    length=50818468,  # Taken from wikipedia, but should really be based on GRCh38.
    mean_mutation_rate=1e-8,  # WRONG!
    mean_recombination_rate=1e-8)  # WRONG!


# ETC
