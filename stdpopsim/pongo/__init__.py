"""
Pongo models, recombination and mutation rates.
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

# Maybe good to package this part differently if wanting to use
# some information but not all. 
# E.g. If you want mean mutation rate and recomb rates, but user-def length
chrom = Chromosome(
    length=50818468,  # from H.sapiens chr22, Wikipedia
    mean_mutation_rate=2e-8, # assumption from Locke etal
    mean_recombination_rate=.95*1e-8)  # from Locke etal


# ETC
