"""
Infrastructure for defining basic information about the genomes of
species.
"""

import stdpopsim.genetic_maps as genetic_maps
import msprime
import warnings


class Genome(object):
    """
    Class representing the genome for a species.

    .. todo:: Define the facilities that this object provides.
    """
    def __init__(self, species, chromosomes, default_genetic_map=None):
        self.species = species
        self.default_genetic_map = default_genetic_map
        self.chromosomes = {}
        self.length = 0
        for chromosome in chromosomes:
            self.chromosomes[chromosome.name] = chromosome
            chromosome.default_genetic_map = default_genetic_map
            chromosome.species = species
            self.length += chromosome.length

    def __str__(self):
        s = "Genome for {}:\n".format(self.species)
        s += "Chromosomes:\n"
        length_sorted = sorted(self.chromosomes.values(), key=lambda x: -x.length)
        for chrom in length_sorted:
            s += "\t{}\n".format(chrom)
        return s

    @property
    def mean_recombination_rate(self):
        """
        This method return the weighted mean recombination rate
        across all chomosomes in the genome.
        :rtype: float
        """
        mean_recombination_rate = 0
        for chrom in self.chromosomes.values():
            normalized_weight = chrom.length / self.length
            cont = chrom.default_recombination_rate*normalized_weight
            mean_recombination_rate += cont
        return mean_recombination_rate


class Chromosome(object):
    """
    Class representing a single chromosome for a species.

    .. todo:: Define the facilities that this object provides.
    """
    def __init__(self, name, length, default_recombination_rate, default_mutation_rate):
        self.name = name
        self.length = length
        self.default_recombination_rate = default_recombination_rate
        self.default_mutation_rate = default_mutation_rate
        self.species = None
        self.default_genetic_map = None

    def __repr__(self):
        return (
            "{{'name': {}, 'length': {}, "
            "'default_recombination_rate': {}, "
            "'default_mutation_rate': {}}}".format(
                self.name, self.length, self.default_recombination_rate,
                self.default_mutation_rate))

    def __str__(self):
        return repr(self)

    def recombination_map(self, map_name=None):
        """
        Returns an :class:`msprime.RecombinationMap` instance representing the
        recombination map for this chromosome. If ``map_name`` is provided,
        return the corresponding recombination map; if not, use the default
        recombination map for this species.
        """
        if map_name is None:
            map_name = self.default_genetic_map
        genetic_map = genetic_maps.get_genetic_map(self.species, map_name)
        if genetic_map.contains_chromosome_map(self.name):
            ret = genetic_map.get_chromosome_map(self.name)
        else:
            warnings.warn(
                "Warning: recombination map not found for chromosome: '{}'"
                " on map: '{}', substituting a zero"
                "-recombination map.".format(self.name, map_name))
            ret = msprime.RecombinationMap.uniform_map(self.length, 0)
        return ret
