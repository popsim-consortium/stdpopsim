"""
Infrastructure for defining basic information about the genomes of
species.
"""

import stdpopsim.genetic_maps as genetic_maps


class Genome(object):
    """
    Class representing the genome for a species.
    """
    def __init__(self, species, chromosomes, default_genetic_map=None):
        self.species = species
        self.default_genetic_map = default_genetic_map
        self.chromosomes = {}
        for chromosome in chromosomes:
            self.chromosomes[chromosome.name] = chromosome
            chromosome.default_genetic_map = default_genetic_map
            chromosome.species = species

    def __str__(self):
        s = "Genome for {}:\n".format(self.species)
        s += "Chromosomes:\n"
        length_sorted = sorted(self.chromosomes.values(), key=lambda x: -x.length)
        for chrom in length_sorted:
            s += "\t{}\n".format(chrom)
        return s


class Chromosome(object):
    """
    Class representing a single chromosome for a species.
    """

    def __init__(self, name, length, mean_recombination_rate, mean_mutation_rate):
        self.name = name
        self.length = length
        self.mean_recombination_rate = mean_recombination_rate
        self.mean_mutation_rate = mean_mutation_rate
        self.species = None
        self.default_genetic_map = None

    def __repr__(self):
        return (
            "{{'name': {}, 'length': {}, "
            "'mean_recombination_rate': {}, "
            "'mean_mutation_rate': {}}}".format(
                self.name, self.length, self.mean_recombination_rate,
                self.mean_mutation_rate))

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
        return genetic_map.get_chromosome_map(self.name)
