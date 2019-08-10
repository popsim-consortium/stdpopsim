"""
Infrastructure for defining basic information about the genomes of
species.
"""
import warnings
import logging

import msprime

import stdpopsim.genetic_maps as genetic_maps

logger = logging.getLogger(__name__)

registered_genomes = {}


def get_genome(species):
    """
    Returns the genome definition for the specified species.
    Raises a ValueError if the species has not been registered.
    """
    if species not in registered_genomes:
        raise ValueError("Unknown species '{}'".format(species))
    return registered_genomes[species]


def register_genome(genome):
    """
    Registers the genome definition that it can be loaded.
    """
    key = genome.species
    logger.debug("Registering genome '{}'".format(key))
    registered_genomes[key] = genome


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

    @property
    def mean_mutation_rate(self):
        """
        This method return the weighted mean mutation rate
        across all chomosomes in the genome.
        :rtype: float
        """
        mean_mutation_rate = 0
        for chrom in self.chromosomes.values():
            normalized_weight = chrom.length / self.length
            cont = chrom.default_mutation_rate*normalized_weight
            mean_mutation_rate += cont
        return mean_mutation_rate


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


class Contig(object):
    """
    Class representing a contiguous region of genome that is to be
    simulated. This contains all information required to simulate
    the region, including mutation, recombination rates, etc.

    See :func:`.contig_factory` function for information on how to
    populate this
    """
    def __init__(self):
        self.recombination_map = None
        self.mutation_rate = None


def contig_factory(species, chromosome, genetic_map=None, length_multiplier=1):
    """
    Returns a :class:`.Contig` instance describing a section of genome that
    is to be simulated based on empirical information for a given species
    and chromosome.

    :param str species: The name of the species to simulate.
    :param str chromosome: The name of the chromosome to simulate.
    :param str genetic_map: If specified, obtain recombination rate information
        from the genetic map with the specified name. If None, simulate
        a flat recombination rate on a region with the length of the specified
        chromosome. (Default: None)
    :param float length_multiplier: If specified simulate a contig of length
        length_multiplier times the length of the specified chromosome.
    :rtype: Contig
    :return: A Contig describing a simulation of the section of genome.
    """
    genome = get_genome(species)
    chrom = genome.chromosomes[chromosome]
    if genetic_map is None:
        logger.debug(f"Making flat chromosome {length_multiplier} * {chrom.name}")
        recomb_map = msprime.RecombinationMap.uniform_map(
            chrom.length * length_multiplier, chrom.default_recombination_rate)
    else:
        if length_multiplier != 1:
            raise ValueError("Cannot use length multiplier with empirical maps")
        logger.debug(f"Getting map for {chrom.name} from {genetic_map}")
        recomb_map = chrom.recombination_map(genetic_map)

    ret = Contig()
    ret.recombination_map = recomb_map
    ret.mutation_rate = chrom.default_mutation_rate
    return ret
