"""
Infrastructure for defining basic information about the genomes of
species.
"""
import warnings
import logging

import msprime

import stdpopsim.genetic_maps as genetic_maps

logger = logging.getLogger(__name__)

registered_species = {}


def register_species(species):
    """
    Registers the specified species.
    """
    logger.debug(f"Registering species '{species.name}'")
    registered_species[species.name] = species


def get_species(name):
    return registered_species[name]


class Species(object):
    """
    Class representing a single species.
    """
    def __init__(
            self, name, genome, generation_time=None, population_size=None):
        self.name = name
        self.genome = genome
        self.generation_time = generation_time
        self.population_size = population_size
        self.models = []
        self.genetic_maps = []

    def get_contig(self, chromosome, genetic_map=None, length_multiplier=1):
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
        chrom = self.genome.get_chromosome(chromosome)
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

    def get_model(self, kind, num_populations=1):
        """
        Returns a model with the specified name with the specified number of
        populations. Please see the documentation [] for a list of available
        models.

        - TODO we can add functionality here
        """
        for model in self.models:
            if model.kind == kind and len(model.populations) == num_populations:
                return model
        raise ValueError("Model not found")

    def add_model(self, model):
        self.models.append(model)

    def add_genetic_map(self, genetic_map):
        self.genetic_maps.append(genetic_map)


class Genome(object):
    """
    Class representing the genome for a species.

    .. todo:: Define the facilities that this object provides.
    """
    def __init__(self, species, chromosomes, default_genetic_map=None):
        self.species = species
        self.default_genetic_map = default_genetic_map
        self.chromosomes = chromosomes
        self.length = 0
        for chromosome in chromosomes:
            self.length += chromosome.length

    def __str__(self):
        s = "Genome for {}:\n".format(self.species)
        s += "Chromosomes:\n"
        length_sorted = sorted(self.chromosomes.values(), key=lambda x: -x.length)
        for chrom in length_sorted:
            s += "\t{}\n".format(chrom)
        return s

    def get_chromosome(self, name):
        """
        Returns the chromosome with the specified name.
        """
        for chrom in self.chromosomes:
            if chrom.name == name:
                return chrom
        raise ValueError("Chromosome not found")

    @property
    def mean_recombination_rate(self):
        """
        This method return the weighted mean recombination rate
        across all chromosomes in the genome.
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
        # TODO remove the default here --- default for what? The chromosome object
        # models the actual chromsomes, it shouldn't need to know about other bits
        # of the software mode. We should either change it to
        # recombination_rate or mean_recombination_rate, if we want to be
        # picky.
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
