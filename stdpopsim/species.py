"""
Infrastructure for defining basic information about species and
organising the species catalog.
"""
import logging

import attr
import msprime

from . import genomes

logger = logging.getLogger(__name__)

registered_species = {}


def register_species(species):
    """
    Registers the specified species.
    """
    logger.debug(f"Registering species '{species.id}'")
    registered_species[species.id] = species


def get_species(id_):
    if id_ not in registered_species:
        # TODO we should probably have a custom exception here and standardise
        # on using these for all the catalog search functions.
        raise ValueError("species not found")
    return registered_species[id_]


# Convenience methods for getting all the species/genetic maps/models
# we have defined in the catalog.

def all_species():
    """
    Returns an iterator over all species in the catalog.
    """
    for species in registered_species.values():
        yield species


def all_genetic_maps():
    for species in all_species():
        for genetic_map in species.genetic_maps:
            yield genetic_map


def all_models():
    for species in all_species():
        for model in species.models:
            yield model


@attr.s(frozen=True)
class Species(object):
    """
    Class representing a species in the catalog.

    :ivar id: The unique identifier for this species. The species ID is
        usually an abbreviation of the species name, which does not
        contain any spaces or puncutation.. The usual scheme is to
        use the first three letters of the genus and species (similar to the
        approach used in the UCSC genome browser), e.g., "homsap"
        is the ID for Homo Sapiens.
    :vartype id: str
    :ivar name: The informal name for this species as it would
        be used in written text, e.g., "Homo sapiens"
    :vartype informal_name: str
    :ivar genome: The :class:`.Genome` instance describing the details
        of this species' genome.
    :vartype genome: stdpopsim.Genome
    :ivar generation_time: The current best estimate for the generation
        time of this species in years. Note that individual population
        models in the catalog may or may not use this estimate: each
        model uses the generation time that was used in the original
        publication(s).
    :vartype generation_time: float
    :ivar generation_time_citations: A list of :class:`.Citation` objects
        providing justification for the genertion time estimate.
    :vartype generation_time_citations: list
    :ivar population_size: The current best estimate for the population
        size of this species. Note that individual population
        models in the catalog may or may not use this estimate: each
        model uses the populations sizes defined in the original
        publication(s).
    :vartype population_size: float
    :ivar population_size_citations: A list of :class:`.Citation` objects
        providing justification for the population size estimate.
    :vartype population_size_citations: list
    :ivar models: This list of :class:`Model` instances in the catalog
        for this species.
    :vartype models: list()
    """

    id = attr.ib(type=str, kw_only=True)
    name = attr.ib(type=str, kw_only=True)
    genome = attr.ib(type=int, kw_only=True)
    generation_time = attr.ib(default=1, kw_only=True)
    generation_time_citations = attr.ib(factory=list, kw_only=True)
    population_size = attr.ib(default=1, kw_only=True)
    population_size_citations = attr.ib(factory=list, kw_only=True)
    models = attr.ib(factory=list, kw_only=True)
    genetic_maps = attr.ib(factory=list, kw_only=True)

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
                chrom.length * length_multiplier, chrom.recombination_rate)
        else:
            if length_multiplier != 1:
                raise ValueError("Cannot use length multiplier with empirical maps")
            logger.debug(f"Getting map for {chrom.name} from {genetic_map}")
            gm = self.get_genetic_map(genetic_map)
            recomb_map = gm.get_chromosome_map(chrom.name)

        ret = genomes.Contig(
            recombination_map=recomb_map, mutation_rate=chrom.mutation_rate)
        return ret

    def get_model(self, id_):
        """
        Returns a model with the specified id.

        - TODO explain where we find models from the catalog.
        """
        for model in self.models:
            if model.id == id_:
                return model
        raise ValueError("Model not found")

    def add_model(self, model):
        self.models.append(model)

    def add_genetic_map(self, genetic_map):
        genetic_map.species = self
        self.genetic_maps.append(genetic_map)

    def get_genetic_map(self, name):
        for gm in self.genetic_maps:
            if gm.name == name:
                return gm
        raise ValueError("Genetic map not found")
