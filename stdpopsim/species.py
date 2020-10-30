"""
Infrastructure for defining basic information about species and
organising the species catalog.
"""
import logging
import warnings

import attr
import msprime

import stdpopsim

logger = logging.getLogger(__name__)

registered_species = {}


def register_species(species):
    """
    Registers the specified ``species``.

    :param species: The species to be registered.
    :type: :class:`Species`
    """
    if species.id in registered_species:
        raise ValueError(f"{species.id} already registered.")
    logger.debug(f"Registering species '{species.id}'")
    registered_species[species.id] = species


def get_species(id):
    """
    Returns a :class:`Species` object for the specified ``id``.

    :param str id: The string identifier for the requested species. E.g. "HomSap".
         A complete list of species, and their IDs, can be found in the
         :ref:`sec_catalog`.
    :return: An object containing the species definition.
    :rtype: :class:`Species`
    """
    if id not in registered_species:
        # TODO we should probably have a custom exception here and standardise
        # on using these for all the catalog search functions.
        raise ValueError(f"Species '{id}' not in catalog")
    return registered_species[id]


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


def all_demographic_models():
    for species in all_species():
        for model in species.demographic_models:
            yield model


def all_annotations():
    for species in all_species():
        for an in species.annotations:
            yield an


@attr.s()
class Species:
    """
    Class representing a species in the catalog.

    :ivar ~.id: The unique identifier for this species. The species ID is
        the three first letters of the genus name followed by the first
        three letters of the species name, and does not
        contain any spaces or punctuation. The usual scheme is to
        use the first three letters of the genus and species (similar to the
        approach used in the UCSC genome browser), e.g., "HomSap"
        is the ID for Homo Sapiens.
    :vartype ~.id: str
    :ivar name: The full name of this species in binominal nomenclature as
        it would be used in written text, e.g., "Homo sapiens".
    :vartype name: str
    :ivar common_name: The name of this species as it would most often be
        used informally in written text, e.g., "human", or "Orang-utan".
        Where no common name for the species exist, use the most common
        abbreviation, e.g., "E. Coli".
    :vartype common_name: str
    :ivar genome: The :class:`.Genome` instance describing the details
        of this species' genome.
    :vartype genome: stdpopsim.Genome
    :ivar generation_time: The current best estimate for the generation
        time of this species in years. Note that individual demographic
        models in the catalog may or may not use this estimate: each
        model uses the generation time that was used in the original
        publication(s).
    :vartype generation_time: float
    :ivar generation_time_citations: A list of :class:`.Citation` objects
        providing justification for the genertion time estimate.
    :vartype generation_time_citations: list
    :ivar population_size: The current best estimate for the population
        size of this species. Note that individual demographic
        models in the catalog may or may not use this estimate: each
        model uses the population sizes defined in the original
        publication(s).
    :vartype population_size: float
    :ivar population_size_citations: A list of :class:`.Citation` objects
        providing justification for the population size estimate.
    :vartype population_size_citations: list
    :ivar demographic_models: This list of :class:`DemographicModel`
        instances in the catalog for this species.
    :vartype demographic_models: list
    :ivar ensembl_id: The ensembl id for the species' genome assembly,
        which will be used by maintenance scripts to query ensembl's database.
        This parameter will be automatically populated from the species name,
        and should not be set directly unless a non-default assembly is used
        for the species definition (e.g. see E. coli).
    :vartype ensembl_id: str
    """

    id = attr.ib(type=str, kw_only=True)
    name = attr.ib(type=str, kw_only=True)
    common_name = attr.ib(type=str, kw_only=True)
    genome = attr.ib(type=int, kw_only=True)
    generation_time = attr.ib(default=1, kw_only=True)
    generation_time_citations = attr.ib(factory=list, kw_only=True)
    population_size = attr.ib(default=1, kw_only=True)
    population_size_citations = attr.ib(factory=list, kw_only=True)
    demographic_models = attr.ib(factory=list, kw_only=True)
    ensembl_id = attr.ib(type=str, kw_only=True)
    # A list of genetic maps. This is undocumented as the parameter is not
    # intended to be used when the Species is initialsed.
    # Use add_genetic_map() instead.
    genetic_maps = attr.ib(factory=list, kw_only=True)
    annotations = attr.ib(factory=list, kw_only=True)

    @ensembl_id.default
    def _default_ensembl_id(self):
        """
        Returns the ID of this species for the Ensembl REST API.
        This is the species name, underscore delimited and in lowercase.
        """
        return self.name.lower().replace(" ", "_")

    def get_contig(
        self, chromosome, genetic_map=None, length_multiplier=1, sequence_length=None
    ):
        """
        Returns a :class:`.Contig` instance describing a section of genome that
        is to be simulated based on empirical information for a given species
        and chromosome.

        :param str chromosome: The ID of the chromosome to simulate.
             A complete list of chromosome IDs for each species can be found in the
             "Genome" subsection for the species in the :ref:`sec_catalog`.
        :param str genetic_map: If specified, obtain recombination rate information
            from the genetic map with the specified ID. If None, simulate
            using a default uniform recombination rate on a region with the length of
            the specified chromosome. The default rates are species- and chromosome-
            specific, and can be found in the :ref:`sec_catalog`. (Default: None)
        :param float length_multiplier: If specified, simulate a region of length
            `length_multiplier` times the length of the specified chromosome with the
            same chromosome-specific mutation and recombination rates.
            This option cannot currently be used in conjunction with the
            ``genetic_map`` argument.
        :param float sequence_length: Used with a "generic" contig, specifies the
            length of genome sequence for this contig. For a generic contig, mutation
            and recombination rates are equal to the genome-wide average across all
            autosomal chromosomes.
        :rtype: :class:`.Contig`
        :return: A :class:`.Contig` describing the section of the genome.
        """
        # TODO: add non-autosomal support
        if chromosome is not None and chromosome.lower() in (
            "x",
            "y",
            "m",
            "mt",
            "chrx",
            "chry",
            "chrm",
        ):
            warnings.warn(
                stdpopsim.NonAutosomalWarning(
                    "Non-autosomal simulations are not yet supported. See "
                    "https://github.com/popsim-consortium/stdpopsim/issues/383 and "
                    "https://github.com/popsim-consortium/stdpopsim/issues/406"
                )
            )
        if chromosome == "generic":
            if genetic_map is not None:
                raise ValueError("Cannot use genetic map with generic contic")
            if length_multiplier != 1:
                raise ValueError("Cannot use length multiplier for generic contig")
            if sequence_length is None:
                raise ValueError("Must specify sequence_length of generic contig")
            non_autosomal_lower = ["x", "y", "m", "mt", "chrx", "chry", "chrm"]
            L_tot = 0
            r_tot = 0
            u_tot = 0
            for chrom_data in self.genome.chromosomes:
                if chrom_data.id.lower() not in non_autosomal_lower:
                    L_tot += chrom_data.length
                    r_tot += chrom_data.length * chrom_data.recombination_rate
                    u_tot += chrom_data.length * chrom_data.mutation_rate
            u = u_tot / L_tot
            r = r_tot / L_tot
            recomb_map = msprime.RecombinationMap.uniform_map(sequence_length, r)
            ret = stdpopsim.Contig(recombination_map=recomb_map, mutation_rate=u)
        else:
            if sequence_length is not None:
                raise ValueError("Can only use sequence_length with generic contig")
            chrom = self.genome.get_chromosome(chromosome)
            if genetic_map is None:
                logger.debug(f"Making flat chromosome {length_multiplier} * {chrom.id}")
                gm = None
                recomb_map = msprime.RecombinationMap.uniform_map(
                    chrom.length * length_multiplier, chrom.recombination_rate
                )
            else:
                if length_multiplier != 1:
                    raise ValueError("Cannot use length multiplier with empirical maps")
                logger.debug(f"Getting map for {chrom.id} from {genetic_map}")
                gm = self.get_genetic_map(genetic_map)
                recomb_map = gm.get_chromosome_map(chrom.id)

            ret = stdpopsim.Contig(
                recombination_map=recomb_map,
                mutation_rate=chrom.mutation_rate,
                genetic_map=gm,
            )
        return ret

    def get_demographic_model(self, id):
        """
        Returns a demographic model with the specified ``id``.

        :param str id: The string identifier for the demographic model.
             A complete list of IDs for each species can be found in the
             "Demographic Models" subsection for the species in the
             :ref:`sec_catalog`.
        :rtype: :class:`DemographicModel`
        :return: A :class:`DemographicModel` that defines the requested model.
        """
        for model in self.demographic_models:
            if model.id == id:
                return model
        raise ValueError(f"DemographicModel '{self.id}/{id}' not in catalog")

    def add_demographic_model(self, model):
        if model.id in [m.id for m in self.demographic_models]:
            raise ValueError(
                f"DemographicModel '{self.id}/{model.id}' already in catalog."
            )
        self.demographic_models.append(model)

    def add_genetic_map(self, genetic_map):
        if genetic_map.id in [gm.id for gm in self.genetic_maps]:
            raise ValueError(
                f"Genetic map '{self.id}/{genetic_map.id}' " "already in catalog."
            )
        genetic_map.species = self
        self.genetic_maps.append(genetic_map)

    def get_genetic_map(self, id):
        """
        Returns a genetic map (AKA. recombination map) with the specified ``id``.

        :param str id: The string identifier for the genetic map.
             A complete list of IDs for each species can be found in the
             "Genetic Maps" subsection for the species in the :ref:`sec_catalog`.
        :rtype: :class:`GeneticMap`
        :return: A :class:`GeneticMap` that defines the frequency of
            recombinations across the genome.
        """
        for gm in self.genetic_maps:
            if gm.id == id:
                return gm
        raise ValueError(f"Genetic map '{self.id}/{id}' not in catalog")

    def add_annotations(self, annotations):
        if annotations.id in [an.id for an in self.annotations]:
            raise ValueError(
                f"Annotations '{self.id}/{annotations.id}' " "already in catalog."
            )
        annotations.species = self
        self.annotations.append(annotations)

    def get_annotations(self, id):
        """
        Returns a set of annotations with the specified ``id``.

        :param str id: The string identifier for the set of annotations
            A complete list of IDs for each species can be found in the
            "Annotations" subsection for the species in the :ref:`sec_catalog`.
        :rtype: :class:`Annotation`
        :return: A :class:`Annotation` that holds genome annotation
            information from Ensembl
        """
        for an in self.annotations:
            if an.id == id:
                return an
        raise ValueError(f"Annotations '{self.id}/{id}' not in catalog")
