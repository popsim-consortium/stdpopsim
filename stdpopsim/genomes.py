"""
Infrastructure for defining information about species' genomes and genomic
regions to be simulated.
"""
import logging
import warnings
import attr

import numpy as np
import msprime

import stdpopsim

logger = logging.getLogger(__name__)


@attr.s
class Genome:
    """
    Class representing the genome for a species.

    :ivar chromosomes: A list of :class:`.Chromosome` objects.
    :vartype chromosomes: list
    :ivar citations: A list of :class:`.Citation` objects
        providing the source for the genome assembly,
        mutation rate and recombination rate estimates.
    :vartype citations: list
    :ivar length: The total length of the genome.
    :vartype length: int
    """

    # TODO document the assembly_name and accession

    chromosomes = attr.ib(factory=list)
    assembly_name = attr.ib(default=None, kw_only=True)
    assembly_accession = attr.ib(default=None, kw_only=True)
    citations = attr.ib(factory=list, kw_only=True)

    @staticmethod
    def from_data(genome_data, *, recombination_rate, mutation_rate, citations):
        """
        Construct a Genome object from the specified dictionary of
        genome information from Ensembl, recombination_rate and
        mutation_rate dictionaries.

        This method is for internal use only.
        """
        chr_names = set(genome_data["chromosomes"].keys())
        assert set(recombination_rate.keys()) == chr_names
        assert set(mutation_rate.keys()) == chr_names
        chromosomes = []
        for name, data in genome_data["chromosomes"].items():
            chromosomes.append(
                Chromosome(
                    id=name,
                    length=data["length"],
                    synonyms=data["synonyms"],
                    mutation_rate=mutation_rate[name],
                    recombination_rate=recombination_rate[name],
                )
            )
        return Genome(
            chromosomes=chromosomes,
            assembly_name=genome_data["assembly_name"],
            assembly_accession=genome_data["assembly_accession"],
            citations=citations,
        )

    @property
    def length(self):
        return sum(chrom.length for chrom in self.chromosomes)

    def __str__(self):
        s = "Chromosomes:\n"
        length_sorted = sorted(self.chromosomes, key=lambda x: -x.length)
        for chrom in length_sorted:
            s += "\t{}\n".format(chrom)
        return s

    def get_chromosome(self, id):
        """
        Returns the chromosome with the specified ``id``.

        :param str id: The string ID of the chromosome.
            A complete list of chromosome IDs for each species can be found in the
            "Genome" subsection for the species in the :ref:`sec_catalog`.
        :rtype: :class:`Chromosome`
        :return: A :class:`Chromosome` that defines properties of the chromosome
            such as length, mutation rate, and recombination rate.
        """
        for chrom in self.chromosomes:
            if chrom.id == id or id in chrom.synonyms:
                return chrom
        raise ValueError("Chromosome not found")

    @property
    def mean_recombination_rate(self):
        """
        The length-weighted mean recombination rate across all chromosomes.
        """
        length = self.length
        mean_recombination_rate = 0
        for chrom in self.chromosomes:
            normalized_weight = chrom.length / length
            cont = chrom.recombination_rate * normalized_weight
            mean_recombination_rate += cont
        return mean_recombination_rate

    @property
    def mean_mutation_rate(self):
        """
        The length-weighted mean mutation rate across all chromosomes.
        """
        length = self.length
        mean_mutation_rate = 0
        for chrom in self.chromosomes:
            normalized_weight = chrom.length / length
            cont = chrom.mutation_rate * normalized_weight
            mean_mutation_rate += cont
        return mean_mutation_rate


@attr.s
class Chromosome:
    """
    Class representing a single chromosome for a species.

    :ivar str ~.id: The string identifier for the chromosome.
    :ivar int length: The length of the chromosome.
    :ivar float mutation_rate: The mutation rate used when simulating this
        chromosome.
    :ivar float recombination_rate: The recombination rate used when simulating
        this chromosome (if not using a genetic map).
    :ivar synonyms: List of synonyms that may be used when requesting this
        chromosome by ID, e.g. from the command line interface.
    :vartype synonyms: list of str
    """

    id = attr.ib(type=str, kw_only=True)
    length = attr.ib(kw_only=True)
    recombination_rate = attr.ib(type=float, kw_only=True)
    mutation_rate = attr.ib(type=float, kw_only=True)
    synonyms = attr.ib(factory=list, kw_only=True)


@attr.s(kw_only=True)
class GenomicElementType(object):
    """
    Class that represents a map between genomic intervals and a group of
    mutation types and the rate in which they occur within those intervals.

    :ivar mutation_type_ids: Indices for the mutation types
        (:class:`.ext.MutationType`) that can occur in genomic elements of this
        type.
    :vartype mutation_type_ids: list
    :ivar proportions: Proportions (out of the mutation rate) in which each
        mutation type can happen in genomic intervals of this type.
    :vartype: proportions: numpy.array, dtype=np.float32
    :ivar intervals: (n, 2) numpy array with [left, right) genomic intervals
        that will be of this GenomicElementType. The intervals are used to
        specify SLiM GenomicElements. Intervals must be non-overlapping.
    :vartype intervals: numpy.array, dtype=np.int32
    """

    mutation_type_ids = attr.ib(factory=lambda: [], type=list)
    proportions = attr.ib(factory=lambda: np.array([]), type=np.ndarray)
    intervals = attr.ib(default=None, type=np.array)

    def __attrs_post_init__(self):
        _check_mut_proportions(self.proportions, self.mutation_type_ids)
        stdpopsim.utils.check_intervals_validity(self.intervals)


@attr.s(kw_only=True)
class Contig:
    """
    Class representing a contiguous region of genome that is to be
    simulated. This contains the information about mutation rates
    and recombination rates that are needed to simulate this region.
    Genomic element types can be provided.

    :ivar mutation_rate: The rate of mutation per base per generation.
    :vartype mutation_rate: float
    :ivar recombination_map: The recombination map for the region. See the
        :class:`msprime.RateMap` for details.
    :vartype recombination_map: msprime.RateMap
    :ivar mask_intervals: Intervals to keep in simulated tree sequence, as a list
        of (left_position, right_position), such that intervals are non-overlapping
        and in ascending order. Should have shape Nx2, where N is the number of
        intervals.
    :vartype mask_intervals: array-like (?)
    :ivar exclude: If True, ``mask_intervals`` specify regions to exclude. If False,
        ``mask_intervals`` specify regions in keep.
    :vartype exclude: bool
    :ivar genomic_element_types: a list of :class:`.GenomicElementType` objects.
        By default, the only element is a single genomic element type that spans
        the entire contiguous region.
    :vartype genomic_element_types: list
    :ivar mutation_types: a list of :class:`.ext.MutationType` object. By default,
        the only mutation type is neutral. See :class:`.ext.MutationType` for
        more details.

    .. note::
        To run stdpopsim simulations with alternative, user-specified mutation
        or recombination rates, a new contig can be created based on an existing
        one. For instance, the following will create a ``new_contig`` that,
        when simulated, will have double the mutation rate of the ``old_contig``:

        .. code-block:: python

            new_contig = stdpopsim.Contig(
                mutation_rate=old_contig.mutation_rate * 2,
                recombination_map=old_contig.recombination_map,
                genetic_map=old_contig.genetic_map,
            )
    """

    recombination_map = attr.ib()
    mutation_rate = attr.ib(type=float)
    genetic_map = attr.ib(default=None)
    inclusion_mask = attr.ib(default=None)
    exclusion_mask = attr.ib(default=None)
    genomic_element_types = attr.ib(factory=list)
    mutation_types = attr.ib(factory=list)

    def __attrs_post_init__(self):
        self.fully_neutral()

    @staticmethod
    def basic_contig(*, length, mutation_rate=0, recombination_rate=0):
        recomb_map = msprime.RateMap.uniform(length, recombination_rate)
        return Contig(recombination_map=recomb_map, mutation_rate=mutation_rate)

    @staticmethod
    def species_contig(
        *,
        species,
        chromosome=None,
        genetic_map=None,
        length_multiplier=1,
        length=None,
        mutation_rate=None,
        inclusion_mask=None,
        exclusion_mask=None,
    ):
        """
        Build a Contig for a species.
        """
        # TODO: add non-autosomal support
        non_autosomal_lower = ["x", "y", "m", "mt", "chrx", "chry", "chrm"]
        if chromosome is not None and chromosome.lower() in non_autosomal_lower:
            warnings.warn(
                stdpopsim.NonAutosomalWarning(
                    "Non-autosomal simulations are not yet supported. See "
                    "https://github.com/popsim-consortium/stdpopsim/issues/383 and "
                    "https://github.com/popsim-consortium/stdpopsim/issues/406"
                )
            )
        if chromosome is None:
            if genetic_map is not None:
                raise ValueError("Cannot use genetic map with generic contic")
            if length_multiplier != 1:
                raise ValueError("Cannot use length multiplier for generic contig")
            if inclusion_mask is not None or exclusion_mask is not None:
                raise ValueError("Cannot use mask with generic contig")
            if length is None:
                raise ValueError("Must specify sequence length of generic contig")
            L_tot = 0
            r_tot = 0
            u_tot = 0
            for chrom_data in species.genome.chromosomes:
                if chrom_data.id.lower() not in non_autosomal_lower:
                    L_tot += chrom_data.length
                    r_tot += chrom_data.length * chrom_data.recombination_rate
                    u_tot += chrom_data.length * chrom_data.mutation_rate
            if mutation_rate is None:
                mutation_rate = u_tot / L_tot
            r = r_tot / L_tot
            contig = Contig.basic_contig(
                length=length,
                mutation_rate=mutation_rate,
                recombination_rate=r,
            )
        else:
            if length is not None:
                raise ValueError("Cannot specify sequence length for named contig")
            if inclusion_mask is not None and exclusion_mask is not None:
                raise ValueError("Cannot specify both inclusion and exclusion masks")
            chrom = species.genome.get_chromosome(chromosome)
            if genetic_map is None:
                logger.debug(f"Making flat chromosome {length_multiplier} * {chrom.id}")
                gm = None
                recomb_map = msprime.RateMap.uniform(
                    round(chrom.length * length_multiplier), chrom.recombination_rate
                )
            else:
                if length_multiplier != 1:
                    raise ValueError("Cannot use length multiplier with empirical maps")
                logger.debug(f"Getting map for {chrom.id} from {genetic_map}")
                gm = species.get_genetic_map(genetic_map)
                recomb_map = gm.get_chromosome_map(chrom.id)

            inclusion_intervals = None
            exclusion_intervals = None
            if inclusion_mask is not None:
                if length_multiplier != 1:
                    raise ValueError("Cannot use length multiplier with mask")
                if isinstance(inclusion_mask, str):
                    inclusion_intervals = stdpopsim.utils.read_bed(
                        inclusion_mask, chromosome
                    )
                else:
                    inclusion_intervals = inclusion_mask
            if exclusion_mask is not None:
                if length_multiplier != 1:
                    raise ValueError("Cannot use length multiplier with mask")
                if isinstance(exclusion_mask, str):
                    exclusion_intervals = stdpopsim.utils.read_bed(
                        exclusion_mask, chromosome
                    )
                else:
                    exclusion_intervals = exclusion_mask

            if mutation_rate is None:
                mutation_rate = chrom.mutation_rate

            contig = stdpopsim.Contig(
                recombination_map=recomb_map,
                mutation_rate=mutation_rate,
                genetic_map=gm,
                inclusion_mask=inclusion_intervals,
                exclusion_mask=exclusion_intervals,
            )

        return contig

    @property
    def length(self):
        return self.recombination_map.sequence_length

    @property
    def all_intervals_array(self):
        """
        Returns an (n, 3) numpy array with the intervals across all genomic
        element types in the form (left, right, type).
        """
        all_intervals = np.concatenate(
            [
                np.column_stack(
                    (g.intervals[:, :2], np.full((g.intervals.shape[0]), i))
                )
                for i, g in enumerate(self.genomic_element_types)
            ],
            axis=0,
        )
        all_intervals = stdpopsim.utils.build_intervals_array(all_intervals)
        return all_intervals

    def add_mutation_types(self, mutation_types, proportions, genomic_element_type_id):
        """
        Adds mutation types with their respective proportions to the genomic
        element type provided.
        """
        if not (0 <= genomic_element_type_id < len(self.genomic_element_types)):
            raise ValueError(
                f"Genomic element type of id {genomic_element_type_id} does not exist."
            )

        if len(mutation_types) != len(proportions):
            raise ValueError(
                f"number of mutation_types ({len(mutation_types)}) doesn't match "
                f"the number of proportions ({len(proportions)})"
            )

        def _is_index(y, X):
            for i, e in enumerate(X):
                if e is y:
                    return i
            return -1

        for mut_type, prop in zip(mutation_types, proportions):
            mid = _is_index(mut_type, self.mutation_types)
            if mid == -1:
                self.mutation_types.append(mut_type)
                mid = len(self.mutation_types) - 1
            self.genomic_element_types[
                genomic_element_type_id
            ].mutation_type_ids.append(mid)
            self.genomic_element_types[genomic_element_type_id].proportions = np.append(
                self.genomic_element_types[genomic_element_type_id].proportions, [prop]
            )
        genomic_element_type = self.genomic_element_types[genomic_element_type_id]
        _check_mut_proportions(
            genomic_element_type.proportions, genomic_element_type.mutation_type_ids
        )

    def add_genomic_element_type(self, intervals, mutation_types, proportions):
        stdpopsim.utils.check_intervals_validity(intervals)
        ge_type = stdpopsim.GenomicElementType(intervals=intervals)
        self.genomic_element_types.append(ge_type)
        getid = len(self.genomic_element_types) - 1
        self.add_mutation_types(
            mutation_types, proportions, genomic_element_type_id=getid
        )

    def clear_genomic_mutation_types(self):
        self.genomic_element_types = []
        self.mutation_types = []

    def fully_neutral(self, slim_mutations=False, convert_to_substitution=True):
        self.clear_genomic_mutation_types()
        self.add_genomic_element_type(
            intervals=np.array([[0, int(self.length)]]),
            mutation_types=[
                stdpopsim.ext.MutationType(
                    convert_to_substitution=convert_to_substitution
                )
            ],
            proportions=[1 if slim_mutations else 0],
        )

    def __str__(self):
        gmap = "None" if self.genetic_map is None else self.genetic_map.id
        s = (
            "Contig(length={:.2G}, recombination_rate={:.2G}, "
            "mutation_rate={:.2G}, genetic_map={}, genomic_element_types={}, "
            "mutation_types={})"
        ).format(
            self.length,
            self.recombination_map.mean_rate,
            self.mutation_rate,
            gmap,
            self.genomic_element_types,
            self.mutation_types,
        )
        return s


def _check_mut_proportions(proportions, mutation_type=None, atol=1e-15):
    if mutation_type is not None and len(proportions) != len(mutation_type):
        raise ValueError(
            "The list of proportions and mutation types should be of the same length."
        )
    if any([not (0 <= prop <= 1) for prop in proportions]):
        raise ValueError("Proportions must lie within the [0,1] interval.")
    if sum(proportions) > 1 + atol:
        raise ValueError(
            "The sum of the proportions within genomic element types cannot"
            " exceed 1."
        )
