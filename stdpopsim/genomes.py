"""
Infrastructure for defining information about species' genomes and genomic
regions to be simulated.
"""
import logging
import attr
import numpy as np
import stdpopsim

logger = logging.getLogger(__name__)


@attr.s
class Genome:
    """
    Class representing the genome for a species.

    :ivar chromosomes: A list of :class:`.Chromosome` objects.
    :vartype chromosomes: list
    :ivar mutation_rate_citations: A list of :class:`.Citation` objects
        providing justification for the mutation rate estimate.
    :vartype mutation_rate_citations: list
    :ivar recombination_rate_citations: A list of :class:`.Citation` objects
        providing justification for the recombination rate estimate.
    :vartype recombination_rate_citations: list
    :ivar assembly_citations: A list of :class:`.Citation` objects
        providing reference to the source of the genome assembly.
    :vartype assembly_citations: list
    :ivar length: The total length of the genome.
    :vartype length: int
    """

    chromosomes = attr.ib(factory=list)
    mutation_rate_citations = attr.ib(factory=list, kw_only=True)
    recombination_rate_citations = attr.ib(factory=list, kw_only=True)
    assembly_citations = attr.ib(factory=list, kw_only=True)
    assembly_name = attr.ib(default=None, kw_only=True)
    assembly_accession = attr.ib(default=None, kw_only=True)
    length = attr.ib(default=0, init=False)

    def __attrs_post_init__(self):
        for chromosome in self.chromosomes:
            self.length += chromosome.length

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
        mean_recombination_rate = 0
        for chrom in self.chromosomes:
            normalized_weight = chrom.length / self.length
            cont = chrom.recombination_rate * normalized_weight
            mean_recombination_rate += cont
        return mean_recombination_rate

    @property
    def mean_mutation_rate(self):
        """
        The length-weighted mean mutation rate across all chromosomes.
        """
        mean_mutation_rate = 0
        for chrom in self.chromosomes:
            normalized_weight = chrom.length / self.length
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

    :ivar mutation_type_ids: Indices for the mutation types that can occur in
        genomic elements of this type.
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


@attr.s
class Contig:
    """
    Class representing a contiguous region of genome that is to be
    simulated. This contains the information about mutation rates
    and recombination rates that are needed to simulate this region.
    Genomic element types can be provided

    :ivar mutation_rate: The rate of mutation per base per generation.
    :vartype mutation_rate: float
    :ivar recombination_map: The recombination map for the region. See the
        `msprime documentation
        <https://msprime.readthedocs.io/en/stable/api.html#msprime.RecombinationMap>`_
        for more details.
    :vartype recombination_map: msprime.simulations.RecombinationMap
    :ivar mask_intervals: Intervals to keep in simulated tree sequence, as a list
        of (left_position, right_position), such that intervals are non-overlapping
        and in ascending order. Should have shape Nx2, where N is the number of
        intervals.
    :vartype mask_intervals: array-like (?)
    :ivar exclude: If True, ``mask_intervals`` specify regions to exclude. If False,
        ``mask_intervals`` specify regions in keep.
    :vartype exclude: bool
    :ivar genomic_element_types: a list of stdpopsim.GenomicElementType
        objects. By default, the only element is a single genomic element type
        that spans the entire contiguous region.
    :vartype genomic_element_types: list
    :ivar mutation_types: a list of stdpopsim.ext.MutationType object. By
        default, the only mutation type is neutral. See
        stdpopsim.MutationType for more details.

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

    recombination_map = attr.ib(default=None, kw_only=True)
    mutation_rate = attr.ib(default=None, type=float, kw_only=True)
    genetic_map = attr.ib(default=None, kw_only=True)
    inclusion_mask = attr.ib(default=None, kw_only=True)
    exclusion_mask = attr.ib(default=None, kw_only=True)
    genomic_element_types = attr.ib(factory=lambda: [], type=list, kw_only=True)
    mutation_types = attr.ib(factory=lambda: [], type=list, kw_only=True)

    @property
    def length(self):
        if self.recombination_map is None:
            return None
        return self.recombination_map.get_length()

    @property
    def slim_fractions(self):
        map(lambda g: _check_mut_proportions(g.proportions), self.genomic_element_types)
        props = np.array([sum(g.proportions) for g in self.genomic_element_types])
        return props

    @property
    def genomic_intervals(self):
        """
        Returns an (n, 3) numpy array with n intervals (left, right, type)
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

    @property
    def msp_mutation_rate_map(self):
        """
        breaks = [0,10,20] means [0,10) and [10,20)
        """
        breaks = [0]
        rates = []
        for (start, end, t) in self.genomic_intervals:
            if start not in breaks:
                breaks.append(start)
                rates.append(self.mutation_rate)
            breaks.append(end)
            rates.append((1 - self.slim_fractions[t]) * self.mutation_rate)
        if not np.isclose(breaks[-1], int(self.length)):
            breaks.append(int(self.length))
            rates.append(self.mutation_rate)
        assert len(breaks) == len(rates) + 1
        return (breaks, rates)

    @property
    def slim_mutation_rate_map(self):
        """
        breaks = [10,20] means [0,10] and (10,20]
        """
        # making sure the mut types are still concordant with the ge types
        breaks, rates = self.msp_mutation_rate_map
        assert breaks[0] == 0
        slim_breaks = [b - 1 for b in breaks[1:]]
        slim_rates = [self.mutation_rate - r for r in rates]
        return (slim_breaks, slim_rates)

    def add_mutation_types(self, mutation_types, proportions, genomic_element_type_id):
        """
        Adds mutation types with their respective proportions to the genomic
        element type provided.
        """
        _check_mut_proportions(proportions, mutation_types)
        if genomic_element_type_id >= len(self.genomic_element_types):
            raise ValueError(
                f"Genomic element type of id {genomic_element_type_id} does not exist."
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

    def add_genomic_element_type(self, intervals, mutation_types, proportions):
        if self.recombination_map is None:
            raise ValueError(
                "You cannot add genomic element types to a "
                "contig that has no recombination map."
            )
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
            self.recombination_map.mean_recombination_rate,
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
