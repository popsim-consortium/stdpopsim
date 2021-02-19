"""
Infrastructure for defining information about species' genomes.
"""
import attr


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


@attr.s
class Contig:
    """
    Class representing a contiguous region of genome that is to be
    simulated. This contains the information about mutation rates
    and recombination rates that are needed to simulate this region.

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

    def __str__(self):
        gmap = "None" if self.genetic_map is None else self.genetic_map.id
        s = (
            "Contig(length={:.2G}, recombination_rate={:.2G}, "
            "mutation_rate={:.2G}, genetic_map={})"
        ).format(
            self.recombination_map.sequence_length,
            self.recombination_map.mean_rate,
            self.mutation_rate,
            gmap,
        )
        return s
