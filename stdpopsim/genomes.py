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
    :ivar bacterial_recombination: Whether recombination is via horizontal gene
        transfer (if this is True) or via crossover and possibly gene
        conversion (if this is False). Default: False.
    :vartype bacterial_recombination: bool
    :ivar citations: A list of :class:`.Citation` objects
        providing the source for the genome assembly,
        mutation rate, recombination rate, and gene conversion estimates.
    :vartype citations: list
    :ivar length: The total length of the genome.
    :vartype length: int
    """

    # TODO document the assembly_name and accession

    chromosomes = attr.ib(factory=list)
    assembly_name = attr.ib(default=None, kw_only=True)
    assembly_accession = attr.ib(default=None, kw_only=True)
    bacterial_recombination = attr.ib(type=bool, default=False, kw_only=True)
    citations = attr.ib(factory=list, kw_only=True)

    @staticmethod
    def from_data(
        genome_data,
        *,
        recombination_rate,
        mutation_rate,
        citations,
        bacterial_recombination=False,
        gene_conversion_fraction=None,
        gene_conversion_length=None,
    ):
        """
        Construct a Genome object from the specified dictionary of
        genome information from Ensembl, recombination_rate and
        mutation_rate dictionaries.

        This method is for internal use only.
        """
        chr_names = set(genome_data["chromosomes"].keys())
        assert set(recombination_rate.keys()) == chr_names
        assert set(mutation_rate.keys()) == chr_names
        if gene_conversion_fraction is None:
            gene_conversion_fraction = {k: None for k in chr_names}
        assert set(gene_conversion_fraction.keys()) == chr_names
        if gene_conversion_length is None:
            gene_conversion_length = {k: None for k in chr_names}
        assert set(gene_conversion_length.keys()) == chr_names
        chromosomes = []
        for name, data in genome_data["chromosomes"].items():
            chromosomes.append(
                Chromosome(
                    id=name,
                    length=data["length"],
                    synonyms=data["synonyms"],
                    mutation_rate=mutation_rate[name],
                    recombination_rate=recombination_rate[name],
                    gene_conversion_fraction=gene_conversion_fraction[name],
                    gene_conversion_length=gene_conversion_length[name],
                )
            )
        return Genome(
            chromosomes=chromosomes,
            assembly_name=genome_data["assembly_name"],
            assembly_accession=genome_data["assembly_accession"],
            bacterial_recombination=bacterial_recombination,
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
    def mean_gene_conversion_fraction(self):
        """
        The genetic length-weighted mean gene conversion fraction across all
        chromosomes.
        """
        if self.bacterial_recombination is True:
            return 0.0
        else:
            total_length = 0
            total_gc_fraction = 0
            for chrom in self.chromosomes:
                length = chrom.length * chrom.recombination_rate
                total_length += length
                if chrom.gene_conversion_fraction is None:
                    frac = 0.0
                else:
                    frac = chrom.gene_conversion_fraction
                total_gc_fraction += frac * length
            return total_gc_fraction / total_length

    @property
    def range_gene_conversion_lengths(self):
        """
        The range of gene conversion tract lengths across chromosomes.
        """
        prelim_lengths = [chrom.gene_conversion_length for chrom in self.chromosomes]
        gc_tract_lengths = [0 if gctl is None else gctl for gctl in prelim_lengths]
        min_gctl = min(gc_tract_lengths)
        max_gctl = max(gc_tract_lengths)
        return f"{min_gctl}" if min_gctl == max_gctl else f"{min_gctl} - {max_gctl}"

    @property
    def mean_recombination_rate(self):
        """
        The length-weighted mean recombination rate across all chromosomes.
        (Note that this is the rate of crossovers, not all double-stranded breaks.)
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
    Class representing a single chromosome for a species. Note that although
    the recombination rate for a Chromosome is the rate of crossovers,
    the recombination map for a Contig gives the rate of both crossovers
    and gene conversion, so will differ if the gene conversion fraction is nonzero.

    :ivar str ~.id: The string identifier for the chromosome.
    :ivar int length: The length of the chromosome.
    :ivar float mutation_rate: The mutation rate used when simulating this
        chromosome.
    :ivar float recombination_rate: The rate of crossovers used when simulating
        this chromosome (if not using a genetic map).
    :ivar float gene_conversion_fraction: The gene conversion fraction used when
        simulating this chromosome.
    :ivar float gene_conversion_length: The mean tract length of gene
        conversion events used when simulating this chromosome.
    :ivar synonyms: List of synonyms that may be used when requesting this
        chromosome by ID, e.g. from the command line interface.
    :vartype synonyms: list of str
    """

    id = attr.ib(type=str, kw_only=True)
    length = attr.ib(kw_only=True)
    recombination_rate = attr.ib(type=float, kw_only=True)
    gene_conversion_fraction = attr.ib(default=None, type=float, kw_only=True)
    gene_conversion_length = attr.ib(default=None, type=float, kw_only=True)
    mutation_rate = attr.ib(type=float, kw_only=True)
    synonyms = attr.ib(factory=list, kw_only=True)


@attr.s(kw_only=True)
class Contig:
    """
    Class representing a contiguous region of genome that is to be simulated.
    This contains the information about mutation rates, distributions of
    fitness effects (DFEs), gene conversion rates and recombination rates
    that are needed to simulate this region.

    Information about targets of selection are contained in the ``dfe_list``
    and ``interval_list``. These must be of the same length, and the k-th DFE
    applies to the k-th interval; see :meth:`.add_dfe` for more information.

    The contig may be a segment of a named chromosome. If so, the original
    coordinate system is used for :meth:`.add_dfe` and :meth:`add_single_site`.

    :ivar mutation_rate: The rate of mutation per base per generation.
    :vartype mutation_rate: float
    :ivar recombination_map: The recombination map for the region, that gives
        the rates of double-stranded breaks. Note that if gene conversion
        fraction is nonzero, then this map can have a larger rate than the
        corresponding chromosome's recombination rate, since the chromosome's
        rate describes only crossovers. A :class:`msprime.RateMap` object.
    :vartype recombination_map: msprime.RateMap
    :ivar bacterial_recombination: Whether the model of recombination is
        by horizontal gene transfer or not (default is False, i.e.,
        recombination is by crossovers and possibly gene conversion).
    :ivar gene_conversion_fraction: The fraction of recombinations that resolve
        as gene conversion events. The recombination map gives the rates of
        *both* crossovers and gene conversions, with a fraction of gene
        conversions given by this value. Defaults to None, i.e., no gene
        conversion. Must be None or between 0 and 1. Must be None if
        bacterial_recombination is True.
    :vartype gene_conversion_fraction: float
    :ivar gene_conversion_length: The mean tract length of gene conversions (if
        bacterial_recombination is False), or horizontally tranferred segments (if
        it is True). Must be present, and above 1 if either gene_conversion_fraction
        is given or bacterial_recombination is True.
    :vartype gene_conversion_length: float
    :ivar mask_intervals: Intervals to keep in simulated tree sequence, as a list
        of (left_position, right_position), such that intervals are non-overlapping
        and in ascending order. Should have shape Nx2, where N is the number of
        intervals.
    :vartype mask_intervals: array-like (?)
    :ivar exclude: If True, ``mask_intervals`` specify regions to exclude. If False,
        ``mask_intervals`` specify regions in keep.
    :vartype exclude: bool
    :ivar dfe_list: A list of :class:`.DFE` objects.
        By default, the only DFE is completely neutral.
    :vartype dfe_list: list
    :ivar interval_list: A list of :class:`np.array` objects containing integers.
        By default, the inital interval list spans the whole chromosome with the
        neutral DFE.
    :ivar original_coordinates: The location of the contig on a named chromosome,
        as a tuple of the form `(chromosome, left, right)`. If `None`, the contig
        is assumed to be generic (i.e. it does not inherit a coordinate system
        from a larger chromosome).
    :vartype original_coordinates: tuple

    .. note::
        To run stdpopsim simulations with alternative, user-specified mutation,
        recombination, or gene conversion rates, a new contig can be created
        based on an existing one. For instance, the following will create a
        ``new_contig`` that, when simulated, will have double the mutation rate
        of the ``old_contig``:

        .. code-block:: python

            new_contig = stdpopsim.Contig(
                mutation_rate=old_contig.mutation_rate * 2,
                recombination_map=old_contig.recombination_map,
                genetic_map=old_contig.genetic_map,
            )
    """

    recombination_map = attr.ib()
    mutation_rate = attr.ib(type=float)
    bacterial_recombination = attr.ib(default=False, type=bool, kw_only=True)
    gene_conversion_fraction = attr.ib(default=None, type=float, kw_only=True)
    gene_conversion_length = attr.ib(default=None, type=float, kw_only=True)
    genetic_map = attr.ib(default=None)
    inclusion_mask = attr.ib(default=None)
    exclusion_mask = attr.ib(default=None)
    dfe_list = attr.ib(factory=list)
    interval_list = attr.ib(factory=list)
    original_coordinates = attr.ib(default=None, type=tuple)

    def __attrs_post_init__(self):
        if self.original_coordinates is None:
            self.original_coordinates = (None, 0, int(self.length))
        _, left, right = self.original_coordinates
        self.add_dfe(
            [[left, right]],
            stdpopsim.dfe.neutral_dfe(),
        )

    @staticmethod
    def basic_contig(
        *,
        length,
        mutation_rate=0,
        recombination_rate=0,
        bacterial_recombination=False,
        gene_conversion_fraction=None,
        gene_conversion_length=None,
    ):
        if bacterial_recombination:
            if gene_conversion_fraction is not None:
                raise ValueError(
                    "Cannot set gene conversion fraction for bacterial recombination."
                )
            if gene_conversion_length is None:
                raise ValueError(
                    "Must set gene conversion length for bacterial recombination."
                )
        else:
            if gene_conversion_fraction is not None:
                if gene_conversion_length is None:
                    raise ValueError(
                        "Cannot set gene conversion fraction without setting "
                        "gene conversion length"
                    )
            if gene_conversion_length is not None:
                if gene_conversion_fraction is None:
                    raise ValueError(
                        "Cannot set gene conversion length without setting "
                        "gene conversion fraction"
                    )
        if (gene_conversion_fraction is not None) and (
            gene_conversion_fraction < 0 or gene_conversion_fraction > 1
        ):
            raise ValueError("Gene conversion fraction must be between 0 and 1.")
        if (gene_conversion_length is not None) and (gene_conversion_length < 1):
            # values of even 0 are OK for SLiM but not msprime
            raise ValueError("Gene conversion length must be greater than 1.")
        recomb_map = msprime.RateMap.uniform(length, recombination_rate)
        return Contig(
            recombination_map=recomb_map,
            mutation_rate=mutation_rate,
            bacterial_recombination=bacterial_recombination,
            gene_conversion_fraction=gene_conversion_fraction,
            gene_conversion_length=gene_conversion_length,
        )

    @staticmethod
    def species_contig(
        *,
        species,
        chromosome=None,
        genetic_map=None,
        length_multiplier=1,
        length=None,
        mutation_rate=None,
        use_species_gene_conversion=False,
        inclusion_mask=None,
        exclusion_mask=None,
        left=None,
        right=None,
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
        if (
            species.genome.bacterial_recombination is True
            and use_species_gene_conversion is True
        ):
            raise ValueError(
                "Cannot use species gene conversion with bacterial recombination."
            )
        gene_conversion_fraction = None
        gene_conversion_length = None
        if chromosome is None:
            if left is not None or right is not None:
                raise ValueError(
                    "Cannot use left or right coordinates with generic contig"
                )
            if genetic_map is not None:
                raise ValueError("Cannot use genetic map with generic contig")
            if length_multiplier != 1:
                raise ValueError("Cannot use length multiplier for generic contig")
            if inclusion_mask is not None or exclusion_mask is not None:
                raise ValueError("Cannot use mask with generic contig")
            if length is None:
                raise ValueError("Must specify sequence length of generic contig")
            L_tot = 0
            r_tot = 0
            u_tot = 0
            gcf_tot = 0
            gcl_tot = 0
            for chrom_data in species.genome.chromosomes:
                if chrom_data.id.lower() not in non_autosomal_lower:
                    L_tot += chrom_data.length
                    r_tot += chrom_data.length * chrom_data.recombination_rate
                    u_tot += chrom_data.length * chrom_data.mutation_rate
                    if chrom_data.gene_conversion_fraction is not None:
                        gcf_tot += (
                            chrom_data.length
                            * chrom_data.recombination_rate
                            * chrom_data.gene_conversion_fraction
                        )
                    if chrom_data.gene_conversion_length is not None:
                        gcl_tot += (
                            chrom_data.length
                            * chrom_data.recombination_rate
                            * chrom_data.gene_conversion_length
                        )
            if mutation_rate is None:
                mutation_rate = u_tot / L_tot
            r = r_tot / L_tot
            if use_species_gene_conversion is True:
                if gcf_tot > 0:
                    gene_conversion_fraction = gcf_tot / r_tot
                    gene_conversion_length = gcl_tot / r_tot
                    r /= 1 - gene_conversion_fraction
            if species.genome.bacterial_recombination is True:
                gene_conversion_length = gcl_tot / r_tot
            contig = Contig.basic_contig(
                length=length,
                mutation_rate=mutation_rate,
                recombination_rate=r,
                bacterial_recombination=species.genome.bacterial_recombination,
                gene_conversion_fraction=gene_conversion_fraction,
                gene_conversion_length=gene_conversion_length,
            )
        else:
            # Named contig:
            if length is not None:
                raise ValueError("Cannot specify sequence length for named contig")
            if inclusion_mask is not None and exclusion_mask is not None:
                raise ValueError("Cannot specify both inclusion and exclusion masks")
            if length_multiplier != 1 and (left is not None or right is not None):
                raise ValueError(
                    "Cannot use length multiplier when specifying left or "
                    "right coordinates of the contig on a named chromosome."
                )

            chrom = species.genome.get_chromosome(chromosome)
            if left is None:
                left = 0
            else:
                left = round(left)
                if not 0 <= left < chrom.length:
                    raise ValueError(
                        f"Left coordinate {left} falls outside chromosome "
                        f"{chrom.id} with length {chrom.length}."
                    )
            if right is None:
                right = round(chrom.length * length_multiplier)
            else:
                right = round(right)
                if not left < right <= chrom.length:
                    raise ValueError(
                        f"Right coordinate {right} falls outside chromosome "
                        f"{chrom.id} with length {chrom.length}."
                    )

            if genetic_map is None:
                gm = None
                if length_multiplier != 1:
                    logger.debug(
                        f"Making flat chromosome {length_multiplier} * {chrom.id}"
                    )
                else:
                    logger.debug(
                        f"Making flat contig of length {right - left} from {chrom.id}"
                    )
                recomb_map = msprime.RateMap.uniform(
                    right - left, chrom.recombination_rate
                )
            else:
                if length_multiplier != 1:
                    raise ValueError("Cannot use length multiplier with empirical maps")
                logger.debug(f"Getting map for {chrom.id} from {genetic_map}")
                gm = species.get_genetic_map(genetic_map)
                recomb_map = gm.get_chromosome_map(chrom.id)
                recomb_map = recomb_map.slice(left=left, right=right, trim=True)

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
                inclusion_intervals = stdpopsim.utils.clip_and_shift_intervals(
                    inclusion_intervals, left, right
                )
            if exclusion_mask is not None:
                if length_multiplier != 1:
                    raise ValueError("Cannot use length multiplier with mask")
                if isinstance(exclusion_mask, str):
                    exclusion_intervals = stdpopsim.utils.read_bed(
                        exclusion_mask, chromosome
                    )
                else:
                    exclusion_intervals = exclusion_mask
                exclusion_intervals = stdpopsim.utils.clip_and_shift_intervals(
                    exclusion_intervals, left, right
                )

            if mutation_rate is None:
                mutation_rate = chrom.mutation_rate

            if species.genome.bacterial_recombination is True:
                gene_conversion_length = chrom.gene_conversion_length

            if use_species_gene_conversion is True:
                gene_conversion_fraction = chrom.gene_conversion_fraction
                gene_conversion_length = chrom.gene_conversion_length
                if (
                    gene_conversion_fraction is not None
                    and gene_conversion_fraction > 0
                ):
                    recomb_map = msprime.RateMap(
                        position=recomb_map.position,
                        rate=recomb_map.rate / (1 - gene_conversion_fraction),
                    )

            contig = stdpopsim.Contig(
                recombination_map=recomb_map,
                mutation_rate=mutation_rate,
                genetic_map=gm,
                bacterial_recombination=species.genome.bacterial_recombination,
                gene_conversion_fraction=gene_conversion_fraction,
                gene_conversion_length=gene_conversion_length,
                inclusion_mask=inclusion_intervals,
                exclusion_mask=exclusion_intervals,
                original_coordinates=(chromosome, left, right),
            )

        return contig

    @property
    def length(self):
        return self.recombination_map.sequence_length

    @property
    def origin(self):
        """
        The location of the contig on a named chromosome as a string with format,
        "chromosome:left-right"; or None if a generic contig.
        """
        chromosome, left, right = self.original_coordinates
        if chromosome is None:
            return None
        else:
            return f"{chromosome}:{left}-{right}"

    def dfe_breakpoints(self, relative_coordinates=True):
        """
        Returns two things: the sorted vector of endpoints of all intervals across
        all DFEs in the contig, and a vector of integer labels for these intervals,
        saying which DFE goes with which interval.
        This provides a complementary method to tell which bit of the contig
        has which DFE attached, which may be more convenient than the list of two-column
        arrays provided by interval_list.

        Suppose there are n+1 unique interval endpoints across all the DFE intervals
        in :attr:`.interval_list`.  (If the intervals in that list
        cover the whole genome, the number of intervals is n.)
        This method returns a tuple of two things: ``breaks, dfe_labels``.
        "breaks" is the array containing those n+1 unique endpoints, in increasing order,
        and "dfe" is the array of length n containing the index of the DFE
        the applies to that interval. So, ``breaks[0]`` is always 0, and
        ``breaks[n+1]`` is always the length of the contig, and
        ``dfe_labels[k] = j`` if
        ``[breaks[k], breaks[k+1]]`` is an interval in ``contig.interval_list[j]``,
        i.e., if ``contig.dfe_list[j]`` applies to the interval starting at
        ``breaks[k]``.  Some intervals may not be covered by a DFE, in which
        case they will have the label ``-1`` (beware of python indexing!).

        :param bool relative_coordinates: If True, the returned breakpoints
            will be relative to the start of the contig, rather than to the
            start of the chromosome to which to contig belongs.

        :return: A tuple (breaks, dfe_labels).
        """
        breaks = np.unique(
            np.vstack(self.interval_list + [[[0, int(self.length)]]])  # also sorted
        )
        dfe_labels = np.full(len(breaks) - 1, -1, dtype="int")
        for j, intervals in enumerate(self.interval_list):
            dfe_labels[np.isin(breaks[:-1], intervals[:, 0], assume_unique=True)] = j
        if not relative_coordinates:
            _, left, _ = self.original_coordinates
            breaks += left
        return breaks, dfe_labels

    def clear_dfes(self):
        """
        Removes all DFEs from the contig (as well as the corresponding list of
        intervals).
        """
        self.dfe_list = []
        self.interval_list = []

    def add_dfe(self, intervals, DFE):
        """
        Adds the provided DFE to the contig, applying to the regions of the
        contig specified in ``intervals``. These intervals are also *removed*
        from the intervals of any previously-present DFEs - in other words,
        more recently-added DFEs take precedence. Each DFE added in this way
        carries its own MutationTypes: no mutation types are shared between
        DFEs added by different calls to ``add_dfe( )``, even if the same DFE
        object is added more than once.

        For instance, if we do
        ```
        a1 = np.array([[0, 100]])
        a2 = np.array([[50, 120]])
        contig.add_dfe(a1, dfe1)
        contig.add_dfe(a2, dfe2)
        ```
        then ``dfe1`` applies to the region from 0 to 50 and ``dfe2`` applies
        to the region from 50 to 120.

        Any of the ``intervals`` that fall outside of the contig will be
        clipped to the contig boundaries. If no ``intervals`` overlap the
        contig, an "empty" DFE will be added with a warning.

        :param array intervals: A valid set of intervals.
        :param DFE dfe: A DFE object.
        """
        _, left, right = self.original_coordinates
        intervals = stdpopsim.utils.clip_and_shift_intervals(intervals, left, right)
        stdpopsim.utils._check_intervals_validity(intervals, start=0, end=self.length)
        for j, ints in enumerate(self.interval_list):
            self.interval_list[j] = stdpopsim.utils.mask_intervals(ints, intervals)
        self.dfe_list.append(DFE)
        self.interval_list.append(intervals)

    def add_single_site(
        self,
        id,
        coordinate,
        description=None,
        long_description=None,
    ):
        """
        Adds a single site mutation class at the provided coordinate.
        The string label of the single site may be referenced by extended
        events passed to the simulation engine.

        For instance, to simulate a selective sweep in population `pop_0`
        starting at 1000 generations in the past:

        .. code-block:: python

            contig.add_single_site(
                id="hard_sweep",
                coordinate=mutation_coordinate,
            )
            extended_events = ext.selective_sweep(
                mutation_generation_ago=1000,
                selection_coeff=0.1,
                population="pop_0",
            )
            engine = stdpopsim.get_engine("slim")
            ts_sweep = engine.simulate(
                ...,
                extended_events=extended_events,
            )

        See :func:`ext.selective_sweep` for more details on the sweep model.

        :param str id: A unique identifier for the single site.
        :param int coordinate: The coordinate of the site on the contig.
        :param str description: A short description of the single site mutation
            model as it would be used in written text, e.g., "Strong selective
            sweep".
        :param str long_description: If necessary, a more detailed summary.
        """
        if description is None:
            description = "Single site mutation"
        if long_description is None:
            long_description = "Added by Contig.add_single_site"

        mut_type = stdpopsim.MutationType(
            distribution_type="f",
            dominance_coeff=0.5,
            distribution_args=[0.0],
            convert_to_substitution=False,
        )

        dfe = stdpopsim.DFE(
            id=id,
            mutation_types=[mut_type],
            proportions=[1.0],
            description=description,
            long_description=long_description,
        )
        self.add_dfe(
            intervals=np.array([[coordinate, coordinate + 1]], dtype="int"),
            DFE=dfe,
        )

    @property
    def is_neutral(self):
        """
        Returns True if the contig has no non-neutral mutation types.
        """
        return all(mt.is_neutral for d in self.dfe_list for mt in d.mutation_types)

    def mutation_types(self):
        """
        Provides information about the MutationTypes assigned to this Contig,
        along with information about which DFE they correspond to. This is
        useful because when simulating with SLiM, the mutation types are
        assigned numeric IDs in the order provided here (e.g., in the order
        encountered when iterating over the DFEs); this method provides an easy
        way to map back from their numeric ID to the DFE that each MutationType
        corresponds to.

        This method returns a list of dictionaries of length equal to the
        number of MutationTypes in all DFEs of the contig, each dictionary
        containing three things: ``"dfe_id"``: the ID of the DFE this
        MutationType comes from; ``mutation_type``: the mutation type; and
        ``id``: the index in the list.

        For instance, if ``muts`` is a list of mutation objects in a SLiM tree
        sequence, then the following code will print the IDs of the DFEs that
        each comes from:

        .. code-block:: python

            mut_types = contig.mutation_types()
            for m in muts:
                dfe_ids = [mut_types[k]["dfe_id"] for md in m.metadata["muation_list"]]
                print(f"Mutation {m.id} has mutations from DFE(s) {','.join(dfe_ids)")

        """
        id = 0
        mut_types = []
        for d in self.dfe_list:
            for mt in d.mutation_types:
                mut_types.append(
                    {
                        "dfe_id": d.id,
                        "mutation_type": mt,
                        "id": id,
                    }
                )
                id += 1
        return mut_types

    def __str__(self):
        gmap = "None" if self.genetic_map is None else self.genetic_map.id
        s = (
            "Contig(length={:.2G}, recombination_rate={:.2G}, "
            "mutation_rate={:.2G}, bacterial_recombination={}, "
            "gene_conversion_fraction={}, gene_conversion_length={}, "
            "genetic_map={}, origin={}, dfe_list={}, "
            "interval_list={})"
        ).format(
            self.length,
            self.recombination_map.mean_rate,
            self.mutation_rate,
            self.bacterial_recombination,
            self.gene_conversion_fraction,
            self.gene_conversion_length,
            gmap,
            self.origin,
            self.dfe_list,
            self.interval_list,
        )
        return s
