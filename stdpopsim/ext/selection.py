import attr
import tskit


@attr.s
class GenerationAfter(object):
    """
    We often want to do things one generation after another event. This class
    is essentially just a marker so that the SLiM engine can calculate the next
    generation by taking into account the scaling parameter Q.
    We call float(generation_after_obj) to get the time value.

    XXX: Rename this? GenerationAfter means the generation after something
         in a forwards-time convention, but times are specified with units of
         generations before the present!
    """

    time = attr.ib(type=float)

    def __attrs_post_init__(self):
        validate_time(self)

    def __float__(self):
        return float(self.time)


def validate_time(time):
    if float(time) < 0:
        raise ValueError("Negative times are invalid.")
    if float(time) == 0 and isinstance(time, GenerationAfter):
        raise ValueError("The generation after time=0 is in the future.")


def validate_time_range(start_time, end_time):
    validate_time(start_time)
    validate_time(end_time)
    if float(start_time) < float(end_time):
        raise ValueError(f"start_time={start_time} < end_time={end_time})")
    if (
        float(start_time) == float(end_time)
        and isinstance(start_time, GenerationAfter)
        and not isinstance(end_time, GenerationAfter)
    ):
        raise ValueError(
            f"start_time=GenerationAfter({start_time}) < end_time={end_time}"
        )


class ExtendedEvent(object):
    """
    ExtendedEvent is analogous to msprime.DemographicEvent, but is used here
    for explicitly manipulating mutations. Currently, these are only used by
    the SLiM engine. Times are all in units of generations before the present,
    just like for msprime.DemographicEvent.
    """

    pass


@attr.s(kw_only=True)
class DrawMutation(ExtendedEvent):
    """
    Introduce a new mutation with the given mutation type in the given
    population at a given single site. The mutation will be added to one
    randomly chosen chromosome in the given population. Only one mutation
    may be drawn per single site.

    :ivar time: The number of generations in the past at which to introduce the
        mutation.
    :vartype start_time: float
    :ivar single_site_id: The string ID of the single site at which to
        introduce the mutation, see :meth:`Contig.add_single_site`.
    :vartype single_site_id: str
    :ivar population: The name of the population in which to introduce
        the mutation.
    :vartype population: str
    """

    time = attr.ib(type=float)
    single_site_id = attr.ib(type=str)
    population = attr.ib(type=str)

    def __attrs_post_init__(self):
        validate_time(self.time)


@attr.s(kw_only=True)
class ChangeMutationFitness(ExtendedEvent):
    """
    Change the selection and/or dominance coefficients for the given mutation
    type in the given population. The mutation this applies to need not exist
    in the population for all generations between start_time and end_time.
    E.g. if the mutation could be transferred to the population via migration.
    This requires that a `DrawMutation` event with the same `single_site_id` is
    also supplied to the simulation engine.

    See SLiM documentation for registerFitnessCallback().

    :ivar start_time: The number of generations in the past at which the change
        in mutation fitness begins.
    :vartype start_time: float
    :ivar end_time: The number of generations in the past at which the change
        in mutation fitness ends.
    :vartype end_time: float
    :ivar single_site_id: The string ID of the single site at which to change
        the mutation fitness, see :meth:`Contig.add_single_site`.
    :vartype single_site_id: str
    :ivar population: The name of the population in which to change
        the mutation fitness. If `None` (the default), then all populations are
        affected.
    :vartype population: str
    :ivar selection_coeff: The new selection coefficient of the mutation.
    :vartype selection_coeff: float
    :ivar dominance_coeff: The new dominance coefficient of the mutation.
    :vartype dominance_coeff: float
    """

    start_time = attr.ib(type=float)
    end_time = attr.ib(type=float)
    single_site_id = attr.ib(type=str)
    population = attr.ib(type=str, default=None)
    selection_coeff = attr.ib(type=float)
    dominance_coeff = attr.ib(type=float)

    def __attrs_post_init__(self):
        validate_time_range(self.start_time, self.end_time)


@attr.s(kw_only=True)
class ConditionOnAlleleFrequency(ExtendedEvent):
    """
    Condition on the allele frequency of a drawn mutation with the given
    mutation type in the given population. The mutation need not have been
    drawn in the population being conditioned upon.  This requires that a
    `DrawMutation` event for the same single site is also supplied to the
    simulation engine.

    :ivar start_time: The number of generations in the past at which to start
        checking the allele frequency condition.
    :vartype start_time: float
    :ivar end_time: The number of generations in the past at which to stop
        checking the allele frequency condition.
    :vartype end_time: float
    :ivar single_site_id: The string ID of the single site at which to check
        the allele frequency condition, see :meth:`Contig.add_single_site`.
    :vartype single_site_id: str
    :ivar population: The name of the population in which to check
        the allele frequency condition.
    :vartype population: str
    :ivar op: The comparison operator used for the allele frequency condition,
        one of `(">", ">=", "<", "<=")`.
    :vartype op: str
    :ivar allele_frequency: The allele frequency to compare against.
    :vartype op: float
    """

    start_time = attr.ib(type=float)
    end_time = attr.ib(type=float)
    single_site_id = attr.ib(type=str)
    population = attr.ib(type=str, default=None)
    op = attr.ib(type=str, default=None)
    allele_frequency = attr.ib(type=float, default=None)

    op_types = ("<", "<=", ">", ">=")

    def __attrs_post_init__(self):
        if self.op not in self.op_types:
            raise ValueError(f"Invalid conditioning op `{self.op}`.")
        if not (0 <= self.allele_frequency <= 1):
            raise ValueError("Must have 0 <= allele_frequency <= 1.")
        if (self.op == "<" and self.allele_frequency == 0) or (
            self.op == ">" and self.allele_frequency == 1
        ):
            raise ValueError(
                f"allele_frequency {self.op} {self.allele_frequency}: "
                "condition is always false."
            )
        if (self.op == ">=" and self.allele_frequency == 0) or (
            self.op == "<=" and self.allele_frequency == 1
        ):
            raise ValueError(
                f"allele_frequency {self.op} {self.allele_frequency}: "
                "condition is always true."
            )
        validate_time_range(self.start_time, self.end_time)

    @classmethod
    def op_id(cls, op):
        return cls.op_types.index(op)


def selective_sweep(
    single_site_id,
    population,
    mutation_generation_ago,
    selection_coeff,
    dominance_coeff=0.5,
    start_generation_ago=None,
    end_generation_ago=None,
    min_freq_at_start=None,
    min_freq_at_end=None,
    globally_adaptive=True,
):
    """
    Creates a list of extended events corresponding to a selective sweep at a
    single site (see :meth:`Contig.add_single_site`) that may subsequently be
    passed to the simulation engine. The sequence of events is:

        (1) In generation `mutation_generation_ago`, a **neutral** mutation is
            introduced into `population`.

        (2) In generation `start_generation_ago`, the mutation becomes
            beneficial (whether globally or locally is controlled by the
            `globally_adaptive` option). If the frequency of the mutation in
            `population` is below `min_freq_at_start` in this generation, the
            simulation will restart from `mutation_generation_ago`.

        (3) In the generation after `end_generation_ago`, the mutation reverts to
            neutrality in `population`. If the frequency of the mutation in
            `population` is below `min_freq_at_end` in this generation, the
            simulation will restart from `mutation_generation_ago`.

    During the sweep, the fitnesses of individuals that are homozygous or
    heterozygous for the beneficial mutation are `1 + selection_coeff` and  `1
    + selection_coeff * dominance_coeff`, respectively. If `globally_adaptive`
    is `True` (the default), then all individuals carrying the mutation are
    impacted regardless of the population they belong to. Otherwise, the
    mutation only influences fitness within `population`.

    If the mutation is lost in `population` at any point between
    `mutation_generation_ago` and `end_generation_ago`, then the simulation
    will restart from `mutation_generation_ago`. The effect of "restarting"
    (aka "rejection sampling") is that the resulting simulations are
    conditioned on the mutation being extant after time
    `mutation_generation_ago` with frequency above `min_freq_at_start` at
    `start_generation_ago` and `min_freq_at_end` at `end_generation_ago`.

    Because the simulation will restart until the allele frequency conditions
    are satisfied, care must be taken to ensure that these are reasonable. For
    example, if `start_generation_ago` is very close to
    `mutation_generation_ago` and `min_freq_at_start` is substantially greater
    than zero, the vast majority of proposed allele frequency trajectories will
    be rejected and the simulation may take an infeasibly long time to
    complete.

    :param str single_site_id: The string ID of the single site at which to introduce
        the mutation, see :meth:`Contig.add_single_site`.
    :param int population: The name of the population in which the mutation is
        introduced and the sweep occurs.
    :param float mutation_generation_ago: The number of generations in the past at which
        to introduce the mutation.
    :param float selection_coeff: The selection coefficient of the beneficial
        mutation during the sweep.  Must be greater than or equal to zero.
    :param float dominance_coeff: The dominance coefficient of the beneficial
        mutation during the sweep.  Defaults to `0.5`.
    :param float start_generation_ago: The number of generations in the past at which the
        sweep starts and the mutation switches from neutral to beneficial. Defaults to
        `mutation_generation_ago` (i.e. a "hard" sweep).
    :param float end_generation_ago: The number of generations in the past at which the
        sweep ends and the mutation switches from beneficial to neutral. Defaults to
        `0` (i.e. the end of the simulation).
    :param float min_freq_at_start: The minimum allowed frequency for the
        mutation in `population` at the start of the sweep (optional).
    :param float min_freq_at_end: The minimum allowed frequency for the
        mutation in `population` at the end of the sweep (optional).
    :param bool globally_adaptive: If `True` (default), the mutation is beneficial in all
        populations. If `False`, the mutation is beneficial only in the population of
        origin.
    """

    if start_generation_ago is None:
        start_generation_ago = mutation_generation_ago
    if end_generation_ago is None:
        end_generation_ago = 0

    validate_time_range(mutation_generation_ago, start_generation_ago)
    validate_time_range(start_generation_ago, end_generation_ago)
    if selection_coeff < 0.0:
        raise ValueError("Selection coefficient must be non-negative.")

    extended_events = [
        DrawMutation(
            single_site_id=single_site_id,
            time=mutation_generation_ago,
            population=population,
        ),
        ConditionOnAlleleFrequency(
            single_site_id=single_site_id,
            start_time=GenerationAfter(mutation_generation_ago),
            end_time=end_generation_ago,
            population=population,
            allele_frequency=0.0,
            op=">",
        ),
    ]

    # Sweep period is [start_generation_ago, end_generation_ago] (endpoints inclusive).
    sweep_restricted_to = None if globally_adaptive else population
    extended_events.extend(
        [
            ChangeMutationFitness(
                single_site_id=single_site_id,
                start_time=start_generation_ago,
                end_time=end_generation_ago,
                population=sweep_restricted_to,
                selection_coeff=selection_coeff,
                dominance_coeff=dominance_coeff,
            ),
        ]
    )

    # Only check sweep AF conditions in start/end generations, so that allele
    # trajectories are otherwise unconstrained.
    if min_freq_at_start is not None:
        if not 0 < min_freq_at_start <= 1:
            raise ValueError(
                "If specified, the minimum allele frequency at the start of the sweep "
                "must be in (0, 1]."
            )
        if start_generation_ago == mutation_generation_ago:
            raise ValueError(
                "If the start time of the sweep coincides with the introduction "
                "of the mutation, then `min_freq_at_start` cannot be set."
            )
        extended_events.extend(
            [
                ConditionOnAlleleFrequency(
                    single_site_id=single_site_id,
                    start_time=start_generation_ago,
                    end_time=start_generation_ago,
                    population=population,
                    allele_frequency=min_freq_at_start,
                    op=">=",
                ),
            ]
        )

    if min_freq_at_end is not None:
        if not 0 < min_freq_at_end <= 1:
            raise ValueError(
                "If specified, the minimum allele frequency at the end of the sweep "
                "must be in (0, 1]."
            )
        extended_events.extend(
            [
                ConditionOnAlleleFrequency(
                    single_site_id=single_site_id,
                    start_time=end_generation_ago,
                    end_time=end_generation_ago,
                    population=population,
                    allele_frequency=min_freq_at_end,
                    op=">=",
                ),
            ]
        )

    return extended_events


def selection_coeff_from_mutation(ts, mutation):
    """
    Extract the selection coefficient from a (possibly stacked) mutation.

    :param ts: A ``tskit.TreeSequence`` containing the mutation.
    :param mutation: A ``tskit.Mutation`` for which to extract the selection coefficient.
    """

    if not isinstance(ts, tskit.TreeSequence):
        raise ValueError("`ts` must be a `tskit.TreeSequence` object")

    if not isinstance(mutation, tskit.Mutation):
        raise ValueError("`mutation` must be a `tskit.Mutation` object")

    if not isinstance(mutation.metadata, dict):
        return 0.0

    selection_coeff = sum(
        [m.get("selection_coeff") for m in mutation.metadata["mutation_list"]]
    )
    if mutation.parent != tskit.NULL:
        parent = ts.mutation(mutation.parent)
        selection_coeff -= sum(
            [m.get("selection_coeff") for m in parent.metadata["mutation_list"]]
        )

    return selection_coeff
