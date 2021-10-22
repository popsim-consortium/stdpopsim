import attr


@attr.s(kw_only=True)
class MutationType(object):
    """
    TODO: docstring.

    See SLiM documentation for initializeMutationType().
    """

    dominance_coeff = attr.ib(default=0.5, type=float)
    distribution_type = attr.ib(default="f", type=str)
    distribution_args = attr.ib(factory=lambda: [0], type=list)
    convert_to_substitution = attr.ib(default=True, type=bool)

    def __attrs_post_init__(self):
        if self.dominance_coeff < 0:
            raise ValueError(f"Invalid dominance coefficient {self.dominance_coeff}.")

        # To add a new distribution type: validate the
        # distribution_args here, and add unit tests.
        if self.distribution_type == "f":
            # Fixed-value selection coefficent.
            if len(self.distribution_args) != 1:
                raise ValueError(
                    "Fixed-value mutation types (distribution_type='f')"
                    "take a single selection-coefficient parameter."
                )
        elif self.distribution_type == "g":
            # Gamma-distributed selection coefficient with (mean, shape)
            # parameterisation. A negative value for the mean is permitted,
            # and indicates a reflection of the horizontal axis.
            # See Eidos documentation for rgamma().
            if len(self.distribution_args) != 2 or self.distribution_args[1] <= 0:
                raise ValueError(
                    "Gamma-distributed sel. coefs. (distribution_type='g') "
                    "use a (mean, shape) parameterisation, requiring shape > 0."
                )
        elif self.distribution_type == "e":
            # An exponentially-distributed fitness effect (mean).
            # See Eidos documentation for rexp().
            if len(self.distribution_args) != 1:
                raise ValueError(
                    "Exponentially-distributed sel. coefs. (distribution_type='e') "
                    "use a (mean) parameterisation."
                )
        elif self.distribution_type == "n":
            # An normally-distributed fitness effect (mean, standard deviation).
            # See Eidos documentation for rnorm().
            if len(self.distribution_args) != 2 or self.distribution_args[1] <= 0:
                raise ValueError(
                    "Normally-distributed sel. coefs. (distribution_type='n') "
                    "use a (mean, sd) parameterisation, requiring sd > 0."
                )
        elif self.distribution_type == "w":
            # A Weibull-distributed fitness effect (scale, shape).
            # See Eidos documentation for rweibull().
            if (
                len(self.distribution_args) != 2
                or self.distribution_args[0] <= 0
                or self.distribution_args[1] <= 0
            ):
                raise ValueError(
                    "Weibull-distributed sel. coef. (distribution_type='w') "
                    "use a (scale, shape) parameterisation, requiring parameters > 0."
                )
        elif self.distribution_type == "l":
            # An lognormally-distributed fitness effect (logmean, sdmean).
            # See Eidos documentation for rlnorm().
            if len(self.distribution_args) != 2 or self.distribution_args[1] <= 0:
                raise ValueError(
                    "Lognormally-distributed sel. coefs. (distribution_type='l') "
                    "use a (logmean, logsd) parameterisation, requiring logsd > 0."
                )
            self.distribution_type = "s"
            # dealing with lognormal distribution
            # (adding instead of multiplying the mean):
            logmean = self.distribution_args[0]
            logsd = self.distribution_args[1]
            self.distribution_args = f'"return rlnorm(1, {logmean} + log(Q), {logsd});"'
        else:
            raise ValueError(
                f"{self.distribution_type} is not a supported distribution type"
            )

        # The index(s) of the param in the distribution_args list that should be
        # multiplied by Q when using --slim-scaling-factor Q.
        self.Q_scaled_index = {
            "e": [0],  # mean
            "f": [0],  # fixed value
            "g": [0],  # mean
            "n": [0, 1],  # mean and sd
            "w": [0],  # scale
            "s": [],  # script types should just printout arguments
        }[self.distribution_type]

    def is_neutral(self):
        """
        tests if the mutation type is strictly neutral
        """
        return self.distribution_type == "f" and self.distribution_args[0] == 0


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
    for SLiM things that msprime doesn't support. Times are all in units of
    generations before the present, just like for msprime.DemographicEvent.
    """

    pass


@attr.s(kw_only=True)
class DrawMutation(ExtendedEvent):
    """
    TODO: docstring.

    Draw a new mutation with the given mutation type in the given population
    at the given genomic coordinate. The mutation will be added to one randomly
    chosen chromosome in the given population.
    If save=True, the simulation state is saved before the mutation is
    introduced, to facilitate the save/restore mechanism used for allele
    frequency conditioning.

    FIXME: Drawing multiple mutations using the same mutation type causes
           erroneous allele frequency calculation in the SLiM code!
           Multiple drawn mutations should be disallowed, or the allele frequency
           calculation should somehow be tied to specific drawn mutations.

    See SLiM documentation for addNewDrawnMutation().
    """

    time = attr.ib(type=float)
    mutation_type_id = attr.ib(type=int)
    population_id = attr.ib(type=int)
    coordinate = attr.ib(type=int)
    save = attr.ib(type=bool, default=False)

    def __attrs_post_init__(self):
        validate_time(self.time)


@attr.s(kw_only=True)
class ChangeMutationFitness(ExtendedEvent):
    """
    TODO: docstring.

    Change the selection and/or dominance coefficients for the given mutation
    type in the given population. The mutation this applies to need not exist
    in the population for all generations between start_time and end_time.
    E.g. if the mutation could be transferred to the population via migration.

    See SLiM documentation for registerFitnessCallback().
    """

    start_time = attr.ib(type=float)
    end_time = attr.ib(type=float)
    mutation_type_id = attr.ib(type=int)
    population_id = attr.ib(type=int)
    selection_coeff = attr.ib(type=float)
    dominance_coeff = attr.ib(type=float)

    def __attrs_post_init__(self):
        validate_time_range(self.start_time, self.end_time)


@attr.s(kw_only=True)
class ConditionOnAlleleFrequency(ExtendedEvent):
    """
    TODO: docstring.

    Condition on the allele frequency of a drawn mutation with the given
    mutation type in the given population. The mutation need not have been
    drawn in the population being conditioned upon.
    If save=True, the simulation state will be saved if the condition is
    met at start_time.

    This uses a save/restore mechanism in the SLiM code to do rejection
    sampling on a simulation (when the condition is not met) without throwing
    away the entire simulation.
    """

    start_time = attr.ib(type=float)
    end_time = attr.ib(type=float)
    mutation_type_id = attr.ib(type=int)
    population_id = attr.ib(type=int)
    op = attr.ib(type=str, default=None)
    allele_frequency = attr.ib(type=float, default=None)
    save = attr.ib(type=bool, default=False)

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
