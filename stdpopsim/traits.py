"""
Methods related to traits and effects of mutations on them,
including environment and fitness (so, this includes DFE machinery).
"""

import copy
import textwrap
import attr
import collections.abc
import numpy as np


def _copy_converter(x):
    if isinstance(x, list):
        x = x.copy()
    return x


def _check_nonoverlapping_intervals(intervals):
    if not isinstance(intervals, list):
        raise ValueError("Intervals must be a list.")
    prev_end = -float("inf")
    if len(intervals) < 1:
        raise ValueError("Cannot supply an empty list of intervals.")
    for interval in sorted(intervals):
        if len(interval) != 2:
            raise ValueError(
                "Each interval in time_interval must be of the form "
                "(start, end), specified backward in time."
            )
        if not isinstance(interval[0], (int, float)):
            raise ValueError("Intervals must be numeric.")
        if not isinstance(interval[1], (int, float)):
            raise ValueError("Intervals must be numeric.")
        if interval[0] < 0:
            raise ValueError(
                "Intervals must start at the present or some more ancient time."
            )
        if interval[0] > interval[1]:
            raise ValueError(
                "Lower end of interval must be strictly less than upper end."
            )
        if interval[0] < prev_end:
            raise ValueError("Intervals must be non-overlapping.")
        prev_end = interval[1]


def _check_trait_ids(trait_ids):
    # this might be too strict but we want to avoid things like trait_ids="foo"
    # which might then be interpreted as three trait IDs: "f", "o", and "o"
    # if we just iterate over it
    if not isinstance(trait_ids, list):
        raise ValueError("Trait IDs must be a list.")
    for tid in trait_ids:
        if not (isinstance(tid, str) and (len(tid) > 0)):
            raise ValueError(f"Each trait ID must be a nonempty string; found {tid}.")
    if (len(trait_ids) == 0) or (len(set(trait_ids)) != len(trait_ids)):
        raise ValueError("Trait IDs must be a nonempty list of unique strings.")


def _find_oldest_event(traits_model, demography_debugger):
    # Returns a tuple containing (time_of_oldest_event, condition_id) where
    # time_of_oldest_event is the time (in generations) of the oldest (finite)
    # event that occurs either in the demography or in one of the fitness
    # functions or environments.  condition_id will be the str ID of the
    # fitness or function triggering the oldest event or will be the empty
    # string if the oldest event is driven by the demography.
    oldest_event = max(e.start_time for e in demography_debugger.epochs)
    oldest_time = oldest_event
    oldest_id = ""
    conditions = traits_model.environments + traits_model.fitness_functions
    for condition in conditions:
        for interval in condition.time_intervals:
            for t in interval:
                if np.isfinite(t) and t > oldest_time:
                    oldest_time = t
                    oldest_id = condition.id
    return (oldest_time, oldest_id)


class TraitsModel(object):
    def __init__(self, traits=None):
        """
        A list of genetically determined ``traits``,
        linked together by possibly shared effects of ``environments``
        and by ``fitness_functions`` that depend on their values.

        A TraitsModel always includes a multiplicative trait called "fitness". This is
        automatically included, so this should *not* be included in ``traits``.

        On initialization ``environments`` and ``fitness_functions``
        are empty and can be added with :meth:`.add_environment`
        and :meth:`.add_fitness_function`.
        These must all have unique IDs.

        If multiple environments are added that overlap in time or space, then
        their effects will combine additively.

        If multiple fitness functions are added that overlap in time or space,
        then their effects will combine multiplicatively.

        :ivar traits: List of :class:`Trait` objects, with unique IDs.
            or ``None``.
        :vartype traits: list
        """
        # we'll put fitness first, BUT DO NOT RELY ON THIS
        self.traits = [Trait(id="fitness", type="multiplicative")]
        if traits is not None:
            self.traits.extend(traits)

        pids = [p.id for p in self.traits]
        if len(set(pids)) != len(pids):
            raise ValueError("Trait IDs must be unique.")

        self.environments = []
        self.fitness_functions = []
        # We *could* take in Environment and FitnessFunction objects
        # to the construtor here, but no need;
        # we'll just add them with the add_X functions.

    def _check_traits_defined(self, trait_ids):
        tids = [t.id for t in self.traits]
        for tid in trait_ids:
            if tid not in tids:
                raise ValueError(f"Unknown trait ID `{tid}'.")

    def add_fitness_function(self, **kwargs):
        """
        Adds a :class:`.FitnessFunction` to the :class:`TraitsModel`.
        The arguments are passed directly to :class:`FitnessFunction`;
        see that documentation for more information.
        The IDs of all traits referred to by this fitness function
        must be present in the :class:`TraitsModel`.
        """
        ff = FitnessFunction(**kwargs)
        self._check_traits_defined(ff.trait_ids)
        fids = [f.id for f in self.fitness_functions]
        if ff.id in fids:
            raise ValueError(
                "FitnessFunction IDs must be unique; "
                f" fitness function with ID `{ff.id}` already exists."
            )
        self.fitness_functions.append(ff)

    def add_environment(self, **kwargs):
        """
        Add random "environmental" (i.e., non-genetic) effects to the specified
        traits. The arguments are passed directly to :class:`Environment`;
        see that documentation for more information.
        The IDs of all traits referred to by this environment
        must be present in the :class:`TraitsModel`.
        """
        env = Environment(**kwargs)
        self._check_traits_defined(env.trait_ids)
        eids = [e.id for e in self.environments]
        if env.id in eids:
            raise ValueError(
                "Environment IDs must be unique; "
                f" environment with ID `{env.id}` already exists."
            )
        self.environments.append(env)


#     def check_model_params(self, model, params):
#         # Check for consistency with a given demographic model.
#         # TODO: is this where we want to do this?
#         #
#         # sort elements by epoch start time
#         # (which is the backwards-in-time generation time)
#         # tuples should look like
#         # (epoch_start, [list of applicable pops], distribution params...)
#         # check populations and generation times make sense?
#         # distributions should propagate backwards in time
#         # so to apply the same distribution everywhere, use the tuple
#         # (0, "all", ...)
#         # also check that num_traits matches distribution params
#         pass
#
#     def print(model):
#         """
#         Prints how the TraitsModel maps on a demographic model.
#
#         Note: this is why environments and fitness functions have IDs:
#         both need to be uniquely specified so we can show things like
#         "this environment applies to these populations for this time period".
#         """
#         # TODO
#         pass


@attr.s(kw_only=True)
class Trait:
    """
    Represents a single trait, something that we measure or observe.
    This class defines how the underlying (latent) value,
    which is a sum of genetic value and environmental deviation(s),
    is "transformed" to the observed value (the phenotype).
    The `transform` is thus analogous to an inverse link function
    from generalized linear models.

    The ``type`` can be either "additive" or "multiplicative",
    and determines whether the per-site genetic effects are added together
    or multiplied to produce the genetic value
    (in other words, the genetic component of the latent value):
    for additive traits, the genetic value is the sum of :math:`s_i`,
    where :math:`s_i` is the effect at site :math:`i` (including effects
    of dominance); for multiplicative traits, the genetic value
    is the product of :math:`1+s_i`.

    The genetic effects themselves, including dominance, are specified through
    :class:`DistributionOfMutationEffects`,
    while environmental deviations are specified through
    :class:`Environment`s.

    Options for "transform" are:

    "identity": the observed value is equal to the latent value.

    "threshold" (parameters: t): the observed value is equal to 1 if the latent
        value is greater than t, and is equal to 0 otherwise.

    "liability" (parameters center, slope): the observed value is equal to 1
        with probability 1 / (1 + exp((x - center) * slope)) for the latent value x,
        and is equal to 0 otherwise.

    TODO: Add "exponential" transform to get log-normal traits?

    TODO: Add Poisson (count) traits?

    TODO: Add "logistic" traits?

    :ivar id: ID of the trait (think of this as the 'name').
    :vartype id: str
    :ivar type: Type of the trait ("additive" or "multiplicative").
    :vartype type: str
    :ivar transform: Type of transformation. (default: "identity")
    :vartype transform: str
    :ivar transform_args: A list of parameters given to the transformation.
    :vartype transform_args: list
    """

    id = attr.ib(type=str)
    type = attr.ib(type=str)
    transform = attr.ib(default="identity", type=str)
    transform_args = attr.ib(default=None, type=list, converter=_copy_converter)

    def __attrs_post_init__(self):
        if not (isinstance(self.id, str) and self.id != ""):
            raise ValueError("id must be a nonempty string")
        if self.transform_args is None:
            self.transform_args = []
        if not (
            isinstance(self.type, str) and self.type in ("multiplicative", "additive")
        ):
            raise ValueError(f"Unknown trait type '{self.type}'.")
        if not isinstance(self.transform_args, list):
            raise ValueError("transform_args must be a list")
        if self.transform == "identity":
            if len(self.transform_args) != 0:
                raise ValueError("identity transform takes no parameters.")
        elif self.transform == "threshold":
            if len(self.transform_args) != 1:
                raise ValueError(
                    "threshold transform requires one parameter (the threshold)"
                )
        elif self.transform == "liability":
            if len(self.transform_args) != 2:
                raise ValueError(
                    "liability transform requires two parameters (center and slope)"
                )
            if self.transform_args[1] <= 0:
                raise ValueError("slope for liability transform must be positive")
        else:
            raise ValueError(f"Transform '{self.transform}' unknown.")


@attr.s(kw_only=True)
class Environment:
    """
    Represents random "environmental" (i.e., non-genetic) effects on traits.
    These are added to genetic values to produce the latent values
    for a :class:`Trait`.
    The Environment may be restricted to apply only to a given
    span of time and/or set of populations.

    Each environment has an ``id``; this is for debugging purposes, and so each
    environment used in the same :class:`.TraitsModel` should have a unique name.

    TODO: does this accept all of the distribution types? Right now it shares code
    with MutationType so it does.

    :ivar id: An ID (i.e., a name) for this environment.
    :vartype id: str
    :ivar trait_ids: List of trait IDs.
    :vartype trait_ids: list
    :ivar distribution_type: A str abbreviation for the distribution
        of environmental efffects (see TODO WHERE).
    :vartype distribution_type: str
    :ivar distribution_args: Arguments to the distribution.
    :vartype distribution_args: list
    :ivar time_intervals: List of tuples defining when (backward-in-time) this
        environment applies. Setting an upper limit of float('inf') will cause
        this environment to apply from the beginning of the simulation
        (including the burn-in). Units are generations.
        Defaults to applying for all of time.
    :vartype time_intervals: list
    :ivar population_list: List of population ids specifying the populations this
        environment applies to. Defaults to applying to all populations. These
        can be specified either as integers (representing population indices)
        or strings with the names of populations.
    :vartype population_list: list
    """

    id = attr.ib(type=str)
    trait_ids = attr.ib(type=list, converter=_copy_converter)  # list of trait IDs
    distribution_type = attr.ib(type=str)
    distribution_args = attr.ib(type=list, converter=_copy_converter)
    time_intervals = attr.ib(default=None, converter=_copy_converter)
    population_list = attr.ib(default=None, converter=_copy_converter)

    def __attrs_post_init__(self):
        if not (isinstance(self.id, str) and self.id != ""):
            raise ValueError("id must be a nonempty string")
        _check_trait_ids(self.trait_ids)
        _check_distribution(
            self.distribution_type, self.distribution_args, len(self.trait_ids)
        )
        if self.population_list is not None:
            for pid in self.population_list:
                if not isinstance(pid, (int, str)):
                    raise ValueError(
                        "population_list entries must be integers "
                        "representing population indices or must be "
                        "strings with population names."
                    )
            if len(self.population_list) != len(set(self.population_list)):
                raise ValueError("population_list contains repeated entries")

        if self.time_intervals is not None:
            _check_nonoverlapping_intervals(self.time_intervals)


@attr.s(kw_only=True)
class FitnessFunction:
    """
    A function that computes a component of fitness:
    the total fitness is obtained by multiplying together
    all fitness functions in the :class:`TraitsModel`.
    Each fitness function operates on a collection of traits
    (the ``trait_ids``), and returns a value that multiplies the fitness.
    The Fitness Function may be restricted to apply only to a given
    span of time and/or set of populations.

    Options for ``function_type``, and corresponding ``function_args``, are:

    "gaussian", arguments (m, s): fitness is given by the Gaussian density
        with mean m and (co)variance s; so if there is a single trait in ``trait_id``,
        then this is f(x) = exp((x - m)**2 / s) / sqrt(2*pi*s),
        while if there is more than one trait then it is the multivariate
        Gaussian density.

    "threshold", arguments (q, a, b): fitness is one of two values, either
        :math:`f(x) = a` if the quantile of :math:`x` among the values in the
        population is less than :math:`q`, and :math:`f(x) = b` otherwise.


    Each fitness function has an ``id``; this is for debugging purposes, and so each
    fitness function used in the same :class:`.TraitsModel` should have a unique name.

    :ivar id: An ID (i.e., a name) for this fitness function.
    :vartype id: str
    :ivar trait_ids: List of trait IDs.
    :vartype trait_ids: list
    :ivar function_type: String corresponding to fitness function type
    :vartype function_type: str
    :ivar function_args: Tuple containing parameters for the fitness function
    :vartype function_args: str
    :ivar time_intervals: List of tuples defining when (backward-in-time) this
        fitness function applies. Units are generations. Setting an upper limit
        of float('inf') will cause this fitness function to apply from the
        beginning of the simulation (including the burn-in).
        Defaults to applying for all of time.
    :vartype time_intervals: list
    :ivar population_list: List of population ids specifying the populations this
        fitness function applies to. Defaults to applying to all populations. These
        can be specified either as integers (representing population indices)
        or strings with the names of populations.

    :vartype population_list: list

    """

    id = attr.ib(type=str)
    trait_ids = attr.ib(type=list, converter=_copy_converter)
    function_type = attr.ib(type=str)
    function_args = attr.ib(type=tuple, converter=_copy_converter)
    time_intervals = attr.ib(default=None, converter=_copy_converter)
    population_list = attr.ib(default=None, converter=_copy_converter)

    def __attrs_post_init__(self):
        if not (isinstance(self.id, str) and self.id != ""):
            raise ValueError("id must be a nonempty string")
        _check_trait_ids(self.trait_ids)
        num_traits = len(self.trait_ids)
        if not isinstance(self.function_type, str):
            raise ValueError("function_type must be a str")
        _check_args_list(
            "function_args",
            self.function_args,
            num_traits,
            arrays=(self.function_type == "gaussian"),
        )

        if self.population_list is not None:
            for pid in self.population_list:
                if not isinstance(pid, (int, str)):
                    raise ValueError(
                        "population_list entries must be integers "
                        "representing population indices or must be "
                        "strings with population names."
                    )
            if len(self.population_list) != len(set(self.population_list)):
                raise ValueError("population_list contains repeated entries")

        if self.time_intervals is not None:
            _check_nonoverlapping_intervals(self.time_intervals)

        if self.function_type == "gaussian":
            _check_gaussian_args(self.function_args, num_traits)
        elif self.function_type == "threshold":
            if len(self.function_args) != 3:
                raise ValueError(
                    "threshold function takes three arguments: "
                    "(quantile, low_fitness, high_fitness)"
                )
            if self.function_args[0] < 0 or self.function_args[0] > 1:
                raise ValueError(
                    "quantile argument to threshold function "
                    "must be between 0 and 1."
                )
            if self.function_args[1] < 0 or self.function_args[2] < 0:
                raise ValueError(
                    "fitness arguments to threshold function must be nonnegative"
                )
        else:
            raise ValueError(f"Unknown function type {self.function_type}.")


@attr.s(kw_only=True)
class MutationType(object):
    """
    Class representing a "type" of mutation, that affects fitness and/or
    a collection of other traits.

    The main thing that mutation types carry is a way of drawing an *effect*
    for each new mutation from a distribution. This ``distribution_type`` should
    be one of:

    - ``f``: fixed, one parameter per trait (a single value for each)
    - ``e``: exponential, one parameter (mean)
    - ``g``: gamma, two parameters (mean, shape)
    - ``n``: normal, two parameters (mean, sd)
    - ``w``: Weibull, two parameters (scale, shape)
    - ``u``: Uniform, two parameters (min, max)
    - ``lp``: positive logNormal, two parameters (mean and sd on log scale; see rlnorm)
    - ``ln``: negative logNormal, two parameters (mean and sd on log scale; see rlnorm)
    - ``mvn``: TODO

    Currently, only "fixed" and "mvn" can apply to more than one trait.

    Type "lp" is always positive, and type "ln" is always negative: both use
    the same log-normal distribution, but "ln" is multiplied by -1.  For
    exponential and gamma, a negative mean can be provided, obtaining always
    negative values.

    Instead of a single dominance coefficient (which would be specified with
    `dominance_coeff`), a discretized relationship between dominance and
    effect can be implemented: if dominance_coeff_list is
    provided, mutations with effect ``s`` for which
    ``dominance_coeff_breaks[k-1] <= s <= dominance_coeff_breaks[k]`` will have
    ``dominance coefficient dominance_coeff[k]``. In other words, the first entry
    of ``dominance_coeff_list`` applies to any mutations with effect
    below the first entry of ``dominance_coeff_breaks``; the second entry of
    ``dominance_coeff_list`` applies to mutations with effect
    between the first and second entries of ``dominance_coeff_breaks``, and so
    forth. The list of breaks must therefore be of length one less than the
    list of dominance coefficients.

    :ivar trait_ids: A list of trait IDs this mutation type affects.
        (default: ["fitness"])
    :vartype trait_ids: list
    :ivar distribution_type: A str abbreviation for the distribution of
        effects that each new mutation of this type draws from (see above).
    :vartype distribution_type: str
    :ivar distribution_args: Arguments for the distribution type.
    :vartype distribution_type: list
    :ivar dominance_coeff: The dominance coefficient (negative = underdominance,
        0 = recessive, 0.5 = additive, 1.0 = completely dominant, > 1.0 = overdominant)
        Default: 0.5.
    :vartype dominance_coeff: float
    :ivar convert_to_substitution: Whether to retain any fixed mutations in the
        simulation: if not, we cannot ask about their frequency once fixed.
        (Either way, they will remain in the tree sequence).  Default: True.
    :vartype convert_to_substitution: bool
    :ivar dominance_coeff_list: Either None (the default) or a list of floats describing
        a list of dominance coefficients, to apply to different effects
        (see details). Cannot be specified along with dominance_coeff.
    :vartype dominance_coeff_list: list of floats
    :ivar dominance_coeff_breaks: Either None (the default) or a list of floats
        describing the intervals of selection coefficient over which each of the entries
        of dominance_coeff_list applies (see details). Must be of length one shorter than
        dominance_coeff_list.
    :vartype dominance_coeff_breaks: list of floats
    """

    trait_ids = attr.ib(default=None, type=list, converter=_copy_converter)
    distribution_type = attr.ib(default="f", type=str)
    distribution_args = attr.ib(default=None, type=list, converter=_copy_converter)
    dominance_coeff = attr.ib(default=None, type=float)
    convert_to_substitution = attr.ib(default=True, type=bool)
    dominance_coeff_list = attr.ib(default=None, type=list, converter=_copy_converter)
    dominance_coeff_breaks = attr.ib(default=None, type=list, converter=_copy_converter)

    def __attrs_post_init__(self):
        if self.trait_ids is None:
            self.trait_ids = ["fitness"]
        _check_trait_ids(self.trait_ids)

        if self.distribution_args is None:
            self.distribution_args = [0 for _ in self.trait_ids]

        if self.dominance_coeff is None and self.dominance_coeff_list is None:
            self.dominance_coeff = 0.5

        if self.dominance_coeff is not None:
            if (self.dominance_coeff_list is not None) or (
                self.dominance_coeff_breaks is not None
            ):
                raise ValueError(
                    "Cannot specify both dominance_coeff and dominance_coeff_list."
                )
            if not isinstance(self.dominance_coeff, (float, int)):
                raise ValueError("dominance_coeff must be a number.")
            if not np.isfinite(self.dominance_coeff):
                raise ValueError(
                    f"Invalid dominance coefficient {self.dominance_coeff}."
                )

        if self.dominance_coeff_list is not None:
            if len(self.trait_ids) != 1 or self.trait_ids[0] != "fitness":
                raise ValueError(
                    "Cannot specify dominance_coeff_list for non-fitness traits."
                )
            # disallow the inefficient and annoying length-one case
            if len(self.dominance_coeff_list) < 2:
                raise ValueError("dominance_coeff_list must have at least 2 elements.")
            for h in self.dominance_coeff_list:
                if not isinstance(h, (float, int)):
                    raise ValueError("dominance_coeff_list must be a list of numbers.")
                if not np.isfinite(h):
                    raise ValueError(f"Invalid dominance coefficient {h}.")
            if self.dominance_coeff_breaks is None:
                raise ValueError(
                    "A list of dominance coefficients provided but no breaks."
                )
            if len(self.dominance_coeff_list) != len(self.dominance_coeff_breaks) + 1:
                raise ValueError(
                    "len(dominance_coeff_list) must be equal "
                    "to len(dominance_coeff_breaks) + 1"
                )
            lb = -1 * np.inf
            for b in self.dominance_coeff_breaks:
                if not isinstance(b, (float, int)):
                    raise ValueError(
                        "dominance_coeff_breaks must be a list of numbers."
                    )
                if not np.isfinite(b):
                    raise ValueError(f"Invalid dominance coefficient break {b}.")
                if b < lb:
                    raise ValueError("dominance_coeff_breaks must be nondecreasing.")
                lb = b

        if not isinstance(self.convert_to_substitution, bool):
            raise ValueError("convert_to_substitution must be bool.")

        _check_distribution(
            self.distribution_type, self.distribution_args, len(self.trait_ids)
        )

        # rewrite some of these for Eidos
        # TODO: this should probably happen downstream, in slim_engine.py?
        if self.distribution_type in ("lp", "ln"):
            # lognormal distribution:
            logmean, logsd = self.distribution_args
            sign = "" if self.distribution_type == "lp" else "-1 *"
            self.distribution_args = [
                f"return {sign}rlnorm(1, {logmean} + log(Q), {logsd});"
            ]
            self.distribution_type = "s"
        elif self.distribution_type == "u":
            umin, umax = self.distribution_args
            self.distribution_args = [f"return runif(1, Q * {umin}, Q * {umax});"]
            self.distribution_type = "s"

        # The index(s) of the param in the distribution_args list that should be
        # multiplied by Q when using --slim-scaling-factor Q.
        # Note that "u", "lp", and "ln" got remapped to "s" above,
        # which is why they do not appear here.
        scaling_factor_index_lookup = {
            "f": [0],  # fixed value
            "g": [0],  # mean
            "e": [0],  # mean
            "n": [0, 1],  # mean and sd
            "w": [0],  # scale
            "s": [],  # script types should just printout arguments
            "mvn": [],  # TODO: how to do scaling for multivariate traits
        }
        assert self.distribution_type in scaling_factor_index_lookup
        self.Q_scaled_index = scaling_factor_index_lookup[self.distribution_type]

    @property
    def is_neutral(self):
        """
        Tests whether the mutation type is strictly neutral. This is defined here to
        be:
        - only affecting "fitness";
        - of type "f";
        - and with fitness effect 0.0,
        and so excludes other situations that also produce only neutral
        mutations (e.g., exponential with mean 0, or affecting some other trait
        with no effect on fitness).

        TODO: make a TraitsModel method that looks at whether a trait affects
        fitness and so can decide whether additional mutation types are neutral.
        """
        neutral = (
            (len(self.trait_ids) == 1)
            and (self.trait_ids[0] == "fitness")
            and (self.distribution_type == "f")
            and (self.distribution_args[0] == 0)
        )
        return neutral


# at least conceptually a superclass of DFE, so we call it DME
@attr.s(kw_only=True)
class DistributionOfMutationEffects(object):
    """
    Class representing all mutations that affect a given segment of genome,
    and hence contains a list of :class:`.MutationType`
    and corresponding list of proportions,
    that gives the proportions of mutations falling in this region
    that are of the corresponding mutation type.

    ``proportions`` and ``mutation_types`` must be lists of the same length,
    and ``proportions`` should be nonnegative numbers summing to 1.

    :ivar ~.mutation_types: A list of :class:`.MutationType`
        objects associated with the DME. Defaults to an empty list.
    :vartype ~.mutation_types: list
    :ivar ~.proportions: A list of the proportions of new mutations that
        fall in to each of the mutation types (must sum to 1).
    :vartype ~.proportions: list
    :ivar ~.id: The unique identifier for this model. DME IDs should be
        short and memorable, and conform to the stdpopsim
        :ref:`naming conventions <sec_development_naming_conventions>`
        for DME models.
    :vartype ~.id: str
    :ivar ~.description: A short description of this model as it would be used in
        written text, e.g., "Lognormal DME". This should
        describe the DME itself and not contain author or year information.
    :vartype ~.description: str
    :ivar long_description: A concise, but detailed, summary of the DME model.
    :vartype long_description: str
    """

    # TODO: what about id and description and stuff???
    # TODO: implement __str__??
    mutation_types = attr.ib(default=None)
    proportions = attr.ib(default=None)

    def __attrs_post_init__(self):
        self.mutation_types = [] if self.mutation_types is None else self.mutation_types
        if self.proportions is None and len(self.mutation_types) == 0:
            self.proportions = []
        elif self.proportions is None:
            # will error below if this doesn't make sense
            self.proportions = [1]

        if not (isinstance(self.proportions, (collections.abc.Sequence, np.ndarray))):
            raise ValueError("proportions must be a list or numpy array.")

        if not (isinstance(self.mutation_types, list)):
            raise ValueError("mutation_types must be a list.")

        if not (len(self.proportions) == len(self.mutation_types)):
            raise ValueError(
                "proportions and mutation_types must be lists of the same length."
            )

        for p in self.proportions:
            if not isinstance(p, (float, int)) or p < 0:
                raise ValueError("proportions must be nonnegative numbers.")

        if len(self.proportions) > 0:
            sum_p = sum(self.proportions)
            if not np.isclose(sum_p, 1):
                raise ValueError("proportions must sum to 1.0.")

        for m in self.mutation_types:
            if not isinstance(m, MutationType):
                raise ValueError(
                    "mutation_types must be a list of MutationType objects."
                )

    @property
    def is_neutral(self):
        # TODO: implement me
        return False


@attr.s(kw_only=True)
class DFE(DistributionOfMutationEffects):
    """
    Class representing a "Distribution of Fitness Effects", i.e., a DFE.
    The class records the different *mutation types*, and the *proportions*
    with which they occur. The overall rate of mutations will be determined
    by the Contig to which the DFE is applied (see :meth:`.Contig.add_dme`).

    This is a specialization of :class:`.DistributionOfMutationEffects`
    to distributions that only affect fitness, and have associated publications
    (and hence citations).

    Instances of this class are constructed by DFE implementors, following the
    :ref:`developer documentation <sec_development_dfe_model>`. To instead
    obtain a pre-specified model as listed in the :ref:`sec_catalog`,
    see :meth:`Species.get_dfe`.

    ``proportions`` and ``mutation_types`` must be lists of the same length,
    and ``proportions`` should be nonnegative numbers summing to 1.

    :ivar ~.mutation_types: A list of :class:`.MutationType` objects associated
        with the DFE. Defaults to an empty list.
    :vartype ~.mutation_types: list
    :ivar ~.proportions: A list of the proportions of new mutations that
        fall in to each of the mutation types (must sum to 1).
    :vartype ~.proportions: list
    :ivar ~.id: The unique identifier for this model. DFE IDs should be
        short and memorable, and conform to the stdpopsim
        :ref:`naming conventions <sec_development_naming_conventions>`
        for DFE models.
    :vartype ~.id: str
    :ivar ~.description: A short description of this model as it would be used in
        written text, e.g., "Lognormal DFE". This should
        describe the DFE itself and not contain author or year information.
    :vartype ~.description: str
    :ivar long_description: A concise, but detailed, summary of the DFE model.
    :vartype long_description: str
    :ivar citations: A list of :class:`Citations <.Citation>`, that describe the primary
        reference(s) for the DFE model.
    :vartype citations: list of :class:`Citation`
    """

    id = attr.ib()
    description = attr.ib()
    long_description = attr.ib()
    citations = attr.ib(default=None)
    qc_dfe = attr.ib(default=None)

    def __attrs_post_init__(self):
        super().__attrs_post_init__()
        self.citations = [] if self.citations is None else self.citations

    @property
    def is_neutral(self):
        return all([m.is_neutral for m in self.mutation_types])

    def __str__(self):
        long_desc_lines = [
            line.strip()
            for line in textwrap.wrap(textwrap.dedent(self.long_description))
        ]
        long_desc = "\n║                     ".join(long_desc_lines)
        s = (
            "DFE:\n"
            f"║  id               = {self.id}\n"
            f"║  description      = {self.description}\n"
            f"║  long_description = {long_desc}\n"
            f"║  citations        = {[cite.doi for cite in self.citations]}\n"
        )
        return s

    def register_qc(self, qc_dfe):
        """
        Register a QC model implementation for this DFE.
        """
        if not isinstance(qc_dfe, self.__class__):
            raise ValueError(f"Cannot register non-DFE '{qc_dfe}' as QC DFE.")
        if self.qc_dfe is not None:
            raise ValueError(f"QC DFE already registered for {self.id}.")
        self.qc_dfe = qc_dfe


def neutral_dfe(convert_to_substitution=True):
    id = "neutral"
    description = "neutral DFE"
    long_description = "strictly neutral mutations"
    neutral = MutationType(convert_to_substitution=convert_to_substitution)
    return DFE(
        id=id,
        description=description,
        long_description=long_description,
        mutation_types=[neutral],
        proportions=[1.0],
    )


def _check_args_list(argname, args, dim, arrays):
    # TODO: do we really want to require arrays?
    if not isinstance(args, list):
        raise ValueError(f"{argname} must be list.")
    for i in range(len(args)):
        if not arrays:
            if not isinstance(args[i], (float, int)):
                raise ValueError(f"{argname}[{i}] is not a number.")
            if not np.isfinite(args[i]):
                raise ValueError(f"{argname}[{i}] is an invalid parameter.")
        else:
            if not isinstance(args[i], np.ndarray):
                raise ValueError(f"{argname}[{i}] is not a numpy array.")
            # TODO: check that entries are finite numbers


def _check_gaussian_args(args, dim):
    if len(args) != 2:
        raise ValueError(
            "Normal distribution requires two parameters: "
            "a mean (or mean vector) and a variance (or covariance matrix)."
        )
    if len(args[0].shape) != 1 or args[0].shape[0] != dim:
        raise ValueError(
            f"Multivariate normal mean vector must be 1 dimensional of length {dim}."
        )
    if len(args[1].shape) != 2:
        raise ValueError("Multivariate normal covariance matrix must be 2 dimensional.")
    if args[1].shape != (dim, dim):
        raise ValueError(
            "Multivariate normal covariance matrix must be square, "
            f"with dimensions ({dim}, {dim})."
        )
    if not np.allclose(args[1], args[1].T):
        raise ValueError("Multivariate normal covariance matrix must be symmetric.")
    try:
        np.linalg.cholesky(args[1])
    except np.linalg.LinAlgError as ve:
        raise ValueError(
            "Problem with multivariate normal covariance matrix: " + str(ve)
        )


def _check_distribution(distribution_type, distribution_args, dim):
    if not isinstance(distribution_type, str):
        raise ValueError("distribution_type must be str.")

    _check_args_list(
        "distribution_args", distribution_args, dim, arrays=(distribution_type == "mvn")
    )

    # shortcut, so we don't have to validate dim in all sub-cases below
    if (dim > 1) and (distribution_type not in ["f", "mvn"]):
        raise ValueError(
            f"Distribution type '{distribution_type}' is not "
            "implemented as a multivariate distribution."
        )

    # To add a new distribution type: validate the
    # distribution_args here, and add unit tests.
    if distribution_type == "f":
        # Fixed-value (non-random)
        if len(distribution_args) != dim:
            raise ValueError(
                "Fixed-value mutation type argument must be a list of "
                f"length {dim}, the number of traits."
            )
    elif distribution_type == "g":
        # Gamma distribution with (mean, shape)
        # parameterization. A negative value for the mean is permitted,
        # and indicates a reflection of the horizontal axis.
        # See Eidos documentation for rgamma().
        if len(distribution_args) != 2:
            raise ValueError(
                "Gamma distribution (distribution_type='g') "
                "uses a (mean, shape) parameterisation."
            )
        if distribution_args[1] <= 0:
            raise ValueError("The shape parameter must be positive.")
    elif distribution_type == "e":
        # An exponential distribution(mean).
        # See Eidos documentation for rexp().
        if len(distribution_args) != 1:
            raise ValueError(
                "Exponential distribution (distribution_type='e') "
                "uses a (mean) parameterisation."
            )
    elif distribution_type == "n":
        # A normal distribution (mean, standard deviation).
        # See Eidos documentation for rnorm().
        # TODO: combine with _check_gaussian_args?
        if len(distribution_args) != 2:
            raise ValueError(
                "Normal distribution (distribution_type='n') "
                "uses a (mean, sd) parameterisation."
            )
        if distribution_args[1] < 0:
            raise ValueError("The sd parameter must be nonnegative.")
    elif distribution_type == "w":
        # A Weibull-distributed fitness effect (scale, shape).
        # See Eidos documentation for rweibull().
        if len(distribution_args) != 2:
            raise ValueError(
                "Weibull distribution (distribution_type='w') "
                "uses a (scale, shape) parameterisation."
            )
        if distribution_args[0] <= 0:
            raise ValueError("The scale parameter must be positive.")
        if distribution_args[1] <= 0:
            raise ValueError("The shape parameter must be positive.")
    elif distribution_type in ("lp", "ln"):
        # A lognormal distribution (meanlog, sdlog),
        # either positive or negative.
        # See Eidos documentation for rlnorm().
        if len(distribution_args) != 2:
            raise ValueError(
                "Lognormal distribution (distribution_type='lp'/'ln') "
                "uses a (meanlog, sdlog) parameterisation, requiring sdlog > 0."
            )
        if distribution_args[1] < 0:
            raise ValueError("The sdlog parameter must be nonnegative.")
    elif distribution_type == "u":
        # Uniform
        if len(distribution_args) != 2 or distribution_args[0] > distribution_args[1]:
            raise ValueError(
                "Uniform distribution (distribution_type='u') "
                "uses a (min, max) parameterisation, with min <= max."
            )
    elif distribution_type == "mvn":
        # Multivariate Normal distribution with
        #   (mean, covariance, indices) parameterization.
        _check_gaussian_args(distribution_args, dim)
    else:
        raise ValueError(f"{distribution_type} is not a supported distribution type.")


def _resolve_population_ids(traits_model, demographic_model):
    # Rewrite each condition's population_list, in place, so that
    # populations are integer indices in the demographic model's
    # population list.
    pop_names = [pop.name for pop in demographic_model.model.populations]
    for event in traits_model.environments + traits_model.fitness_functions:
        if event.population_list is None:
            continue
        pop_id_list = []
        for population in event.population_list:
            if isinstance(population, int):
                if population >= len(pop_names) or population < 0:
                    raise ValueError("Population index out of bounds.")
                pop_id_list.append(population)
            else:
                try:
                    pop_id_list.append(pop_names.index(population))
                except ValueError:
                    raise ValueError(
                        "Population label supplied not in demographic model."
                    )
        if len(pop_id_list) != len(set(pop_id_list)):
            raise ValueError("Repeated population indices.")
        event.population_list = pop_id_list


def _collect_valid_population_intervals(demographic_model, debugger=None):
    # Return a dict mapping population id to a list of [start, end) time
    # intervals (in generations, backward in time) during which that
    # population is active. The oldest interval ends at infinity.
    if debugger is None:
        debugger = demographic_model.model.debug()
    valid_times = {p.id: [] for p in demographic_model.model.populations}
    for epoch in debugger.epochs:
        for pop in epoch.populations:
            if not pop.active:
                continue
            intervals = valid_times[pop.id]
            if intervals and intervals[-1][1] == epoch.start_time:
                intervals[-1][1] = epoch.end_time
            else:
                intervals.append([epoch.start_time, epoch.end_time])
    return valid_times


def _standardize_condition(condition, valid_intervals):
    # Split a condition into one copy per population, with time intervals
    # made explicit: intervals ending at infinity are clipped to the
    # times the population is active, and finite intervals are required
    # to fall entirely within such times. A condition without populations
    # applies to every population that is active during its times, so it
    # is expanded into per-population copies whose intervals are clipped
    # to each population's activity; populations with no overlap are
    # dropped rather than being an error.
    if condition.population_list is None:
        intervals = condition.time_intervals
        if intervals is None:
            intervals = [[0, float("inf")]]
        new_conditions = []
        for p, activity in valid_intervals.items():
            clipped = []
            for a, b in intervals:
                for c, d in activity:
                    lo, hi = max(a, c), min(b, d)
                    if lo < hi:
                        clipped.append([lo, hi])
            if clipped:
                p_copy = copy.deepcopy(condition)
                p_copy.population_list = [p]
                p_copy.time_intervals = clipped
                new_conditions.append(p_copy)
        return new_conditions
    new_conditions = []
    for p in condition.population_list:
        if p not in valid_intervals:
            raise ValueError("Population index out of bounds.")
        standardized = []
        if condition.time_intervals is None:
            standardized.extend(valid_intervals[p])
        else:
            for interval in sorted(condition.time_intervals):
                is_valid = False
                if interval[1] != float("inf"):
                    # If interval is finite, we must identify
                    # a model interval that contains it
                    for demo_interval in valid_intervals[p]:
                        if (
                            interval[0] >= demo_interval[0]
                            and interval[1] <= demo_interval[1]
                        ):
                            is_valid = True
                            standardized.append(interval)
                else:
                    # If interval is infinite, we must identify
                    # all overlapping model intervals
                    for demo_interval in valid_intervals[p]:
                        if interval[0] < demo_interval[0]:
                            standardized.append(demo_interval)
                            is_valid = True
                        elif (
                            interval[0] >= demo_interval[0]
                            and interval[0] < demo_interval[1]
                        ):
                            standardized.append([interval[0], demo_interval[1]])
                            is_valid = True
                if not is_valid:
                    raise ValueError(
                        "An environment or a fitness function was "
                        "specified for a population with a time interval "
                        "during which that population does not exist."
                    )
        p_copy = copy.deepcopy(condition)
        p_copy.population_list = [p]
        p_copy.time_intervals = standardized
        new_conditions.append(p_copy)

    return new_conditions


def _align_traits_model_demography(traits_model, demographic_model, debugger=None):
    # Return a copy of ``traits_model`` in which every environment and
    # fitness function has a single integer population id and explicit
    # time intervals consistent with the demographic model.
    traits_model = copy.deepcopy(traits_model)
    if not traits_model.environments and not traits_model.fitness_functions:
        return traits_model
    _resolve_population_ids(traits_model, demographic_model)
    valid_intervals = _collect_valid_population_intervals(demographic_model, debugger)
    new_env = []
    for env in traits_model.environments:
        new_env.extend(_standardize_condition(env, valid_intervals))
    traits_model.environments = new_env

    new_ff = []
    for ff in traits_model.fitness_functions:
        new_ff.extend(_standardize_condition(ff, valid_intervals))
    traits_model.fitness_functions = new_ff
    return traits_model


# The next two functions are copied from msprime/core.py (GPL-3), so
# that our tables render in the same style as msprime's
# DemographyDebugger.


def _text_table_row(data, alignments, widths):
    num_lines = max(len(item) for item in data)
    for item in data:
        assert isinstance(item, list)
        item.extend([""] * (num_lines - len(item)))
        assert len(item) == num_lines
    s = ""
    for line in range(num_lines):
        out_line = "│"
        for value, align, width in zip(data, alignments, widths):
            out_line += f"{value[line]:{align}{width - 1}}│"
        out_line += "\n"
        s += out_line
    return s


def _text_table(caption, column_titles, column_alignments, data):
    N = len(column_titles)
    assert len(column_alignments) == N
    widths = np.array([len(title) for title in column_titles], dtype=int)
    for row in data + [column_titles]:
        assert N == len(row)
        for j in range(N):
            widths[j] = max(widths[j], max((len(line) for line in row[j]), default=0))
    widths += 3

    hline = "─" * (sum(widths) - 1)
    out = f"{caption}\n"
    out += f"┌{hline}┐\n"
    out += f"{_text_table_row(column_titles, column_alignments, widths)}"
    out += f"├{hline}┤\n"
    for split_row in data:
        out += f"{_text_table_row(split_row, column_alignments, widths)}"
    out += f"└{hline}┘\n"
    return out


class TraitsDebugger:
    """
    Shows how a :class:`.TraitsModel` lines up with a demographic model,
    in the style of ``msprime.DemographyDebugger``. One table is printed
    per time interval within which both the demography and the set of
    applicable environments and fitness functions are constant, so
    epochs of the demographic model are split wherever an environment or
    a fitness function starts or ends.

    The environments and fitness functions shown are the standardized
    ones, i.e. after populations are resolved and time intervals are
    made consistent with the demographic model. This is how the
    simulation engine interprets them.

    :param demographic_model: A :class:`.DemographicModel`.
    :param traits_model: A :class:`.TraitsModel`, or None for a plain
        demography table.
    """

    def __init__(self, demographic_model, traits_model=None):
        self.demographic_model = demographic_model
        if traits_model is None:
            traits_model = TraitsModel()
        self._demography_debugger = demographic_model.model.debug()
        self.traits_model = _align_traits_model_demography(
            traits_model, demographic_model, self._demography_debugger
        )
        # After alignment, every condition has a single population and
        # explicit time intervals.
        for condition in self._conditions():
            assert len(condition.population_list) == 1
            assert condition.time_intervals is not None
        self._make_epochs()

    def _conditions(self):
        return self.traits_model.environments + self.traits_model.fitness_functions

    def _make_epochs(self):
        dd = self._demography_debugger
        breaks = {epoch.start_time for epoch in dd.epochs}
        for condition in self._conditions():
            for start, end in condition.time_intervals:
                breaks.add(start)
                if np.isfinite(end):
                    breaks.add(end)
        breaks = sorted(breaks)
        sizes = dd.population_size_trajectory(breaks)
        self._epochs = []
        for i, start in enumerate(breaks):
            end = breaks[i + 1] if i + 1 < len(breaks) else float("inf")
            parent = next(e for e in dd.epochs if e.start_time <= start < e.end_time)
            # Sizes at the parent's own boundaries come from the parent,
            # so that instantaneous size changes display as in msprime.
            start_sizes = {}
            end_sizes = {}
            for pop in parent.populations:
                if start == parent.start_time:
                    start_sizes[pop.id] = pop.start_size
                else:
                    start_sizes[pop.id] = sizes[i][pop.id]
                if end == parent.end_time:
                    end_sizes[pop.id] = pop.end_size
                else:
                    end_sizes[pop.id] = sizes[i + 1][pop.id]
            self._epochs.append(
                dict(
                    start=start,
                    end=end,
                    parent=parent,
                    start_sizes=start_sizes,
                    end_sizes=end_sizes,
                )
            )

    def _active_condition_ids(self, conditions, pop_id, start, end):
        ids = []
        for condition in conditions:
            if pop_id not in condition.population_list:
                continue
            for a, b in condition.time_intervals:
                if a <= start and end <= b:
                    ids.append(condition.id)
                    break
        return ids

    def _populations_text(self, epoch):
        parent = epoch["parent"]
        column_titles = [
            [""],
            ["start"],
            ["end"],
            ["growth_rate"],
            ["environments"],
            ["fitness_functions"],
        ]
        data = []
        for pop in parent.populations:
            if not pop.active:
                continue
            envs = self._active_condition_ids(
                self.traits_model.environments, pop.id, epoch["start"], epoch["end"]
            )
            ffs = self._active_condition_ids(
                self.traits_model.fitness_functions,
                pop.id,
                epoch["start"],
                epoch["end"],
            )
            data.append(
                [
                    [pop.name],
                    [f"{epoch['start_sizes'][pop.id]: .1f}"],
                    [f"{epoch['end_sizes'][pop.id]: .1f}"],
                    [f"{pop.growth_rate: .3g}"],
                    envs if envs else [""],
                    ffs if ffs else [""],
                ]
            )
        caption = (
            f"Populations (total={len(parent.populations)} "
            f"active={parent.num_active_populations})"
        )
        return _text_table(caption, column_titles, ">>><^^", data)

    def _burn_in_note(self):
        oldest_time, oldest_id = _find_oldest_event(
            self.traits_model, self._demography_debugger
        )
        if oldest_id == "":
            return oldest_id
        oldest_event = max(e.start_time for e in self._demography_debugger.epochs)
        return (
            f"Burn-in ends {oldest_time:.3g} generations ago, set by "
            f"'{oldest_id}' (the oldest demographic event is "
            f"{oldest_event:.3g} generations ago).\n"
        )

    def print_history(self, output=None):
        """
        Prints the decorated demography table to the given file object,
        or to stdout if none is given.
        """
        print(self, file=output, end="")

    def __str__(self):
        # The box-drawing glue below is copied from
        # msprime.DemographyDebugger.__str__ (GPL-3).
        def indent(table, header_char="╟", depth=4):
            lines = table.splitlines()
            s = header_char + (" " * depth) + lines[0] + "\n"
            for line in lines[1:]:
                s += "║" + (" " * depth) + line + "\n"
            return s

        def box(title):
            N = len(title) + 2
            top = "╠" + ("═" * N) + "╗"
            bottom = "╠" + ("═" * N) + "╝"
            return f"{top}\n║ {title} ║\n{bottom}\n"

        out = "TraitsDebugger\n"
        out += self._burn_in_note()
        for i, epoch in enumerate(self._epochs):
            parent = epoch["parent"]
            if epoch["start"] > 0 and epoch["start"] == parent.start_time:
                title = f"Events @ generation {epoch['start']:.3g}"
                out += indent(
                    self.demographic_model.model._events_text(parent.events, title)
                )
            title = (
                f"Epoch[{i}]: [{epoch['start']:.3g}, {epoch['end']:.3g}) " "generations"
            )
            if np.isinf(epoch["end"]):
                title += " (includes burn-in)"
            out += box(title)
            out += indent(self._populations_text(epoch))
        return out
