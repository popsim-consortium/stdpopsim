"""
Methods related to traits and effects of mutations on them,
including fitness (so, this includes DFE machinery).
"""

import textwrap
import attr
import collections.abc
import numpy as np


def _copy_converter(x):
    if isinstance(x, list):
        x = x.copy()
    return x


class TraitsModel(object):
    def __init__(self, traits):
        """
        A collection of genetically determined traits,
        which are a list of ``traits``,
        linked together by possibly shared effects of ``environments``
        and by ``fitness_functions`` that depend on their values.

        On initialization ``environments`` and ``fitness_functions``
        are empty and can be added with :meth:`.add_environment`
        and :meth:`.add_fitness_function`.

        :ivar traits: List of :class:`Trait` objects.
        :vartype traits: list
        """
        pids = [p.id for p in traits]
        if len(set(pids)) != len(pids):
            raise ValueError("Trait IDs must be unique.")

        self.traits = traits
        self.environments = []
        self.fitness_functions = []
        # We *could* take in Environment and FitnessFunction objects
        # to the construtor here, but no need;
        # we'll just add them with the add_X functions.

    def add_fitness_function(
        self, traits, distribution_type, distribution_args, spatiotemporal=None
    ):
        new_fitness = FitnessFunction(
            traits, distribution_type, distribution_args, spacetime=spatiotemporal
        )
        self._fitness_functions.append(new_fitness)

    def add_environment(self, *, trait_ids, distribution_type, distribution_args):
        """
        Add random "environmental" (i.e., non-genetic) effects to the specified
        ``traits``.  See :class:`Environment` more more detail.
        """
        pids = [p.id for p in self.traits]
        for pid in trait_ids:
            if pid not in pids:
                raise ValueError(f"trait {pid} not in traits.")
        env = Environment(
            trait_ids=trait_ids,
            distribution_type=distribution_type,
            distribution_args=distribution_args,
        )
        self.environments.append(env)

    def check_model_params(self, model, params):
        # Check for consistency with a given demographic model.
        # TODO: is this where we want to do this?
        #
        # sort elements by epoch start time
        # (which is the backwards-in-time generation time)
        # tuples should look like
        # (epoch_start, [list of applicable pops], distribution params...)
        # check populations and generation times make sense?
        # distributions should propagate backwards in time
        # so to apply the same distribution everywhere, use the tuple
        # (0, "all", ...)
        # also check that num_traits matches distribution params
        pass


@attr.s(kw_only=True)
class Trait:
    """
    Represents a single trait, i.e., something that can be measured.
    This only defines how the underlying (latent) value,
    which is a sum of genetic value and environmental deviation,
    is mapped to the observed value.

    Options for "transform" (link function) are:

    "identity": the trait is equal to the latent value.

    "threshold" (parameters: x): the trait is equal to 1 if the latent
        value is less than x, and is equal to 0 otherwise.

    "liability" (parameters center, slope): the trait whose
        latent value is z is equal to 1 with probability
        1 / (1 + exp((x - center) * slope)), and is equal to 0 otherwise.

    TODO: Add "exponential" transform to get log-normal traits?

    TODO: Add Poisson (count) traits?

    TODO: Add "logistic" traits?

    :ivar id: ID of the trait (think of this as the 'name').
    :vartype id: str
    :ivar transform: Type of transformation.
    :vartype transform: str
    :ivar transform_args: A list of parameters given to the transformation.
    :vartype transform_args: list
    """

    id = attr.ib()
    transform = attr.ib(default="identity", type=str)
    transform_args = attr.ib(default=None, type=list, converter=_copy_converter)

    def __attrs_post_init__(self):
        if self.transform_args is None:
            self.transform_args = []
        if not isinstance(self.transform, str):
            raise ValueError("transform must be a str")
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
    These are all added to genetic values and may depend on time and/or population.

    TODO: should this be public?

    :ivar trait_ids: List of trait IDs.
    :vartype trait_ids: list
    :ivar distribution_type: A str abbreviation for the distribution
        of environmental efffects (see TODO WHERE).
    :vartype distribution_type: str
    :ivar distribution_type: A str abbreviation for the distribution
        of environmental efffects (see TODO WHERE).
    :vartype distribution_type: str
    """

    trait_ids = attr.ib(type=list, converter=_copy_converter)  # list of trait IDs
    distribution_type = attr.ib(type=str)
    distribution_args = attr.ib(type=list)
    # TODO: add later
    # start_time = attr.ib(default=None)
    # end_time = attr.ib(default=None)
    # populations = attr.ib(default=None)

    def __attrs_post_init__(self):
        dim = len(self.trait_ids)
        if dim < 1:
            raise ValueError("Must have at least one trait.")
        _check_distribution(self.distribution_type, self.distribution_args, dim)


@attr.s(kw_only=True)
class FitnessFunction:
    """
    TODO WRITE THIS BETTER
    Class to store a model of a component of fitness:
    each such component maps a collection of traits
    to a value that multiplies the fitness.
    Also contained here is when and where this component applies.

    Options for ``function_type``, and corresponding ``function_args``, are:

    TODO

    :ivar trait_ids: List of trait IDs.
    :vartype trait_ids: list
    :ivar function_type: String corresponding to fitness function type
    :vartype function_type: str
    :ivar function_args: Tuple containing parameters for the fitness function
    :vartype function_args: str
    """

    # :ivar spacetime: Generations and populations
    #    for which this fitness function applies
    # :vartype spacetime: list of tuples (?)

    # TODO check function_args depending on function_type
    # TODO check dimensions of traits against dimensions of function_args,
    # depending on function_type - plus check dimensions >=1
    # TODO much later - check spacetime is formatted correctly

    trait_ids = attr.ib(type=list)
    function_type = attr.ib(type=str)
    function_args = attr.ib(type=tuple)
    # spacetime = attr.ib(type=list)

    def __attrs_post_init__(self):
        if len(self.trait_ids) < 1:
            raise ValueError("At least one trait must be specified.")
        # _check_function_types(self.function_type, self.function_args,
        # len(self.trait_ids))


def _check_distribution(distribution_type, distribution_args, dim):
    if not isinstance(distribution_type, str):
        raise ValueError("distribution_type must be str.")

    if not isinstance(distribution_args, list):
        raise ValueError("distribution_args must be list.")

    for i in range(len(distribution_args)):
        if not isinstance(distribution_args[i], (float, int)):
            raise ValueError(f"distribution_args[{i}] is not a number.")
        if not np.isfinite(distribution_args[i]):
            raise ValueError(f"distribution_args[{i}] is an invalid parameter.")

    # shortcut, so we don't have to validate dim in all sub-cases below
    if (dim > 1) and (distribution_type not in ["f", "mvn"]):
        raise ValueError(
            f"Distribution type '{distribution_type}' is not "
            " implemented as a multivariate distribution."
        )

    # To add a new distribution type: validate the
    # distribution_args here, and add unit tests.
    if distribution_type == "f":
        # Fixed-value (non-random)
        if len(distribution_args) != dim:
            raise ValueError(
                "Fixed-value mutation type argument must be a list of "
                "length equal to number of traits."
            )
    elif distribution_type == "g":
        # Gamma distribution with (mean, shape)
        # parameterization. A negative value for the mean is permitted,
        # and indicates a reflection of the horizontal axis.
        # See Eidos documentation for rgamma().
        if len(distribution_args) != 2:
            raise ValueError(
                "Gamma-distributed sel. coefs. (distribution_type='g') "
                "use a (mean, shape) parameterisation."
            )
        if distribution_args[1] <= 0:
            raise ValueError("The shape parameter must be positive.")
    elif distribution_type == "e":
        # An exponential distribution(mean).
        # See Eidos documentation for rexp().
        if len(distribution_args) != 1:
            raise ValueError(
                "Exponentially-distributed sel. coefs. (distribution_type='e') "
                "use a (mean) parameterisation."
            )
    elif distribution_type == "n":
        # A normal distribution (mean, standard deviation).
        # See Eidos documentation for rnorm().
        if len(distribution_args) != 2:
            raise ValueError(
                "Normally-distributed sel. coefs. (distribution_type='n') "
                "use a (mean, sd) parameterisation."
            )
        if distribution_args[1] < 0:
            raise ValueError("The sd parameter must be nonnegative.")
    elif distribution_type == "w":
        # A Weibull-distributed fitness effect (scale, shape).
        # See Eidos documentation for rweibull().
        if len(distribution_args) != 2:
            raise ValueError(
                "Weibull-distributed sel. coef. (distribution_type='w') "
                "use a (scale, shape) parameterisation."
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
                "Lognormally-distributed sel. coefs. (distribution_type='lp'/'ln') "
                "use a (meanlog, sdlog) parameterisation, requiring sdlog > 0."
            )
        if distribution_args[1] < 0:
            raise ValueError("The sdlog parameter must be nonnegative.")
    elif distribution_type == "u":
        # Uniform
        if len(distribution_args) != 2 or distribution_args[0] > distribution_args[1]:
            raise ValueError(
                "Uniformly-distributed sel. coefs. (distribution_type='u') "
                "use a (min, max) parameterisation, with min <= max."
            )
    elif distribution_type == "mvn":
        # Multivariate Normal distribution with
        #   (mean, covariance, indices) parameterization.
        if len(distribution_args) != 2:
            raise ValueError(
                "multivariate normal requires two parameters "
                "in distribution_args: "
                "a mean vector and a covariance matrix."
            )
        if not isinstance(distribution_args[0], np.ndarray):
            raise ValueError(
                "mvn mean vector must be a numpy array" f"of length equal to {dim}."
            )
        if not isinstance(distribution_args[1], np.ndarray):
            raise ValueError("mvn covariance matrix must be specified as numpy array.")
        if len(distribution_args[0].shape) != 1 or distribution_args[0].shape[0] != dim:
            raise ValueError(
                "mvn mean vector must be 1 dimensional " f"of length {dim}."
            )
        if len(distribution_args[1].shape) != 2:
            raise ValueError("mvn covariance matrix must be 2 dimensional.")
        if distribution_args[1].shape != (dim, dim):
            raise ValueError(
                "mvn covariance matrix must be square, "
                f"with dimensions ({dim}, {dim})."
            )
        if not np.allclose(distribution_args[1], distribution_args[1].T):
            raise ValueError("mvn covariance matrix must be symmetric.")
        try:
            np.linalg.cholesky(distribution_args[1])
        except np.LinAlgError:
            raise ValueError("mvn covariance matrix is not positive definite.")
    else:
        raise ValueError(f"{distribution_type} is not a supported distribution type.")


@attr.s(kw_only=True)
class MutationType(object):
    """
    Class representing a "type" of mutation, allowing the mutation to affect
    fitness and/or trait(s). This design closely mirrors :class: MutationType.

    The main thing that mutation types carry is a way of drawing a selection
    coefficient for each new mutation. This ``distribution_type`` should be one
    of (see the SLiM manual for more information on these):

    - ``f``: fixed, one parameter (the selection coefficient)
    - ``e``: exponential, one parameter (mean)
    - ``g``: gamma, two parameters (mean, shape)
    - ``n``: normal, two parameters (mean, SD)
    - ``w``: Weibull, two parameters (scale, shape)
    - ``u``: Uniform, two parameters (min, max)
    - ``lp``: positive logNormal, two parameters (mean and sd on log scale; see rlnorm)
    - ``ln``: negative logNormal, two parameters (mean and sd on log scale; see rlnorm)
    - ``mvn``: TODO

    Type "lp" is always positive, and type "ln" is always negative: both use
    the same log-normal distribution, but "ln" is multiplied by -1.  For
    exponential and gamma, a negative mean can be provided, obtaining always
    negative values.

    Instead of a single dominance coefficient (which would be specified with
    `dominance_coeff`), a discretized relationship between dominance and
    selection coefficient can be implemented: if dominance_coeff_list is
    provided, mutations with selection coefficient s for which
    dominance_coeff_breaks[k-1] <= s <= dominance_coeff_breaks[k] will have
    dominance coefficient dominance_coeff[k]. In other words, the first entry
    of dominance_coeff_list applies to any mutations with selection coefficient
    below the first entry of dominance_coeff_breaks; the second entry of
    dominance_coeff_list applies to mutations with selection coefficient
    between the first and second entries of dominance_coeff_breaks, and so
    forth. The list of breaks must therefore be of length one less than the
    list of dominance coefficients.


    TODO: is "dominance_coeff_list" still the way we want to do things?
    SLiM is more flexible in this now.

    :ivar trait_ids: A list of IDs of traits this mutation type affects.
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
        a list of dominance coefficients, to apply to different selection coefficients
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

        if not isinstance(self.trait_ids, list):
            # Note: it'd be nice to also accept tuples, but note we can't
            # just check for being a collections.abc.Sequence since
            # then `trait_id = "fitness"` would pass and be interpreted
            # as seven trait IDs, named "f", "i", "t", etcetera.
            # So, just require a list.
            raise ValueError("Trait IDs must be a list.")
        for pid in self.trait_ids:
            if not (isinstance(pid, str) and (len(pid) > 0)):
                raise ValueError(
                    "Each trait ID must be a nonempty string; " f"found {pid}."
                )
        if (len(self.trait_ids) == 0) or (
            len(set(self.trait_ids)) != len(self.trait_ids)
        ):
            raise ValueError("Trait IDs must be a nonempty list of unique strings.")

        if self.distribution_args is None:
            self.distribution_args = [0 for _ in self.trait_ids]

        if self.dominance_coeff is None and self.dominance_coeff_list is None:
            self.dominance_coeff = 0.5

        if self.dominance_coeff is not None:
            if (self.dominance_coeff_list is not None) or (
                self.dominance_coeff_breaks is not None
            ):
                raise ValueError(
                    "Cannot specify both dominance_coeff " "and dominance_coeff_list."
                )
            if not isinstance(self.dominance_coeff, (float, int)):
                raise ValueError("dominance_coeff must be a number.")
            if not np.isfinite(self.dominance_coeff):
                raise ValueError(
                    f"Invalid dominance coefficient {self.dominance_coeff}."
                )

        if self.dominance_coeff_list is not None:
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
        objects associated with the DFE. Defaults to an empty list.
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
    """

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


@attr.s(kw_only=True)
class DFE(DistributionOfMutationEffects):
    """
    Class representing a "Distribution of Fitness Effects", i.e., a DFE.
    The class records the different *mutation types*, and the *proportions*
    with which they occur. The overall rate of mutations will be determined
    by the Contig to which the DFE is applied (see :meth:`.Contig.add_dfe`).

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
