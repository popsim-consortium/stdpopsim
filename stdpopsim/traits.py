"""
TODO: docstring goes here (preferably, descriptive)
"""

import attr
import numpy as np


# TODO: somewhere we need to check that the number of traits
# is consistent between Phenotypes, Environments, FitnessFunctions,
# MultivariateMutationTypes, DistributionofMutationEffects (DMEs), etc.

# TODO: this was moved from dfe.py so technically
# this would break backward compatibility, but
# I believe this method is no longer used
# anywhere in dfe.py, and it's private so this should be chill
# PLR: thanks for thinking about this, but moving/changing/deleting
# private methods doesn't break backwards compatibility; they're private
def _copy_converter(x):
    if isinstance(x, list):
        x = x.copy()
    return x


class Traits(object):
    def __init__(self, phenotypes):
        """
        A collection of genetically determined traits,
        which are a list of ``phenotypes``,
        linked together by possibly shared effects of ``environments``
        and by ``fitness_functions`` that depend on their values.

        On initialization ``environments`` and ``fitness_functions``
        are empty and can be added with :meth:`.add_environment`
        and :meth:`.add_fitness_function`.

        :ivar phenotypes: List of :class:`Phenotype` objects.
        :vartype phenotypes: list
        """
        pids = [p.id for p in phenotypes]
        if len(set(pids)) != len(pids):
            raise ValueError("Phenotype IDs must be unique.")

        self.phenotypes = phenotypes
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

    def add_environment(self, *, phenotype_ids, distribution_type, distribution_args):
        """
        Add random "environmental" (i.e., non-genetic) effects to the specified
        ``phenotypes``.  See :class:`Environment` more more detail.
        """
        pids = [p.id for p in self.phenotypes]
        for pid in phenotype_ids:
            if pid not in pids:
                raise ValueError(f"Phenotype {pid} not in phenotypes.")
        env = Environment(
            phenotype_ids=phenotype_ids,
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
class Phenotype:
    """
    Represents a single phenotype, i.e., something that can be measured.
    This only defines how the underlying (latent) value,
    which is a sum of genetic value and environmental deviation,
    is mapped to the observed value.

    Options for "transform" (link function) are:

    "identity": the phenotype is equal to the latent value.

    "threshold" (parameters: x): the phenotype is equal to 1 if the latent
        value is less than x, and is equal to 0 otherwise.

    "liability" (parameters center, slope): the phenotype whose
        latent value is z is equal to 1 with probability
        1 / (1 + exp((x - center) * slope)), and is equal to 0 otherwise.

    TODO: Add "exponential" transform to get log-normal phenotypes?

    TODO: Add Poisson (count) phenotypes?

    :ivar id: ID of the trait (think of this as the 'name').
    :vartype id: str
    :ivar transform: Type of transformation.
    :vartype transform: str
    :ivar params: Parameters given to the transformation.
    :vartype params: tuple
    """

    id = attr.ib()
    transform = attr.ib(default="identity")
    params = attr.ib(default=())

    def __attrs_post_init__(self):
        if self.transform == "identity":
            if len(self.params) != 0:
                raise ValueError("identity transform takes no parameters.")
        elif self.transform == "threshold":
            if len(self.params) != 1:
                raise ValueError(
                    "threshold transform requires one parameter (the threshold)"
                )
        elif self.transform == "liability":
            if len(self.params) != 2:
                raise ValueError(
                    "threshold transform requires two parameters " "(center and slope)"
                )
        else:
            raise ValueError(f"Transform {self.transform} unknown.")


@attr.s(kw_only=True)
class Environment:
    """
    Represents random "environmental" (i.e., non-genetic) effects on traits.
    These are all added to genetic values and may depend on time and/or population.

    TODO: should this be public?

    :ivar phenotype_ids: List of phenotype IDs.
    :vartype phenotype_ids: list
    :ivar distribution_type: A str abbreviation for the distribution
        of environmental efffects (see TODO WHERE).
    :vartype distribution_type: str
    :ivar distribution_type: A str abbreviation for the distribution
        of environmental efffects (see TODO WHERE).
    :vartype distribution_type: str
    """

    phenotype_ids = attr.ib()  # list of phenotype IDs
    distribution_type = attr.ib()
    distribution_args = attr.ib()
    # TODO: add later
    # start_time = attr.ib(default=None)
    # end_time = attr.ib(default=None)
    # populations = attr.ib(default=None)

    def __attrs_post_init__(self):
        dim = len(self.phenotype_ids)
        if dim < 1:
            raise ValueError("Must have at least one phenotype.")
        elif dim == 1:
            _check_univariate_distribution(
                self.distribution_type, self.distribution_args
            )
        else:
            _check_multivariate_distribution(
                self.distribution_type, self.distribution_args, dim
            )


class FitnessFunction(object):
    """
    Class to store a fitness function.

    :ivar traits: List of trait names or indices, as initialized in Traits object.
    :vartype traits: list
    :ivar function_type: One-letter string corresponding to fitness function type
    :vartype function_type: str
    :ivar function_args: Tuple containing parameters for the fitness function
    :vartype function_args: str
    :ivar spacetime: Generations and populations for which this fitness function applies
    :vartype spacetime: list of tuples (?)
    """

    # TODO check function_args depending on function_type
    # TODO check dimensions of traits against dimensions of function_args,
    # depending on function_type - plus check dimensions >=1
    # TODO much later - check spacetime is formatted correctly

    supported_function_types = []  # TODO

    traits = attr.ib(type=list)
    function_type = attr.ib(type=str)
    function_args = attr.ib(type=tuple)
    # spacetime = attr.ib(type=list)

    def __attrs_post_init__(self):
        if len(self.traits) < 1:
            raise ValueError("At least one trait must be specified.")
        if self.function_type not in self.supported_function_types:
            raise ValueError("Proposed fitness function not supported at this time.")


@attr.s(kw_only=True)
class MultivariateEffectSizeDistribution:
    def __attrs_post_init__(self):
        # Make sure fitness is not affected
        if 0 in self.supported_indices:
            raise ValueError(
                "distributions cannot simultaneously directly affect fitness "
                "and traits. that is, if 0 (fitness) is a supported index "
                "no other index can also be supported in supported_indices."
            )


def _check_multivariate_distribution(distribution_type, distribution_args, dim):
    if not isinstance(distribution_type, str):
        raise ValueError("distribution_type must be str.")

    if not isinstance(distribution_args, list):
        raise ValueError("distribution_args must be list.")

    if distribution_type == "mvn":
        _check_multivariate_normal_args(distribution_args, dim)
    else:
        raise ValueError(f"{distribution_type} is not a supported distribution type.")


def _check_multivariate_normal_args(distribution_args, dim):
    # TODO: I don't think we need the "list of dimensions that are not zero"
    # Multivariate Normal distribution with
    #   (mean, covariance, indices) parameterization.
    if len(distribution_args) != 3:
        raise ValueError(
            "multivariate normal requires three parameters "
            "in distribution_args: "
            "a mean vector, covariance matrix, and list of dimensions "
            "that are not zero."
        )
    if not isinstance(distribution_args[0], np.ndarray):
        raise ValueError(
            "mvn mean vector must be a numpy array" f"of length equal to {dim}."
        )
    if not isinstance(distribution_args[1], np.ndarray):
        raise ValueError("mvn covariance matrix must be specified as numpy array.")
    if len(distribution_args[0].shape) != 1 or distribution_args[0].shape[0] != dim:
        raise ValueError("mvn mean vector must be 1 dimensional " f"of length {dim}.")
    if len(distribution_args[1].shape) != 2:
        raise ValueError("mvn covariance matrix must be 2 dimensional.")
    if distribution_args[1].shape != (dim, dim):
        raise ValueError(
            "mvn covariance matrix must be square, " f"with dimensions ({dim}, {dim})."
        )
    if not np.allclose(distribution_args[1], distribution_args[1].T):
        raise ValueError("mvn covariance matrix must be symmetric.")
    try:
        np.linalg.cholesky(distribution_args[1])
    except np.LinAlgError:
        raise ValueError("mvn covariance matrix is not positive definite.")


def _check_univariate_distribution(distribution_type, distribution_args):
    if not isinstance(distribution_type, str):
        raise ValueError("distribution_type must be str.")

    if not isinstance(distribution_args, list):
        raise ValueError("distribution_args must be list.")

    for i in range(len(distribution_args)):
        if not isinstance(distribution_args[i], (float, int)):
            raise ValueError(f"distribution_args[{i}] is not a number.")
        if not np.isfinite(distribution_args[i]):
            raise ValueError(f"distribution_args[{i}] is an invalid parameter.")

    # To add a new distribution type: validate the
    # distribution_args here, and add unit tests.
    if distribution_type == "f":
        # Fixed-value (non-random)
        if len(distribution_args) != 1:
            raise ValueError(
                "Fixed-value mutation types (distribution_type='f') "
                "take a single selection-coefficient parameter."
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
        # dealing with lognormal distribution
        # (adding instead of multiplying the mean):
        logmean = distribution_args[0]
        logsd = distribution_args[1]
        sign = "" if distribution_type == "lp" else "-1 *"
        distribution_args = [f"return {sign}rlnorm(1, {logmean} + log(Q), {logsd});"]
        distribution_type = "s"
    elif distribution_type == "u":
        # Uniform
        if len(distribution_args) != 2 or distribution_args[0] > distribution_args[1]:
            raise ValueError(
                "Uniformly-distributed sel. coefs. (distribution_type='u') "
                "use a (min, max) parameterisation, with min <= max."
            )
        umin, umax = distribution_args
        distribution_args = [f"return runif(1, Q * {umin}, Q * {umax});"]
        distribution_type = "s"
    else:
        raise ValueError(f"{distribution_type} is not a supported distribution type.")


@attr.s(kw_only=True)
class MultivariateMutationType(object):
    """
    Class representing a "type" of mutation, allowing the mutation to affect
    fitness and/or trait(s). This design closely mirrors SLiM's MutationType
    class.

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
        (Either way, they will remain in the tree sequence).
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

    trait_distribution_type = attr.ib(default=None, type=str)
    trait_distribution_args = attr.ib(
        factory=lambda: [0], type=list, converter=_copy_converter
    )

    fitness_dominance_coeff = attr.ib(default=None, type=float)
    fitness_distribution_type = attr.ib(default="f", type=str)
    fitness_distribution_args = attr.ib(
        factory=lambda: [0], type=list, converter=_copy_converter
    )
    # TODO: is this okay if something affects traits
    convert_to_substitution = attr.ib(default=True, type=bool)
    fitness_dominance_coeff_list = attr.ib(
        default=None, type=list, converter=_copy_converter
    )
    fitness_dominance_coeff_breaks = attr.ib(
        default=None, type=list, converter=_copy_converter
    )

    def __attrs_post_init__(self):
        if (
            self.fitness_dominance_coeff is None
            and self.fitness_dominance_coeff_list is None
        ):
            self.fitness_dominance_coeff = 0.5

        if self.fitness_dominance_coeff is not None:
            if (self.fitness_dominance_coeff_list is not None) or (
                self.fitness_dominance_coeff_breaks is not None
            ):
                raise ValueError(
                    "Cannot specify both fitness_dominance_coeff "
                    "and fitness_dominance_coeff_list."
                )
            if not isinstance(self.fitness_dominance_coeff, (float, int)):
                raise ValueError("fitness_dominance_coeff must be a number.")
            if not np.isfinite(self.fitness_dominance_coeff):
                raise ValueError(
                    "Invalid fitness dominance "
                    f"coefficient {self.fitness_dominance_coeff}."
                )

        if self.fitness_dominance_coeff_list is not None:
            # disallow the inefficient and annoying length-one case
            if len(self.fitness_dominance_coeff_list) < 2:
                raise ValueError(
                    "fitness_dominance_coeff_list must have at least 2 elements."
                )
            for h in self.fitness_dominance_coeff_list:
                if not isinstance(h, (float, int)):
                    raise ValueError(
                        "fitness_dominance_coeff_list must be a list of numbers."
                    )
                if not np.isfinite(h):
                    raise ValueError(f"Invalid fitness dominance coefficient {h}.")
            if self.fitness_dominance_coeff_breaks is None:
                raise ValueError(
                    "A list of fitness dominance coefficients provided but no breaks."
                )
            if (
                len(self.fitness_dominance_coeff_list)
                != len(self.fitness_dominance_coeff_breaks) + 1
            ):
                raise ValueError(
                    "len(fitness_dominance_coeff_list) must be equal "
                    "to len(fitness_dominance_coeff_breaks) + 1"
                )
            lb = -1 * np.inf
            for b in self.fitness_dominance_coeff_breaks:
                if not isinstance(b, (float, int)):
                    raise ValueError(
                        "fitness_dominance_coeff_breaks must be a list of numbers."
                    )
                if not np.isfinite(b):
                    raise ValueError(
                        f"Invalid fitness dominance coefficient break {b}."
                    )
                if b < lb:
                    raise ValueError(
                        "fitness_dominance_coeff_breaks must be nondecreasing."
                    )
                lb = b

        _check_univariate_distribution(
            self.fitness_distribution_type, self.fitness_distribution_args
        )

        if not isinstance(self.convert_to_substitution, bool):
            raise ValueError("convert_to_substitution must be bool.")

        _check_multivariate_distribution(
            self.trait_distribution_type, self.trait_distribution_args
        )

        # The index(s) of the param in the distribution_args list that should be
        # multiplied by Q when using --slim-scaling-factor Q.
        self.fitness_Q_scaled_index = {
            "e": [0],  # mean
            "f": [0],  # fixed value
            "g": [0],  # mean
            "n": [0, 1],  # mean and sd
            "w": [0],  # scale
            "s": [],  # script types should just printout arguments
        }[self.fitness_distribution_type]

    @property
    def directly_affects_fitness(self):
        """
        Tests whether the mutation type has direct effects on fitness. This is
        defined here to be of type "f" and with fitness effect 0.0, and so
        excludes other situations that also produce only neutral mutations
        (e.g., exponential with mean 0).
        """
        return (
            self.fitness_distribution_type == "f"
            and self.fitness_distribution_args[0] == 0
        )


# superclass of DFE --> DME
class DistributionOfMutationEffects(object):
    # Remember to make sure none of the components MutationTypes are converting
    # to substitutions unless they only affect fitness
    pass
