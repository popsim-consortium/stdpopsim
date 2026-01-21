"""
TODO: docstring goes here (preferably, descriptive)
"""

import attr
import numpy as np


#TODO: somewhere we need to check that the number of traits
# is consistent between the MultivariateMutationTypes, EffectSizeDistributions,
# Environments, and FitnessFunctions


# TODO: this was moved from dfe.py so technically
# this would break backward compatibility, but
# I believe this method is no longer used
# anywhere in dfe.py, and it's private so this should be chill
def _copy_converter(x):
    if isinstance(x, list):
        x = x.copy()
    return x


class Traits(object):
    def __init__(self, model, environments, genetic_val_transform, fitness_functions, num_traits):
        """
        environments: list of tuples (start_time, [pop_indices], Environment object)
        fitness_functions: list of tuples (start_time, [pop_indices], FitnessFunction object)
        """
        self.demographic_model = model
        self.num_traits = num_traits

        assert self.check_model_params(environments)
        assert self.is_valid_g2p(genetic_val_transform)
        assert self.check_model_params(fitness_functions)

        # make these attributes
        self._genetic_val_transform = genetic_val_transform
        self._fitness_function = fitness_functions

    def is_valid_g2p(self, genetic_val_transform):
        # check that provided transform is valid
        return

    def check_model_params(self, params):
        # sort elements by epoch start time (which is the backwards-in-time generation time)
        # tuples should look like (epoch_start, [list of applicable pops], distribution params...)
        # check populations and generation times make sense?
        # distributions should propagate backwards in time
        # so to apply the same distribution everywhere, use the tuple (0, "all", ...)
        # also check that num_traits matches distribution params
        pass

    # TODO check g2p against fitness_functions -- traits with binary fitness fn should have binary phenos

    def add_fitness_function(self, traits, distribution_type, distribution_args, spatiotemporal=None):
        new_fitness = FitnessFunction(traits, distribution_type, distribution_args, spacetime=spatiotemporal)
        self._fitness_functions.extend(new_fitness)


class Environment(Distribution):
    def __init__(self, first, *args):
        super().__init__()
        pass

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

    supported_function_types = [] # TODO

    traits = attr.ib(type=list)
    function_type = attr.ib(type=str)
    function_args = attr.ib(type=tuple)
    # spacetime = attr.ib(type=list)

    def __attrs_post_init__(self):
        if len(self.traits) < 1:
            raise ValueError(
                "At least one trait must be specified."
            )
        if self.function_type not in self.supported_function_types:
            raise ValueError(
                "Proposed fitness function not supported at this time."
            )


@attr.s(kw_only=True)
class MultivariateEffectSizeDistribution(Distribution):
    def __attrs_post_init__(self):
        # Make sure fitness is not affected
        if 0 in supported_indices:
            raise ValueError(
                "distributions cannot simultaneously directly affect fitness "
                "and traits. that is, if 0 (fitness) is a supported index "
                "no other index can also be supported in supported_indices."
            )


def _check_multivariate_distribution(distribution_type, distribution_args):
    if not isinstance(self.distribution_type, str):
        raise ValueError("distribution_type must be str.")

    if not isinstance(self.distribution_args, list):
        raise ValueError("distribution_args must be list.")

    if distribution_type == "mvn":
        # Multivariate Normal distribution with (mean, covariance, indices) parameterization.
        if len(distribution_args) != 3:
            raise ValueError(
                "multivariate normal requires three parameters "
                "in distribution_args: "
                "a mean vector, covariance matrix, and list of dimensions "
                "that are not zero."
            )
        if not isinstance(distribution_args[0], np.ndarray):
            raise ValueError(
                "mvn mean vector must be specified as numpy array."
            )
        if not isinstance(distribution_args[1], np.ndarray):
            raise ValueError(
                "mvn covariance matrix must be specified as numpy array."
            )
        if len(distribution_args[0].shape) != 1:
            raise ValueError(
                "mvn mean vector must be 1 dimensional."
            )
        if len(distribution_args[1].shape) != 2:
            raise ValueError(
                "mvn covariance matrix must be 2 dimensional."
            )
        if distribution_args[0].shape[0] != distribution_args[1].shape[0]:
            raise ValueError(
                "mvn mean vector and covariance matrix must be conformal."
            )
        if distribution_args[1].shape[0] != distribution_args[1].shape[1]:
            raise ValueError(
                "mvn covariance matrix must be square."
            )
        if len(supported_indices) != distribution_args[0].shape[0]:
            raise ValueError(
                "mvn dimensions must match the number of supported indices."
            )
        if not np.allclose(distribution_args[1], distribution_args[1].T):
            raise ValueError(
                "mvn covariance matrix must be symmetric."
            )
        try:
            np.linalg.cholesky(distribution_args[1])
        except np.LinAlgError:
            raise ValueError(
                "mvn covariance matrix is not positive definite."
            )

        supported_indices = distribution_args[2]
        # TODO: check dtype of supported_indices
        for t_idx in supported_indices:
            if t_idx < 0:
                raise ValueError(
                    "elements of supported_indices must be at least 0 "
                    "and less than num_dimensions."
                )
        # Make sure all indices are unique
        if len(supported_indices) != len(np.unique(supported_indices)):
            raise ValueError(
                "cannot have repeated indices in supported_indices."
            )

    else:
        raise ValueError(
            f"{distribution_type} is not a supported distribution type."
        )


    pass


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
        distribution_args = [
            f"return {sign}rlnorm(1, {logmean} + log(Q), {logsd});"
        ]
        distribution_type = "s"
    elif distribution_type == "u":
        # Uniform
        if (
            len(distribution_args) != 2
            or distribution_args[0] > distribution_args[1]
        ):
            raise ValueError(
                "Uniformly-distributed sel. coefs. (distribution_type='u') "
                "use a (min, max) parameterisation, with min <= max."
            )
        umin, umax = distribution_args
        distribution_args = [f"return runif(1, Q * {umin}, Q * {umax});"]
        distribution_type = "s"
    else:
        raise ValueError(
            f"{distribution_type} is not a supported distribution type."
        )


@attr.s(kw_only=True)
class MultivariateMutationType(object):
    """
    Class representing a "type" of mutation, allowing the mutation to affect
    fitness and/or trait(s). This design closely mirrors SLiM's MutationType
    class.

    TODO: write a good docstring
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
                    "Cannot specify both fitness_dominance_coeff and fitness_dominance_coeff_list."
                )
            if not isinstance(self.fitness_dominance_coeff, (float, int)):
                raise ValueError("fitness_dominance_coeff must be a number.")
            if not np.isfinite(self.fitness_dominance_coeff):
                raise ValueError(
                    f"Invalid fitness dominance coefficient {self.fitness_dominance_coeff}."
                )

        if self.fitness_dominance_coeff_list is not None:
            # disallow the inefficient and annoying length-one case
            if len(self.fitness_dominance_coeff_list) < 2:
                raise ValueError("fitness_dominance_coeff_list must have at least 2 elements.")
            for h in self.fitness_dominance_coeff_list:
                if not isinstance(h, (float, int)):
                    raise ValueError("fitness_dominance_coeff_list must be a list of numbers.")
                if not np.isfinite(h):
                    raise ValueError(f"Invalid fitness dominance coefficient {h}.")
            if self.fitness_dominance_coeff_breaks is None:
                raise ValueError(
                    "A list of fitness dominance coefficients provided but no breaks."
                )
            if len(self.fitness_dominance_coeff_list) != len(self.fitness_dominance_coeff_breaks) + 1:
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
                    raise ValueError(f"Invalid fitness dominance coefficient break {b}.")
                if b < lb:
                    raise ValueError("fitness_dominance_coeff_breaks must be nondecreasing.")
                lb = b

        _check_univariate_distribution(
            self.fitness_distribution_type,
            self.fitness_distribution_args
        )

        if not isinstance(self.convert_to_substitution, bool):
            raise ValueError("convert_to_substitution must be bool.")

        _check_multivariate_distribution(
            self.trait_distribution_type,
            self.trait_distribution_args
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
        Tests whether the mutation type has direct effects on fitness. This is defined here to
        be of type "f" and with fitness effect 0.0, and so excludes other situations
        that also produce only neutral mutations (e.g., exponential with mean 0).
        """
        return self.fitness_distribution_type == "f" and self.fitness_distribution_args[0] == 0

# superclass of DFE
class EffectSizeDistribution(object):
    # Remember to make sure none of the components MutationTypes are converting to substitutions unless they only affect fitness
    pass
