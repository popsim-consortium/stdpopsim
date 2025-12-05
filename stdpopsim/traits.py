"""
TODO: docstring goes here (preferably, descriptive)
"""

import attr
import numpy as np


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

class Environment(Distribution):
    def __init__(self, first, *args):
        super().__init__()
        pass

class FitnessFunction(Distribution):
    def __init__(self, first, *args):
        super().__init__()
        # TODO check that distribution is either stabilizing or truncating
        pass

    # TODO figure out structure of distributions

    # For SLiM, we want to pass a list of fitness pseudo-callbacks that each have information on
    # start and end time (in generations), pop ID(s), stabilizing trait indices, stabilizing covariance,
    # truncating trait indices, truncation params.


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



# TODO: convert this to multivariate distribution class
class ProductMultivariateMutationType(object):
    pass

    @property
    def directly_affects_fitness(self):
        pass


# superclass of mutationtype
@attr.s(kw_only=True)
class Distribution(object):
    """

    Class representing a "type" of mutation.  The design closely mirrors SLiM's
    MutationType class.

    The main thing that mutation types carry is a way of drawing effect sizes
    for each new mutation. This ``distribution_type`` should be one
    of (see the SLiM manual for more information on these):

    - ``mvn``: multivariate normal, two parameters (mean, covariance matrix)

    TODO edit docstring - this needs to be a general class
    TODO generalize checks here -- they are specific to mvn
    TODO implement default for affected_trait_idx
    TODO implement truncating distribution (at minimum) + mixtures of gaussians (maybe)

    :ivar distribution_type: A str abbreviation for the parametric family of
        distributions (see below).
    :vartype distribution_type: str
    :ivar num_dimensions: Total number of dimensions considered.
    :vartype num_dimensions: int
    :ivar supported_indices: List of indices for which the distribution can be
        nonzero.
    :vartype supported_indices: list
    """

    distribution_type = attr.ib(type=str)
    distribution_args = attr.ib(
        factory=lambda: [0], type=list, converter=_copy_converter
    )



    def __attrs_post_init__(self):
        # Make sure all indices are within range, and nothing funny is going on with negative indices
        for t_idx in supported_indices:
            if t_idx < 0 or t_idx >= num_dimensions:
                raise ValueError(
                    "elements of supported_indices must be at least 0 "
                    "and less than num_dimensions."
                )
        # Make sure all indices are unique
        if len(supported_indices) != len(np.unique(supported_indices)):
            raise ValueError(
                "cannot have repeated indices in supported_indices."
            )

        if not isinstance(self.distribution_type, str):
            raise ValueError("distribution_type must be str.")

        if not isinstance(self.distribution_args, list):
            raise ValueError("distribution_args must be list.")

        if self.distribution_type == "mvn":
            # Multivariate Normal distribution with (mean, covariance) parameterization.
            if len(self.distribution_args) != 2:
                raise ValueError(
                    "multivariate normal requires two parameters "
                    "in distribution_args: "
                    "a mean vector and covariance matrix."
                )
            if not isinstance(self.distribution_args[0], np.ndarray):
                raise ValueError(
                    "mvn mean vector must be specified as numpy array."
                )
            if not isinstance(self.distribution_args[1], np.ndarray):
                raise ValueError(
                    "mvn covariance matrix must be specified as numpy array."
                )
            if len(self.distribution_args[0].shape) != 1:
                raise ValueError(
                    "mvn mean vector must be 1 dimensional."
                )
            if len(self.distribution_args[1].shape) != 2:
                raise ValueError(
                    "mvn covariance matrix must be 2 dimensional."
                )
            if self.distribution_args[0].shape[0] != self.distribution_args[1].shape[0]:
                raise ValueError(
                    "mvn mean vector and covariance matrix must be conformal."
                )
            if self.distribution_args[1].shape[0] != self.distribution_args[1].shape[1]:
                raise ValueError(
                    "mvn covariance matrix must be square."
                )
            if len(supported_indices) != self.distribution_args[0].shape[0]:
                raise ValueError(
                    "mvn dimensions must match the number of supported indices."
                )
            if not np.allclose(self.distribution_args[1], self.distribution_args[1].T):
                raise ValueError(
                    "mvn covariance matrix must be symmetric."
                )
            try:
                np.linalg.cholesky(self.distribution_args[1])
            except np.LinAlgError:
                raise ValueError(
                    "mvn covariance matrix is not positive definite."
                )
        else:
            raise ValueError(
                f"{self.distribution_type} is not a supported distribution type."
            )


# superclass of DFE
class EffectSizeDistribution(object):
    # Remember to make sure none of the components MutationTypes are converting to substitutions unless they only affect fitness
    pass
