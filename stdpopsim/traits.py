"""
TODO: docstring goes here (preferably, descriptive)
"""

import attr
import numpy as np


class Traits:
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


# TODO: convert this to multivariate distribution class
class ProductMultivariateMutationType(object):
    pass

    @property
    def directly_affects_fitness(self):
        pass


# superclass of mutationtype
@attr.s(kw_only=True)
class Distribution:
    """

    Class representing a "type" of mutation.  The design closely mirrors SLiM's
    MutationType class.

    The main thing that mutation types carry is a way of drawing effect sizes
    for each new mutation. This ``distribution_type`` should be one
    of (see the SLiM manual for more information on these):

    - ``mvn``: multivariate normal, two parameters (mean, covariance matrix)

    TODO edit docstring - this needs to be a general class
    TODO change affected_trait_indices to something more general
    TODO generalize checks here -- they are specific to mvn
    TODO implement default for affected_trait_idx
    TODO implement truncating distribution (at minimum) + mixtures of gaussians (maybe)

    :ivar distribution_type: A str abbreviation for the distribution of
        fitness effects that each new mutation of this type draws from (see below).
    :vartype distribution_type: str
    :ivar num_traits: Total number of traits being considered.
    :vartype num_traits: int
    :ivar affected_trait_indices: List of indices of the traits affected by mutations of this type
    :vartype affected_trait_indices: list
    """

    distribution_type = attr.ib(type=str)
    distribution_args = attr.ib(
        factory=lambda: [0], type=list, converter=_copy_converter
    )



    def __attrs_post_init__(self):
        # Make sure this affects at least one trait
        if len(affected_trait_indices) == 0:
            raise ValueError("TODO")
        # Make sure all indices are within range, and nothing funny is going on with negative indices
        for t_idx in affected_trait_indices:
            if t_idx < 0 or t_idx >= num_traits:
                raise ValueError("TODO")
        # Make sure fitness is not affected
        if 0 in affected_trait_indices:
            raise ValueError("TODO")
        # Make sure all indices are unique
        if len(affected_trait_indices) != len(np.unique(affected_trait_indices)):
            raise ValueError("TODO")

        if not isinstance(self.distribution_type, str):
            raise ValueError("distribution_type must be str.")

        if not isinstance(self.distribution_args, list):
            raise ValueError("distribution_args must be list.")

        if self.distribution_type == "mvn":
            # Multivariate Normal distribution with (mean, covariance) parameterization.
            if len(self.distribution_args) != 2:
                raise ValueError("TODO")
            if not isinstance(self.distribution_args[0], np.ndarray):
                raise ValueError("TODO")
            if not isinstance(self.distribution_args[1], np.ndarray):
                raise ValueError("TODO")
            if len(self.distribution_args[0].shape) != 1:
                raise ValueError("TODO")
            if len(self.distribution_args[1].shape) != 2:
                raise ValueError("TODO")
            if self.distribution_args[0].shape[0] != self.distribution_args[1].shape[0]:
                raise ValueError("TODO")
            if self.distribution_args[1].shape[0] != self.distribution_args[1].shape[1]:
                raise ValueError("TODO")
            if len(affected_trait_indices) != self.distribution_args[0].shape[0]:
                raise ValueError("TODO")
            if not np.allclose(self.distribution_args[1], self.distribution_args[1].T):
                raise ValueError("TODO")
            try:
                np.linalg.cholesky(self.distribution_args[1])
            except np.LinAlgError:
                raise ValueError("Covariance matrix is not positive definite.")
        else:
            raise ValueError(
                f"{self.distribution_type} is not a supported distribution type."
            )


# superclass of DFE
class EffectSizeDistribution(object):
    # Remember to make sure none of the components MutationTypes are converting to substitutions unless they only affect fitness
    pass
