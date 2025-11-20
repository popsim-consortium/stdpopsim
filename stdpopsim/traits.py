"""
TODO: docstring goes here (preferably, descriptive)
"""

import attr
import numpy as np


class Traits:
    def __init__(self, model, genotype_to_phenotype_map, fitness_function_list, num_traits):
        self.demographic_model = model
        self.num_traits = num_traits
        
        assert self.is_valid_g2p(genotype_to_phenotype_map)
        assert self.is_valid_fitness_function(fitness_function_list)

        # make these attributes
        self._genotype_to_phenotype_map = genotype_to_phenotype_map
        self._fitness_function = fitness_function_list

    def is_valid_g2p(self, genotype_to_phenotype_map):
        # check that populations/time points are consistent with demographic model
        return

    def is_valid_fitness_function(self, fitness_function):
        # for each fitness function in list: check that populations/time points 
        # are consistent with demographic model. also check that the list is 
        # comprehensive.
        return

class GenotypeToPhenotypeMap:
    def __init__(self, envs, transformation):
        return

class FitnessFunction:
    def __init__(self, start, end, populations, distributions):
        pass

    # TODO figure out structure of distributions

    # For SLiM, we want to pass a list of fitness pseudo-callbacks that each have information on
    # start and end time (in generations), pop ID(s), stabilizing trait indices, stabilizing covariance,
    # truncating trait indices, truncation params.


# TODO: is this too OO?
class ProductMultivariateMutationType(object):
    pass

    @property
    def directly_affects_fitness(self):
        pass


# superclass of mutationtype
@attr.s(kw_only=True)
class MultivariateMutationType(object):
    """

    Class representing a "type" of mutation.  The design closely mirrors SLiM's
    MutationType class.

    The main thing that mutation types carry is a way of drawing effect sizes
    for each new mutation. This ``distribution_type`` should be one
    of (see the SLiM manual for more information on these):

    - ``mvn``: multivariate normal, two parameters (mean, covariance matrix)

    TODO

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
        if len(affected_trait_indices) != len(np.unique(affected_trait_indices):
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
            except LinAlgError:
                raise ValueError("Covariance matrix is not positive definite.")
         else:
             raise ValueError(
                f"{self.distribution_type} is not a supported distribution type."
             )

    
# superclass of DFE
class EffectSizeDistribution(object):
    # Remember to make sure none of the components MutationTypes are converting to substitutions unless they only affect fitness
    pass
