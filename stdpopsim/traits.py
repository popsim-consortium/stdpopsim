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