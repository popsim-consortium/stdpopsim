"""
Genome and demographic model definitions for Escherichia coli.
"""
import msprime

import stdpopsim.models as models
import stdpopsim.genomes as genomes


###########################################################
#
# Genome definition
#
###########################################################

_chromosomes = []
_chromosomes.append(genomes.Chromosome(
        name="chr",
        length=4641652,
        default_mutation_rate=1e-5+2e-4,
        default_recombination_rate=0.0))
# mean_conversion_rate=8.9e-11 # not implemented yet!
# mean_conversion_length=542 # not implemented yet!

#: :class:`stdpopsim.Genome` definition for E. Coli.
# Chromosome length data is based on strain K-12.
genome = genomes.Genome(
    species="e_coli",
    chromosomes=_chromosomes,
    default_genetic_map=None)


###########################################################
#
# Demographic models
#
###########################################################


class LapierreConstant(models.Model):
    """
    The constant population size model from Lapierre et al. 2016
    https://doi.org/10.1093/molbev/msw048
    """

    def __init__(self):

        N_e = 1.8e8
        # Single population
        self.population_configurations = [
            msprime.PopulationConfiguration(initial_size=N_e),
        ]
        self.migration_matrix = [[0]]
        self.demographic_events = []
