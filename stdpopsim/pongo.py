"""
Orangutans.
"""
import math

import msprime

import stdpopsim.models as models
# import stdpopsim.chromosomes as chromosomes


# # Maybe good to package this part differently if wanting to use
# # some information but not all.
# # E.g. If you want mean mutation rate and recomb rates, but user-def length
# chrom = chromosomes.Chromosome(
#     name="chr22",
#     length=50818468,  # from H.sapiens chr22, Wikipedia
#     mean_mutation_rate=2e-8,  # assumption from Locke etal
#     mean_recombination_rate=.95*1e-8)  # from Locke etal


# ETC


class LockeEtAlPongoIM(models.Model):
    '''
    http://doi.org/10.1038/nature09687
    '''

    def __init__(self):

        # Parameters from paper:
        # ancestral size, before split
        Na = 17934

        # time of split
        T_split_years = 403149
        # get split time in units of generations
        generation_time = 20
        T_split = T_split_years / generation_time

        # proportion of ancestral pop to branch B
        s = 0.592

        # sizes at split
        Na_B = 17934*s
        Na_S = 17934*(1-s)

        # present sizes
        N_B = 8805
        N_S = 37661

        # get growth rates
        r_B = -1*math.log(Na_B/N_B)/T_split
        r_S = -1*math.log(Na_S/N_S)/T_split

        # migration rates (UNITS?)
        m_S_B = 0.395
        m_B_S = 0.239

        # Population IDs correspond to their indexes in the population
        # configuration array. Therefore, we have 0=B and 1=S
        # initially.
        self.population_configurations = [
            msprime.PopulationConfiguration(initial_size=N_B, growth_rate=r_B),
            msprime.PopulationConfiguration(initial_size=N_S, growth_rate=r_S),
        ]
        self.migration_matrix = [
            [      0, m_B_S],  # NOQA
            [m_S_B,       0],  # NOQA
        ]
        self.demographic_events = [
            # CEU and CHB merge into B with rate changes at T_EU_AS
            msprime.MassMigration(
                time=T_split, source=1, destination=0, proportion=1.0),
            msprime.MigrationRateChange(time=T_split, rate=0),
            msprime.PopulationParametersChange(
                time=T_split, initial_size=Na, growth_rate=0, population_id=0),
        ]
