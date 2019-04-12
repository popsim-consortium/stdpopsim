# This script is a QC implementation of the two population Tennessen Out Of Africa model
import msprime
import numpy as np
import math
import stdpopsim.models as models

class TennessenTwoPopOutOfAfrica(models.Model):
    def __init__(self):
        super().__init__()
        # Since the Tennessen two population model largely uses parameters from the Gravel et al 2001, 
        # we begin by taking the maximum likelihood value from the table 2 of Gravel et al. 2011 using 
        # the Low-coverage + exons data. We ignore all values related to the asian (AS) population as 
        # it is not present in the Tennessen two population model. Initially we copy over the pre-
        # exponential growth population size estimates, migration rates, and epoch times:
        generation_time = 25
        
        N_A = 7310  # Ancient population size
        N_AF0 = 14474  # Pre-modern african population size (pre and post OOA)
        N_B = 1861  # OOA population size, pre-expansion
        N_EU0 = 1032  # European population size, pre-expansion

        m_AF0_B = 15e-5  # migration rate from pre-expansion africa to pre-expansion OOA
        m_AF0_EU0 = 2.5e-5  # migration rate from pre-expansion africa to 1st-expansion european
        m_AF1_EU1 = 2.5e-5  # migration rate from pre-expansion africa to 2nd-expansion european

        T_AF = 148000 / generation_time  # Epoch transition from ancient to AF0
        T_B = 51000 / generation_time  # OOA time
        # The european asian split time, begins 1st growth period
        T_EU_AS = 23000 / generation_time

        # Next we include the additional parameters from Tennessen et al 2012 which include all 
        # exponential growth rates and the time of the second round of growth in the European 
        # population/first round in the African population. These parameters are copied from the 
        # section titled "Abundance of rare variation explained by human demographic history" in 
        # Tennessen et al.
        r_EU0 = 0.307e-2  # The growth rate for the 1st european expansion
        r_EU1 = 1.95e-2  # The growth rate for the 2nd european expansion
        r_AF0 = 1.66e-2  # The growth rate for the 1st african expansion

        T_AG = 5115 / generation_time  # start of 2nd european growth epoch

        # For the post exponenential growth popuation sizes we can calcuate the population 
        # sizes at the start of the epoch using the formula f(t) = x_0 * exp(r * (t_0-t))
        # European population size after 1st expansion
        N_EU1 = N_EU0 * math.exp(r_EU0 * (T_EU_AS-T_AG))
        # European population size after 2nd expansion
        N_EU2 = N_EU1 * math.exp(r_EU1 * T_AG)
        # African population size after 1st expansion
        N_AF1 = N_AF0 * math.exp(r_AF0 * T_AG)

        # Now we set up the population configurations. The population IDs are 0=CEU and 1=YRI. This
        # includes both the inital sizes, growth rates, and migration rates.
        self.population_configurations = [
            msprime.PopulationConfiguration(
                initial_size=N_AF1, growth_rate=r_AF0),
            msprime.PopulationConfiguration(
                initial_size=N_EU2, growth_rate=r_EU1)
        ]
        self.migration_matrix = [
            [0, m_AF1_EU1],
            [m_AF1_EU1,       0],
        ]

        # Now we add the demographic events working backwards in time. Starting with the growth
        # slowdown in Europeans and the transition to a fixed population size in Africans.
        self.demographic_events = [
            # Reversion to fixed population size in Africans
            msprime.PopulationParametersChange(
                time=T_AG, initial_size=N_AF0, growth_rate=0, population_id=0),
            # Growth slowdown in Europeans
            msprime.PopulationParametersChange(
                time=T_AG, initial_size=N_EU1, growth_rate=r_EU0, population_id=1),
            # Set the migration rate for 1st CEU growth period (for now stays same)
            msprime.MigrationRateChange(
                time=T_AG, rate=m_AF1_EU1, matrix_index=(0, 1)),
            msprime.MigrationRateChange(
                time=T_AG, rate=m_AF1_EU1, matrix_index=(1, 0)),
            # Reversion to fixed population size at the time of the CHB/CEU split
            msprime.PopulationParametersChange(
                time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),
            # Set the migration rate for pre CEU/CHB split
            msprime.MigrationRateChange(
                time=T_EU_AS, rate=m_AF0_B, matrix_index=(0, 1)),
            msprime.MigrationRateChange(
                time=T_EU_AS, rate=m_AF0_B, matrix_index=(1, 0)),
            # Coalescence between the OOA and YRI pops
            msprime.MassMigration(
                time=T_B, source=1, destination=0, proportion=1.0),
            # Change to ancestral population size pre OOA
            msprime.PopulationParametersChange(
                time=T_AF, initial_size=N_A, population_id=0)
        ]
