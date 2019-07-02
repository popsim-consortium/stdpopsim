# This script is a QC implementation of the two population Tennessen Out Of Africa model
import msprime
import numpy as np
import math
import stdpopsim.models as models

class TennessenOnePopAfrica(models.Model):
    def __init__(self):
        super().__init__()
        # This model is the same as the Tennessen two population model except
        # the European population has been removed.

        # Since the Tennessen one population model largely uses parameters from
        # the Gravel et al 2001, we begin by taking the maximum likelihood value
        # from the table 2 of Gravel et al. 2011 using the Low-coverage + exons
        # data. Initially we copy over the pre- exponential growth population
        # size estimates, migration rates, and epoch times:
        generation_time = 25
        
        N_A = 7310  # Ancient population size
        N_AF0 = 14474  # Pre-modern african population size (pre and post OOA)

        T_AF = 148000 / generation_time  # Epoch transition from ancient to AF0
        T_AG = 5115 / generation_time  # start of 2nd european growth epoch

        # Next we include the additional parameters from Tennessen et al 2012
        # which include all exponential growth rates. These parameters are
        # copied from the section titled "Abundance of rare variation explained
        # by human demographic history" in Tennessen et al.

        r_AF0 = 1.66e-2  # The growth rate for the 1st african expansion

        # For the post exponenential growth popuation sizes we can calcuate the
        # population sizes at the start of the epoch using the formula f(t) =
        # x_0 * exp(r * (t_0-t))

        # African population size after 1st expansion
        N_AF1 = N_AF0 * math.exp(r_AF0 * T_AG)

        # Now we set up the population configurations. The population IDs are 0=YRI. This
        # includes both the inital sizes, growth rates, and migration rates.
        self.population_configurations = [
            msprime.PopulationConfiguration(
                initial_size=N_AF1, growth_rate=r_AF0),
        ]

        self.migration_matrix = [[0]]

        # Now we add the demographic events working backwards in time. Starting with the growth
        # slowdown in Europeans and the transition to a fixed population size in Africans.
        self.demographic_events = [
            # Reversion to fixed population size in Africans
            msprime.PopulationParametersChange(
                time=T_AG, initial_size=N_AF0, growth_rate=0, population_id=0),    
            # Change to ancestral population size pre OOA
            msprime.PopulationParametersChange(
                time=T_AF, initial_size=N_A, population_id=0)
        ]

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
            # Set the migration rate for 1st CEU growth period (for now stays same)
            msprime.MigrationRateChange(
                time=T_AG, rate=m_AF1_EU1, matrix_index=(0, 1)),
            msprime.MigrationRateChange(
                time=T_AG, rate=m_AF1_EU1, matrix_index=(1, 0)),
            # Growth slowdown in Europeans
            msprime.PopulationParametersChange(
                time=T_AG, initial_size=N_EU1, growth_rate=r_EU0, population_id=1),
            # Reversion to fixed population size in Africans
            msprime.PopulationParametersChange(
                time=T_AG, initial_size=N_AF0, growth_rate=0, population_id=0),    
            # Set the migration rate for pre CEU/CHB split
            msprime.MigrationRateChange(
                time=T_EU_AS, rate=m_AF0_B, matrix_index=(0, 1)),
            msprime.MigrationRateChange(
                time=T_EU_AS, rate=m_AF0_B, matrix_index=(1, 0)),
            # Reversion to fixed population size at the time of the CHB/CEU split
            msprime.PopulationParametersChange(
                time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),
            # Coalescence between the OOA and YRI pops
            msprime.MassMigration(
                time=T_B, source=1, destination=0, proportion=1.0),
            # Change to ancestral population size pre OOA
            msprime.PopulationParametersChange(
                time=T_AF, initial_size=N_A, population_id=0)
        ]


class BrowningAmerica(models.Model):
    def __init__(self):
        super().__init__()
        # Parameters are taken from the Methods - Simulated data section
        # Population sizes
        N_AF0 = 7310 # Initial african population size
        N_AF1 = 14474 # Second african pop. size
        N_OOA = 1861 # OOA population size
        N_CEU0 = 1032 # European population size at CEU/CHB split
        N_CHB0 = 554 # Asian population size at CEU/CHB split
        N_ADMIX0 = 30000 # Initial size of admixed population

        # Epoch times
        T_AF0_AF1 = 5920 # initial increase in african pop. size
        T_AF1_OOA = 2040 # Time of OOA event
        T_CEU_CHB = 920 # Time of european/asian split
        T_ADMIX0 = 12

        # Migration rates
        m_AF1_OOA = 1.5e-4 # Bidirectional migration rate between african and OOA pops.
        m_AF1_CEU0 = 2.5e-5 # Migration rates between AF1 and CEU0
        m_AF1_CHB0 = 7.8e-6 # Migration rates between AF1 and CHB0
        m_CEU0_CHB0 = 3.11e-5 # Migration rates between CEU0 and CHB0

        # Mass migration to create admixed populations
        mm_AF1 = 1/6
        mm_CEU0 = 2/5 # Adjusted fraction for remaining population after AF migration (5/6 * 2/5 = 1/3)
        mm_CHB0 = 1.0 # Adjusted fraction for remaining population (1/2 * 1 = 1/2)

        # Growth rates
        r_CEU0 =3.8e-3
        r_CHB0 = 4.8e-3
        r_ADMIX0 = 0.05

        # Calculate population sizes at modern (T=0) time
        N_CEU1 = N_CEU0 * math.exp(r_CEU0 * T_CEU_CHB)
        N_CHB1 = N_CHB0 * math.exp(r_CHB0 * T_CEU_CHB)
        N_ADMIX1 = N_ADMIX0 * math.exp(r_ADMIX0 * T_ADMIX0)

        # Set population sizes at T=0
        # pop0 is Africa, pop1 is Europe, pop2 is Asia, pop3 is admixed
        self.population_configurations = [
            msprime.PopulationConfiguration(
                initial_size=N_AF1, growth_rate=0),
            msprime.PopulationConfiguration(
                initial_size=N_CEU1, growth_rate=r_CEU0),
            msprime.PopulationConfiguration(
                initial_size=N_CHB1, growth_rate=r_CHB0),
            msprime.PopulationConfiguration(
                initial_size=N_ADMIX1, growth_rate=r_ADMIX0)
        ]

        # Migration matrix, all migrations to admixed population are 0
        self.migration_matrix = [
            [0, m_AF1_CEU0, m_AF1_CHB0, 0],
            [m_AF1_CEU0, 0, m_CEU0_CHB0, 0],
            [m_AF1_CHB0, m_CEU0_CHB0, 0, 0],
            [0, 0, 0, 0]
        ]

        # Now we add the demographic events working backwards in time. 
        self.demographic_events = [
            # Admixed population recoalesces with origin populations (T_ADMIX0)
            msprime.MassMigration(
                time=T_ADMIX0, source=3, destination=0, proportion=mm_AF1),
            msprime.MassMigration(
                time=T_ADMIX0 + 0.0001, source=3, destination=1, proportion=mm_CEU0),
            msprime.MassMigration(
                time=T_ADMIX0 + 0.0002, source=3, destination=2, proportion=mm_CHB0),
            # Zero out migration rate (desn't matter but added for equality to prod.)
            msprime.MigrationRateChange(
                time=T_CEU_CHB, rate=0.0),
            # CEU and CHB coalesce and set population to OOA size (T_CEU_CHB)
            msprime.MassMigration(
                time=T_CEU_CHB+0.0001, source=2, destination=1, proportion= 1.0),
            msprime.PopulationParametersChange(
                time=T_CEU_CHB+0.0002, initial_size=N_OOA, growth_rate=0.0, population_id=1),
            ## Set OOA <--> AF migration rate (T_CEU_CHB)
            msprime.MigrationRateChange(
                time=T_CEU_CHB+0.0003, rate=m_AF1_OOA, matrix_index=(0, 1)),
            msprime.MigrationRateChange(
                time=T_CEU_CHB+0.0003, rate=m_AF1_OOA, matrix_index=(1, 0)),
            # Zero out migration rate (desn't matter but added for equality to prod.)
            msprime.MigrationRateChange(
                time=T_AF1_OOA, rate=0.0),
            ## OOA and AF1 coalesce (T_OOA)
            msprime.MassMigration(
                time=T_AF1_OOA+0.0001, source=1, destination=0, proportion= 1.0),
            ## AF1 -> AF0 population size change (T_AF0_AF1)
            msprime.PopulationParametersChange(
                time=T_AF0_AF1, initial_size=N_AF0, population_id=0),
        ]


class RagsdaleArchaic(models.Model):
    def __init__(self):
        super().__init__()

        # All parameters were taken from table 1 of Ragsdale et al. (2019)
        generation_time = 29
        
        # Population sizes
        N_0 = 3600 # Size of archaic populations
        N_YRI = 13900 # Fixed size of YRI population
        N_B = 880 # Size of OOA population
        N_CEU0 = 2300 # Size of CEU population at CEU-CHB split
        N_CHB0 = 650 # Size of CHB population at CEU-CHB split

        # Population growth parameters
        r_CEU = 0.125e-2
        r_CHB = 0.372e-2

        # Migration parameters
        m_AF_B = 52.2e-5
        m_YRI_CEU = 2.48e-5
        m_YRI_CHB = 0
        m_CEU_CHB = 11.3e-5
        m_AF_ARCHAF = 1.98e-5
        m_OOA_NEAN =  0.825e-5

        # Epoch times
        T_AF = 300e3/generation_time
        T_OOA = 60.7e3/generation_time
        T_CEU_CHB  = 36e3/generation_time
        T_ARCHAF_split = 499e3/generation_time
        T_ARCHAF_mig = 125e3/generation_time
        T_NEAN_split = 559e3/generation_time
        T_ARCH_ADMIX_end = 18.7e3/generation_time

        # Calculate population sizes at modern (T=0) time
        N_CEU1 = N_CEU0 * math.exp(r_CEU * T_CEU_CHB)
        N_CHB1 = N_CHB0 * math.exp(r_CHB * T_CEU_CHB)

        # Set population sizes at T=0
        # pop0 is Africa, pop1 is Europe, pop2 is Asia, pop3 is Neanderthal, pop4 is
        # archaic african
        self.population_configurations = [
            msprime.PopulationConfiguration(
                initial_size=N_YRI, growth_rate=0),
            msprime.PopulationConfiguration(
                initial_size=N_CEU1, growth_rate=r_CEU),
            msprime.PopulationConfiguration(
                initial_size=N_CHB1, growth_rate=r_CHB),
            msprime.PopulationConfiguration(
                initial_size=N_0, growth_rate=0),
            msprime.PopulationConfiguration(
                initial_size=N_0, growth_rate=0)
        ]

        # Setup initial migration matrix
        self.migration_matrix = [
            [0,         m_YRI_CEU,  m_YRI_CHB, 0, 0], # noqa
            [m_YRI_CEU, 0,          m_CEU_CHB, 0, 0], # noqa
            [m_YRI_CHB, m_CEU_CHB,  0,         0, 0], # noqa
            [0,         0,          0,         0, 0], # noqa
            [0,         0,          0,         0, 0] # noqa
        ]

        self.demographic_events = [
            # Migration between YRI and ARCHAF(E1)
            msprime.MigrationRateChange(
                time=T_ARCH_ADMIX_end, rate=m_AF_ARCHAF, matrix_index=(0, 4)),
            msprime.MigrationRateChange(
                time=T_ARCH_ADMIX_end, rate=m_AF_ARCHAF, matrix_index=(4, 0)),
            # Migration between CEU and NEAN(E1)
            msprime.MigrationRateChange(
                time=T_ARCH_ADMIX_end, rate=m_OOA_NEAN, matrix_index=(1, 3)),
            msprime.MigrationRateChange(
                time=T_ARCH_ADMIX_end, rate=m_OOA_NEAN, matrix_index=(3, 1)),
            # Migration between CHB and NEAN(E1)
            msprime.MigrationRateChange(
                time=T_ARCH_ADMIX_end, rate=m_OOA_NEAN, matrix_index=(2, 3)),
            msprime.MigrationRateChange(
                time=T_ARCH_ADMIX_end, rate=m_OOA_NEAN, matrix_index=(3, 2)),
            # Coalescence of CHB into CEU (E2)
            msprime.MassMigration(
                time=T_CEU_CHB, source=2, dest=1, proportion=1.0),
            # Reset migration rates (E2)(redundant)*
            msprime.MigrationRateChange(
                time=T_CEU_CHB, rate=0.0),
            # Migration rate change between OOA(CEU) and AF(YRI)(E2)
            msprime.MigrationRateChange(
                time=T_CEU_CHB, rate=m_AF_B, matrix_index=(0, 1)),
            msprime.MigrationRateChange(
                time=T_CEU_CHB, rate=m_AF_B, matrix_index=(1, 0)),    
            # Migration between YRI and ARCHAF (E2)(redundant without mig. rate reset)*
            msprime.MigrationRateChange(
                time=T_CEU_CHB, rate=m_AF_ARCHAF, matrix_index=(0, 4)),
            msprime.MigrationRateChange(
                time=T_CEU_CHB, rate=m_AF_ARCHAF, matrix_index=(4, 0)),
            # Migration between CEU and NEAN (E2)(redundant without mig. rate reset)*
            msprime.MigrationRateChange(
                time=T_CEU_CHB, rate=m_OOA_NEAN, matrix_index=(1, 3)),
            msprime.MigrationRateChange(
                time=T_CEU_CHB, rate=m_OOA_NEAN, matrix_index=(3, 1)),
            # CEU change to fixed population size at the time of the CHB/CEU coal. (E2)
            msprime.PopulationParametersChange(
                time=T_CEU_CHB, initial_size=N_B, growth_rate=0, population_id=1),
            # Coalescence between the OOA and AF pops (E3)
            msprime.MassMigration(
                time=T_OOA, source=1, destination=0, proportion=1.0),
            # Reset migration rates (E3)
            msprime.MigrationRateChange(
                time=T_OOA, rate=0.0),
            # Migration between YRI and ARCHAF (E3)
            msprime.MigrationRateChange(
                time=T_OOA, rate=m_AF_ARCHAF, matrix_index=(0, 4)),
            msprime.MigrationRateChange(
                time=T_OOA, rate=m_AF_ARCHAF, matrix_index=(4, 0)),
            # Migration between archaic african and african pop. "ends" (E4)
            msprime.MigrationRateChange(
                time=T_ARCHAF_mig, rate=0),
            # AF reverts to ancestral population size pre OOA (E5)
            msprime.PopulationParametersChange(
                time=T_AF, initial_size=N_0, population_id=0),
            # Archaic AF population coalesces into AF (E6)
            msprime.MassMigration(
                time=T_ARCHAF_split, source = 4, dest = 0, proportion=1.0),
            # NEAN pop. coalesces into AF (E7)
            msprime.MassMigration(
                time=T_NEAN_split, source = 3, dest = 0, proportion=1.0)
        ]