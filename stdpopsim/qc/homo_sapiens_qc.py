import math

import msprime

import stdpopsim
import stdpopsim.models as models


_species = stdpopsim.get_species("HomSap")

# Some generic populations to use for qc
population_sample_0 = models.Population("sampling_0",
                                        "Population that samples at time 0",
                                        0)
population_sample_none = models.Population("sampling_none",
                                           "Population that does not sample",
                                           None)


class TennessenOnePopAfrica(models.DemographicModel):
    populations = [population_sample_0]

    def __init__(self):
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

        # Now we add the demographic events working backwards in time. Starting with the
        # growth slowdown in Europeans and the transition to a fixed population size in
        # Africans.
        self.demographic_events = [
            # Reversion to fixed population size in Africans
            msprime.PopulationParametersChange(
                time=T_AG, initial_size=N_AF0, growth_rate=0, population_id=0),
            # Change to ancestral population size pre OOA
            msprime.PopulationParametersChange(
                time=T_AF, initial_size=N_A, population_id=0)
        ]


_species.get_demographic_model("Africa_1T12").register_qc(TennessenOnePopAfrica())


class TennessenTwoPopOutOfAfrica(models.DemographicModel):
    populations = [population_sample_0] * 2

    def __init__(self):
        # Since the Tennessen two population model largely uses parameters from
        # the Gravel et al 2001, we begin by taking the maximum likelihood
        # value from the table 2 of Gravel et al. 2011 using the Low-coverage +
        # exons data. We ignore all values related to the asian (AS) population
        # as it is not present in the Tennessen two population model. Initially
        # we copy over the pre- exponential growth population size estimates,
        # migration rates, and epoch times:
        generation_time = 25

        N_A = 7310  # Ancient population size
        N_AF0 = 14474  # Pre-modern african population size (pre and post OOA)
        N_B = 1861  # OOA population size, pre-expansion
        N_EU0 = 1032  # European population size, pre-expansion

        m_AF0_B = 15e-5  # migration from pre-expansion africa to pre-expansion OOA
        m_AF1_EU1 = 2.5e-5  # migration from pre-expansion africa to 2nd-expansion euro

        T_AF = 148000 / generation_time  # Epoch transition from ancient to AF0
        T_B = 51000 / generation_time  # OOA time
        # The european asian split time, begins 1st growth period
        T_EU_AS = 23000 / generation_time

        # Next we include the additional parameters from Tennessen et al 2012
        # which include all exponential growth rates and the time of the second
        # round of growth in the European population/first round in the African
        # population. These parameters are copied from the section titled
        # "Abundance of rare variation explained by human demographic history"
        # in Tennessen et al.
        r_EU0 = 0.307e-2  # The growth rate for the 1st european expansion
        r_EU1 = 1.95e-2  # The growth rate for the 2nd european expansion
        r_AF0 = 1.66e-2  # The growth rate for the 1st african expansion

        T_AG = 5115 / generation_time  # start of 2nd european growth epoch

        # For the post exponenential growth popuation sizes we can calcuate the
        # population sizes at the start of the epoch using the formula
        # f(t) = x_0 * exp(r * (t_0-t)) European population size after 1st expansion
        N_EU1 = N_EU0 * math.exp(r_EU0 * (T_EU_AS-T_AG))
        # European population size after 2nd expansion
        N_EU2 = N_EU1 * math.exp(r_EU1 * T_AG)
        # African population size after 1st expansion
        N_AF1 = N_AF0 * math.exp(r_AF0 * T_AG)

        # Now we set up the population configurations. The population IDs are
        # 0=CEU and 1=YRI. This includes both the inital sizes, growth rates,
        # and migration rates.
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

        # Now we add the demographic events working backwards in time. Starting
        # with the growth slowdown in Europeans and the transition to a fixed
        # population size in Africans.
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


_species.get_demographic_model(
        "OutOfAfrica_2T12").register_qc(TennessenTwoPopOutOfAfrica())


class BrowningAmerica(models.DemographicModel):
    populations = [population_sample_0] * 4

    def __init__(self):
        # Parameters are taken from the Methods - Simulated data section
        # Population sizes
        N_AF0 = 7310  # Initial african population size
        N_AF1 = 14474  # Second african pop. size
        N_OOA = 1861  # OOA population size
        N_CEU0 = 1032  # European population size at CEU/CHB split
        N_CHB0 = 554  # Asian population size at CEU/CHB split
        N_ADMIX0 = 30000  # Initial size of admixed population

        # Epoch times
        T_AF0_AF1 = 5920  # initial increase in african pop. size
        T_AF1_OOA = 2040  # Time of OOA event
        T_CEU_CHB = 920  # Time of european/asian split
        T_ADMIX0 = 12

        # Migration rates
        m_AF1_OOA = 1.5e-4  # Bidirectional migration rate between african and OOA pops.
        m_AF1_CEU0 = 2.5e-5  # Migration rates between AF1 and CEU0
        m_AF1_CHB0 = 7.8e-6  # Migration rates between AF1 and CHB0
        m_CEU0_CHB0 = 3.11e-5  # Migration rates between CEU0 and CHB0

        # Mass migration to create admixed populations
        mm_AF1 = 1/6
        # Adjusted fraction for remaining population after AF migration (5/6 * 2/5 = 1/3)
        mm_CEU0 = 2/5
        # Adjusted fraction for remaining population (1/2 * 1 = 1/2)
        mm_CHB0 = 1.0

        # Growth rates
        r_CEU0 = 3.8e-3
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
                time=T_CEU_CHB+0.0001, source=2, destination=1, proportion=1.0),
            msprime.PopulationParametersChange(
                time=T_CEU_CHB+0.0002, initial_size=N_OOA, growth_rate=0.0,
                population_id=1),
            # Set OOA <--> AF migration rate (T_CEU_CHB)
            msprime.MigrationRateChange(
                time=T_CEU_CHB+0.0003, rate=m_AF1_OOA, matrix_index=(0, 1)),
            msprime.MigrationRateChange(
                time=T_CEU_CHB+0.0003, rate=m_AF1_OOA, matrix_index=(1, 0)),
            # Zero out migration rate (desn't matter but added for equality to prod.)
            msprime.MigrationRateChange(
                time=T_AF1_OOA, rate=0.0),
            # OOA and AF1 coalesce (T_OOA)
            msprime.MassMigration(
                time=T_AF1_OOA+0.0001, source=1, destination=0, proportion=1.0),
            # AF1 -> AF0 population size change (T_AF0_AF1)
            msprime.PopulationParametersChange(
                time=T_AF0_AF1, initial_size=N_AF0, population_id=0),
        ]


_species.get_demographic_model(
        "AmericanAdmixture_4B11").register_qc(BrowningAmerica())


class RagsdaleArchaic(models.DemographicModel):
    populations = [population_sample_0] * 3 + [population_sample_none] * 2

    def __init__(self):

        # All parameters were taken from table 1 of Ragsdale et al. (2019)
        generation_time = 29

        # Population sizes
        N_0 = 3600  # Size of archaic populations
        N_YRI = 13900  # Fixed size of YRI population
        N_B = 880  # Size of OOA population
        N_CEU0 = 2300  # Size of CEU population at CEU-CHB split
        N_CHB0 = 650  # Size of CHB population at CEU-CHB split

        # Population growth parameters
        r_CEU = 0.125e-2
        r_CHB = 0.372e-2

        # Migration parameters
        m_AF_B = 52.2e-5
        m_YRI_CEU = 2.48e-5
        m_YRI_CHB = 0
        m_CEU_CHB = 11.3e-5
        m_AF_ARCHAF = 1.98e-5
        m_OOA_NEAN = 0.825e-5

        # Epoch times
        T_AF = 300e3/generation_time
        T_OOA = 60.7e3/generation_time
        T_CEU_CHB = 36e3/generation_time
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
                time=T_ARCHAF_split, source=4, dest=0, proportion=1.0),
            # NEAN pop. coalesces into AF (E7)
            msprime.MassMigration(
                time=T_NEAN_split, source=3, dest=0, proportion=1.0)
        ]


_species.get_demographic_model(
        "OutOfAfricaArchaicAdmixture_5R19").register_qc(RagsdaleArchaic())


class KammAncientSamples(models.DemographicModel):
    """
    Demographic inferred by momi described in Kamm et al. (2019). The model is
    illustrated in Figure 3, with parameters given in Table 2.
    """
    populations = [
        population_sample_0,
        models.Population("", "", 8e3 / 25),
        population_sample_0,
        models.Population("", "", 7.5e3 / 25),
        models.Population("", "", 24e3 / 25),
        population_sample_0,
        models.Population("", "", 45e3 / 25),
        models.Population("", "", 50e3 / 25),
        population_sample_none
    ]

    def __init__(self):

        generation_time = 25

        # population sizes
        N_Losch = 1.92e3
        N_Mbu = 1.73e4
        N_Mbu_Losch = 2.91e4
        N_Han = 6.3e3
        N_Han_Losch = 2.34e3
        N_Nean_Losch = 1.82e4
        N_Nean = 86.9
        N_LBK = 75.7
        N_Sard = 1.5e4
        N_Sard_LBK = 1.2e4

        # unknown population sizes
        # these should be set to ancestral Eurasian population size,
        # but at the moment it's unclear to me if that is N_Han_Losch? N_Losch?
        N_Basal = N_Losch
        N_Ust = N_Basal
        N_MA1 = N_Basal

        # population merge times in years, divided by generation time
        t_Mbu_Losch = 9.58e4 / generation_time
        t_Han_Losch = 5.04e4 / generation_time
        t_Ust_Losch = 5.15e4 / generation_time
        t_Nean_Losch = 6.96e5 / generation_time
        t_MA1_Losch = 4.49e4 / generation_time
        t_LBK_Losch = 3.77e4 / generation_time
        t_Basal_Losch = 7.98e4 / generation_time
        t_Sard_LBK = 7.69e3 / generation_time
        # t_GhostWHG_Losch = 1.56e3 / generation_time

        # pulse admixture times and fractions
        p_Nean_to_Eur = 0.0296
        t_Nean_to_Eur = 5.68e4 / generation_time
        p_Basal_to_EEF = 0.0936
        t_Basal_to_EEF = 3.37e4 / generation_time
        p_GhostWHG_to_Sard = 0.0317
        t_GhostWHG_to_Sard = 1.23e3 / generation_time

        # sample_times (in years), divided by estimated generation time
        t_Mbuti = 0
        t_Han = 0
        t_Sardinian = 0
        t_Loschbour = 7.5e3 / generation_time
        t_LBK = 8e3 / generation_time
        t_MA1 = 24e3 / generation_time
        t_UstIshim = 45e3 / generation_time
        t_Altai = 50e3 / generation_time

        # set up populations
        self.population_configurations = [
            msprime.PopulationConfiguration(  # Mbuti
                initial_size=N_Mbu, growth_rate=0,
                metadata={"name": "Mbuti", "sampling_time": t_Mbuti}),
            msprime.PopulationConfiguration(  # LBK
                initial_size=N_LBK, growth_rate=0,
                metadata={"name": "LBK", "sampling_time": t_LBK}),
            msprime.PopulationConfiguration(  # Sardinian
                initial_size=N_Sard, growth_rate=0,
                metadata={"name": "Sardinian", "sampling_time": t_Sardinian}),
            msprime.PopulationConfiguration(  # Loschbour
                initial_size=N_Losch, growth_rate=0,
                metadata={"name": "Loschbour", "sampling_time": t_Loschbour}),
            msprime.PopulationConfiguration(  # MA1
                initial_size=N_MA1, growth_rate=0,
                metadata={"name": "MA1", "sampling_time": t_MA1}),
            msprime.PopulationConfiguration(  # Han
                initial_size=N_Han, growth_rate=0,
                metadata={"name": "Han", "sampling_time": t_Han}),
            msprime.PopulationConfiguration(  # UstIshim
                initial_size=N_Ust, growth_rate=0,
                metadata={"name": "UstIshim", "sampling_time": t_UstIshim}),
            msprime.PopulationConfiguration(  # Neanderthal
                initial_size=N_Nean, growth_rate=0,
                metadata={"name": "Altai", "sampling_time": t_Altai}),
            msprime.PopulationConfiguration(  # Basal Eurasian
                initial_size=N_Basal, growth_rate=0,
                metadata={"name": "Basal", "sampling_time": None})
        ]

        # no migration rates, only pulse events, so set mig mat to zeros
        num_pops = len(self.population_configurations)
        self.migration_matrix = [[0] * num_pops] * num_pops

        # Compute Neanderthal pop size decline rate
        # I'm assuming that the N_Nean is the size of Neanderthal population
        # at the time of sampling the Altai individual
        r_Nean = -math.log(N_Nean_Losch/N_Nean) / (t_Mbu_Losch-t_Altai)

        # Using columns in figure in Kamm paper as proxies for pop number
        self.demographic_events = [
            msprime.MassMigration(
                time=t_GhostWHG_to_Sard, source=2,
                destination=3, proportion=p_GhostWHG_to_Sard),
            msprime.MassMigration(
                time=t_Sard_LBK, source=2, destination=1,
                proportion=1.),
            msprime.PopulationParametersChange(
                time=t_Sard_LBK, initial_size=N_Sard_LBK,
                population_id=1),
            msprime.MassMigration(
                time=t_Basal_to_EEF, source=1, destination=8,
                proportion=p_Basal_to_EEF),
            msprime.MassMigration(
                time=t_LBK_Losch, source=1, destination=3,
                proportion=1.),
            msprime.MassMigration(
                time=t_MA1_Losch, source=4, destination=3,
                proportion=1.),
            msprime.PopulationParametersChange(
                time=t_Altai, initial_size=N_Nean,
                growth_rate=r_Nean, population_id=7),
            msprime.MassMigration(
                time=t_Han_Losch, source=5, destination=3,
                proportion=1.),
            msprime.PopulationParametersChange(
                time=t_Han_Losch, initial_size=N_Han_Losch,
                population_id=3),
            msprime.MassMigration(
                time=t_Ust_Losch, source=6, destination=3,
                proportion=1.),
            msprime.MassMigration(
                time=t_Nean_to_Eur, source=3, destination=7,
                proportion=p_Nean_to_Eur),
            msprime.MassMigration(
                time=t_Basal_Losch, source=8, destination=3,
                proportion=1.),
            msprime.MassMigration(
                time=t_Mbu_Losch, source=0, destination=3,
                proportion=1.),
            msprime.PopulationParametersChange(
                time=t_Mbu_Losch, initial_size=N_Mbu_Losch,
                population_id=3),
            msprime.PopulationParametersChange(
                time=t_Mbu_Losch, initial_size=N_Nean_Losch,
                growth_rate=0, population_id=7),
            msprime.MassMigration(
                time=t_Nean_Losch, source=7, destination=3,
                proportion=1.),
            msprime.PopulationParametersChange(
                time=t_Nean_Losch, initial_size=N_Nean_Losch,
                population_id=3)
        ]


_species.get_demographic_model(
        "AncientEurasia_9K19").register_qc(KammAncientSamples())


class DenisovanAncestryInPapuans(models.DemographicModel):
    """
    Demographic model from Jacobs et al (2019). The model is
    illustrated on Figure S5, parameters are in Table S5
    """

    populations = [population_sample_0] * 4 + [
        models.Population("", "", 2058),
        models.Population("", "", 2612),
    ] + [population_sample_none] * 4

    def __init__(self):
        self.generation_time = 29  # Just information

        # sizes of populations
        N_YRI = 48433
        N_CEU = 6962
        N_CHB = 9025
        N_Papuan = 8834
        N_DenA = 5083
        t_DenA = 2058  # generations
        N_NeanA = 826
        t_NeanA = 2612  # generations
        N_Nean1 = 13249
        N_Den1 = N_Nean1
        N_Den2 = N_Nean1
        N_Ghost = 8516
        # after coalescences
        N_CEU_CHB = 12971
        N_Human = 41563
        N_DenAnc = 100
        N_A = 32671

        # bottlenecks
        N_CEU_CHB_bot = 2231
        N_GhostA_bot = 1394
        N_Papuan_bot = 243

        # times of coalesces
        t_CEU_CHB = 1293
        t_CEU_Ghost = 1758
        t_Papuan_Ghost = 1784
        t_YRI_GhostA = 2218
        t_Nean1_NeanA = 3375
        t_Den1_DenA = 9750
        t_Den1_Den2 = 12500
        t_Den_Nean = 15090
        t_Human_Den_Nean = 20225

        # times of bottlenecks
        t_CEU_CHB_bot = 1659
        t_Papuan_bot = 1685
        t_GhostA_bot = 2119

        # migrations
        m_YRI_Ghost = 0.000179
        m_Ghost_CEU = 0.000442
        m_CEU_CHB = 3.14e-5
        m_CHB_Papuan = 5.72e-5
        m_CEUCHB_Papua = 0.000572
        m_Ghost_CEUCHB = 0.000442

        # times and proportions of admixtures
        p1 = 0.55
        t_Nean1_to_CHB = 883
        p_Nean1_to_CHB = 0.002
        t_Den2_to_Papuan = 45.7e3 / self.generation_time
        p_Den2_to_Papuan = (1 - p1) * 0.04
        t_Den1_to_Papuan = 29.8e3 / self.generation_time
        p_Den1_to_Papuan = p1 * 0.04
        t_Nean1_to_Papuan = 1412
        p_Nean1_to_Papuan = 0.002
        t_Nean1_to_CEU_CHB = 1566
        p_Nean1_to_CEU_CHB = 0.011
        t_Nean1_to_GhostA = 1853
        p_Nean1_to_GhostA = 0.024

        # set up populations
        self.population_configurations = [
            msprime.PopulationConfiguration(  # 0 YRI
                initial_size=N_YRI, growth_rate=0,
                metadata={"name": "YRI", "sampling_time": 0}),
            msprime.PopulationConfiguration(  # 1 CEU
                initial_size=N_CEU, growth_rate=0,
                metadata={"name": "CEU", "sampling_time": 0}),
            msprime.PopulationConfiguration(  # 2 CHB
                initial_size=N_CHB, growth_rate=0,
                metadata={"name": "CHB", "sampling_time": 0}),
            msprime.PopulationConfiguration(  # 3 Papuan
                initial_size=N_Papuan, growth_rate=0,
                metadata={"name": "Papuan", "sampling_time": 0}),
            msprime.PopulationConfiguration(  # 4 DenA
                initial_size=N_DenA, growth_rate=0,
                metadata={"name": "DenA", "sampling_time": t_DenA}),
            msprime.PopulationConfiguration(  # 5 NeanA
                initial_size=N_NeanA, growth_rate=0,
                metadata={"name": "NeanA", "sampling_time": t_NeanA}),
            msprime.PopulationConfiguration(  # 6 Den1
                initial_size=N_Den1, growth_rate=0,
                metadata={"name": "Den1", "sampling_time": None}),
            msprime.PopulationConfiguration(  # 7 Den2
                initial_size=N_Den2, growth_rate=0,
                metadata={"name": "Den2", "sampling_time": None}),
            msprime.PopulationConfiguration(  # 8 Nean1
                initial_size=N_Nean1, growth_rate=0,
                metadata={"name": "Nean1", "sampling_time": None}),
            msprime.PopulationConfiguration(  # 9 Ghost
                initial_size=N_Ghost, growth_rate=0,
                metadata={"name": "Ghost", "sampling_time": None})
        ]

        self.migration_matrix = [[0]*10 for _ in range(10)]
        self.migration_matrix[0][9] = m_YRI_Ghost
        self.migration_matrix[9][0] = m_YRI_Ghost
        self.migration_matrix[1][9] = m_Ghost_CEU
        self.migration_matrix[9][1] = m_Ghost_CEU
        self.migration_matrix[1][2] = m_CEU_CHB
        self.migration_matrix[2][1] = m_CEU_CHB
        self.migration_matrix[2][3] = m_CHB_Papuan
        self.migration_matrix[3][2] = m_CHB_Papuan

        self.demographic_events = [
            # Coalescence of CEU and CHB into CHB
            msprime.MassMigration(
                time=t_CEU_CHB, source=1,
                destination=2, proportion=1.),
            # Set size of CEU+CHB population
            msprime.PopulationParametersChange(
                time=t_CEU_CHB, initial_size=N_CEU_CHB,
                population_id=2),
            # Change migration matrix
            msprime.MigrationRateChange(
                time=t_CEU_CHB, rate=0, matrix_index=(2, 1)),
            msprime.MigrationRateChange(
                time=t_CEU_CHB, rate=0, matrix_index=(1, 2)),
            msprime.MigrationRateChange(
                time=t_CEU_CHB, rate=0, matrix_index=(3, 2)),
            msprime.MigrationRateChange(
                time=t_CEU_CHB, rate=0, matrix_index=(2, 3)),
            msprime.MigrationRateChange(
                time=t_CEU_CHB, rate=0, matrix_index=(9, 1)),
            msprime.MigrationRateChange(
                time=t_CEU_CHB, rate=0, matrix_index=(1, 9)),
            msprime.MigrationRateChange(
                time=t_CEU_CHB, rate=m_CEUCHB_Papua, matrix_index=(2, 3)),
            msprime.MigrationRateChange(
                time=t_CEU_CHB, rate=m_CEUCHB_Papua, matrix_index=(3, 2)),
            msprime.MigrationRateChange(
                time=t_CEU_CHB, rate=m_Ghost_CEUCHB, matrix_index=(2, 9)),
            msprime.MigrationRateChange(
                time=t_CEU_CHB, rate=m_Ghost_CEUCHB, matrix_index=(9, 2)),
            # Set bottleneck size of CEU+CHB population
            msprime.PopulationParametersChange(
                time=t_CEU_CHB_bot, initial_size=N_CEU_CHB_bot,
                population_id=2),
            # Change migration matrix
            msprime.MigrationRateChange(time=t_CEU_CHB_bot, rate=0),
            # Set bottleneck size of Papuan population
            msprime.PopulationParametersChange(
                time=t_Papuan_bot, initial_size=N_Papuan_bot,
                population_id=3),
            # Coalescence of CEU+CHB and Ghost to Ghost
            msprime.MassMigration(
                time=t_CEU_Ghost, source=2,
                destination=9, proportion=1.),
            # Coalescence of Papuan and Ghost to GhostA
            msprime.MassMigration(
                time=t_Papuan_Ghost, source=3,
                destination=9, proportion=1.),
            # Set bottleneck size of GhostA population
            msprime.PopulationParametersChange(
                time=t_GhostA_bot, initial_size=N_GhostA_bot,
                population_id=9),
            # Coalescence of Ghost and YRI into Human (YRI)
            msprime.MassMigration(
                time=t_YRI_GhostA, source=9,
                destination=0, proportion=1.),
            # Set size of Human population
            msprime.PopulationParametersChange(
                time=t_YRI_GhostA, initial_size=N_Human,
                population_id=0),
            # Coalescence of NeanA and Nean1 into NeanAnc
            msprime.MassMigration(
                time=t_Nean1_NeanA, source=8,
                destination=5, proportion=1.),
            # Set size of NeanAnc population
            msprime.PopulationParametersChange(
                time=t_Nean1_NeanA, initial_size=N_Nean1,
                population_id=5),
            # Coalescence of Den1 and DenA into DenAnc (DenA)
            msprime.MassMigration(
                time=t_Den1_DenA, source=6,
                destination=4, proportion=1.),
            # Set size of DenA population
            msprime.PopulationParametersChange(
                time=t_Den1_DenA, initial_size=N_DenAnc,
                population_id=4),
            # Coalescence of DenAnc and Den2 into DenAnc
            msprime.MassMigration(
                time=t_Den1_Den2, source=7,
                destination=4, proportion=1.),
            # Set size of DenA population
            msprime.PopulationParametersChange(
                time=t_Den1_Den2, initial_size=N_DenAnc,
                population_id=4),
            # Coalescence of DenAnc and NeanAnc into Den_Nean (Nean1)
            msprime.MassMigration(
                time=t_Den_Nean, source=5,
                destination=4, proportion=1.),
            # Set size of Den_Nean population
            msprime.PopulationParametersChange(
                time=t_Den_Nean, initial_size=N_Nean1,
                population_id=4),
            # Coalescence of Den_Nean and Human into Anc (YRI)
            msprime.MassMigration(
                time=t_Human_Den_Nean, source=4,
                destination=0, proportion=1.),
            # Set ancestral size of population
            msprime.PopulationParametersChange(
                time=t_Human_Den_Nean, initial_size=N_A,
                population_id=0),

            # Admixture events
            # Admixture from Den1 to Papuans
            msprime.MassMigration(
                time=t_Den1_to_Papuan, source=3,
                destination=6, proportion=p_Den1_to_Papuan),
            # Admixture from Den2 to Papuans
            msprime.MassMigration(
                time=t_Den2_to_Papuan, source=3,
                destination=7, proportion=p_Den2_to_Papuan),
            # Admixture from Nean1 to GhostA
            msprime.MassMigration(
                time=t_Nean1_to_GhostA, source=9,
                destination=8, proportion=p_Nean1_to_GhostA),
            # Admixture from Nean1 to CEU+CHB
            msprime.MassMigration(
                time=t_Nean1_to_CEU_CHB, source=2,
                destination=8, proportion=p_Nean1_to_CEU_CHB),
            # Admixture from Nean1 to Papuans
            msprime.MassMigration(
                time=t_Nean1_to_Papuan, source=3,
                destination=8, proportion=p_Nean1_to_Papuan),
            # Admixture from Neandertal to East Asia population
            msprime.MassMigration(
                time=t_Nean1_to_CHB, proportion=p_Nean1_to_CHB, source=2,
                destination=8),
        ]
        self.demographic_events.sort(key=lambda x: x.time)


_species.get_demographic_model(
        "PapuansOutOfAfrica_10J19").register_qc(DenisovanAncestryInPapuans())


class GutenkunstOOA(models.DemographicModel):
    """
    From Gutenkunst et al (2009). Parameters taken from maximum likelihood values
    in Table 1 in that paper.
    """
    populations = [population_sample_0] * 3

    def __init__(self):
        generation_time = 25

        # Population sizes
        N_A = 7300
        N_AF = 12300
        N_B = 2100
        N_EU0 = 1000
        N_AS0 = 510

        # Growth rates per generation
        r_EU = 0.4e-2
        r_AS = 0.55e-2

        # Migration rates
        m_AF_B = 25e-5
        m_AF_EU = 3e-5
        m_AF_AS = 1.9e-5
        m_EU_AS = 9.6e-5

        # Epoch times
        T_AF = 220e3/generation_time
        T_B = 140e3/generation_time
        T_EU_AS = 21.2e3/generation_time

        # Calculate population sizes at modern (T=0) time
        N_EUF = N_EU0 * math.exp(r_EU * T_EU_AS)
        N_ASF = N_AS0 * math.exp(r_AS * T_EU_AS)

        self.population_configurations = [
            msprime.PopulationConfiguration(
                initial_size=N_AF, growth_rate=0),
            msprime.PopulationConfiguration(
                initial_size=N_EUF, growth_rate=r_EU),
            msprime.PopulationConfiguration(
                initial_size=N_ASF, growth_rate=r_AS)
        ]

        # Setup initial migration matrix
        self.migration_matrix = [
            [0, m_AF_EU, m_AF_AS],
            [m_AF_EU, 0, m_EU_AS],
            [m_AF_AS, m_EU_AS, 0]
        ]

        self.demographic_events = [
            # CEU and CHB merge into B, reset migration rates to Af-B, change pop size
            msprime.MassMigration(time=T_EU_AS, source=2, dest=1, proportion=1),
            msprime.MigrationRateChange(time=T_EU_AS, rate=0),
            msprime.MigrationRateChange(time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
            msprime.MigrationRateChange(time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
            msprime.PopulationParametersChange(time=T_EU_AS, initial_size=N_B,
                                               population_id=1, growth_rate=0),
            # B and AF merge, turn migration off, reset population size
            msprime.MassMigration(time=T_B, source=1, dest=0, proportion=1),
            msprime.MigrationRateChange(time=T_B, rate=0),
            # Ancestral size change, reset population size
            msprime.PopulationParametersChange(time=T_AF, initial_size=N_A,
                                               population_id=0)
        ]


_species.get_demographic_model("OutOfAfrica_3G09").register_qc(GutenkunstOOA())


class ZigZag(models.DemographicModel):
    """
    Model from Schiffels et al (2014) used to test inference methods on a "zigzag"
    demography. Single population, repeated growth and decay of pop size.
    """
    populations = [population_sample_0]

    def __init__(self):
        # from section 7 of the supplement
        self.generation_time = 30
        mu = 1.25e-8
        L = 10000000
        theta = 7156
        Ne = theta / 4 / mu / L

        self.population_configurations = [
            msprime.PopulationConfiguration(
                initial_size=Ne, growth_rate=0,
                metadata={"sampling_time": 0})
        ]

        rate_changes = {
            0.000582262: 1318.18,
            0.00232905: -329.546,
            0.00931619: 82.3865,
            0.0372648: -20.5966,
            0.149059: 5.14916,
            # 0.596236: 0
        }

        de = []
        for time, rate in rate_changes.items():
            de.append(
                msprime.PopulationParametersChange(
                    time=time * 4 * Ne,
                    growth_rate=rate / 4 / Ne,
                    population_id=0
                )
            )

        de.append(
            msprime.PopulationParametersChange(
                time=0.596236 * 4 * Ne,
                growth_rate=0,
                population_id=0,
                initial_size=0.1 * Ne
                )
            )

        self.demographic_events = de

        self.migration_matrix = [[0]]


_species.get_demographic_model("Zigzag_1S14").register_qc(ZigZag())
