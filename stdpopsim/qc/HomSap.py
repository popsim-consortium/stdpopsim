import math

import msprime

import stdpopsim


_species = stdpopsim.get_species("HomSap")

# Some generic populations to use for qc
population_sample_0 = stdpopsim.Population(
    "qc_sampling_0", "Population that samples at time 0", 0
)
population_sample_none = stdpopsim.Population(
    "qc_sampling_none", "Population that does not sample", None
)


def TennessenOnePopAfrica():
    # This model is the same as the Tennessen two population model except
    # the European population has been removed.
    id = "QC-Africa_1T12"
    populations = [population_sample_0]

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

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=generation_time,
        populations=populations,
        # Now we set up the population configurations. The population IDs are 0=YRI. This
        # includes both the inital sizes, growth rates, and migration rates.
        population_configurations=[
            msprime.PopulationConfiguration(initial_size=N_AF1, growth_rate=r_AF0),
        ],
        # Now we add the demographic events working backwards in time. Starting with the
        # growth slowdown in Europeans and the transition to a fixed population size in
        # Africans.
        demographic_events=[
            # Reversion to fixed population size in Africans
            msprime.PopulationParametersChange(
                time=T_AG, initial_size=N_AF0, growth_rate=0, population_id=0
            ),
            # Change to ancestral population size pre OOA
            msprime.PopulationParametersChange(
                time=T_AF, initial_size=N_A, population_id=0
            ),
        ],
        # Mapping of per-epoch population names in the production model to
        # population IDs in this model.
        population_id_map=[
            {"AFR": 0},
            {"AFR": 0},
            {"AFR": 0},
        ],
        # I'm guessing the took the mutation rate to be the same as in Gravel et al 2011,
        # (doi: 10.1073/pnas.1019276108), where they got the demographic model from
        mutation_rate=2.36e-8,
    )


_species.get_demographic_model("Africa_1T12").register_qc(TennessenOnePopAfrica())


def TennessenTwoPopOutOfAfrica():
    id = "QC-OutOfAfrica_2T12"
    populations = [population_sample_0] * 2

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
    N_EU1 = N_EU0 * math.exp(r_EU0 * (T_EU_AS - T_AG))
    # European population size after 2nd expansion
    N_EU2 = N_EU1 * math.exp(r_EU1 * T_AG)
    # African population size after 1st expansion
    N_AF1 = N_AF0 * math.exp(r_AF0 * T_AG)

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=generation_time,
        populations=populations,
        # Now we set up the population configurations. The population IDs are
        # 0=CEU and 1=YRI. This includes both the inital sizes, growth rates,
        # and migration rates.
        population_configurations=[
            msprime.PopulationConfiguration(initial_size=N_AF1, growth_rate=r_AF0),
            msprime.PopulationConfiguration(initial_size=N_EU2, growth_rate=r_EU1),
        ],
        migration_matrix=[
            [0, m_AF1_EU1],
            [m_AF1_EU1, 0],
        ],
        # Now we add the demographic events working backwards in time. Starting
        # with the growth slowdown in Europeans and the transition to a fixed
        # population size in Africans.
        demographic_events=[
            # Set the migration rate for 1st CEU growth period (for now stays same)
            msprime.MigrationRateChange(time=T_AG, rate=m_AF1_EU1, matrix_index=(0, 1)),
            msprime.MigrationRateChange(time=T_AG, rate=m_AF1_EU1, matrix_index=(1, 0)),
            # Growth slowdown in Europeans
            msprime.PopulationParametersChange(
                time=T_AG, initial_size=N_EU1, growth_rate=r_EU0, population_id=1
            ),
            # Reversion to fixed population size in Africans
            msprime.PopulationParametersChange(
                time=T_AG, initial_size=N_AF0, growth_rate=0, population_id=0
            ),
            # Set the migration rate for pre CEU/CHB split
            msprime.MigrationRateChange(
                time=T_EU_AS, rate=m_AF0_B, matrix_index=(0, 1)
            ),
            msprime.MigrationRateChange(
                time=T_EU_AS, rate=m_AF0_B, matrix_index=(1, 0)
            ),
            # Reversion to fixed population size at the time of the CHB/CEU split
            msprime.PopulationParametersChange(
                time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1
            ),
            # Coalescence between the OOA and YRI pops
            msprime.MassMigration(time=T_B, source=1, destination=0, proportion=1.0),
            msprime.MigrationRateChange(time=T_B, rate=0),
            # Change to ancestral population size pre OOA
            msprime.PopulationParametersChange(
                time=T_AF, initial_size=N_A, population_id=0
            ),
        ],
        # Mapping of per-epoch population names in the production model to
        # population IDs in this model.
        population_id_map=[
            {"AFR": 0, "EUR": 1},
            {"AFR": 0, "EUR": 1},
            {"AFR": 0, "EUR": 1},
            {"AFR": 0, "EUR": 1},
            {"AFR": 0, "EUR": 1},
        ],
        # I'm guessing the took the mutation rate to be the same as in Gravel et al 2011,
        # (doi: 10.1073/pnas.1019276108), where they got the demographic model from
        mutation_rate=2.36e-8,
    )


_species.get_demographic_model("OutOfAfrica_2T12").register_qc(
    TennessenTwoPopOutOfAfrica()
)


def BrowningAmerica():
    id = "QC-AmericanAdmixture_4B11"
    populations = [population_sample_0] * 4

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
    mm_AF1 = 1 / 6
    # Adjusted fraction for remaining population after AF migration (5/6 * 2/5 = 1/3)
    mm_CEU0 = 2 / 5
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

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=25,
        populations=populations,
        # Set population sizes at T=0
        # pop0 is Africa, pop1 is Europe, pop2 is Asia, pop3 is admixed
        population_configurations=[
            msprime.PopulationConfiguration(initial_size=N_AF1, growth_rate=0),
            msprime.PopulationConfiguration(initial_size=N_CEU1, growth_rate=r_CEU0),
            msprime.PopulationConfiguration(initial_size=N_CHB1, growth_rate=r_CHB0),
            msprime.PopulationConfiguration(
                initial_size=N_ADMIX1, growth_rate=r_ADMIX0
            ),
        ],
        # Migration matrix, all migrations to admixed population are 0
        migration_matrix=[
            [0, m_AF1_CEU0, m_AF1_CHB0, 0],
            [m_AF1_CEU0, 0, m_CEU0_CHB0, 0],
            [m_AF1_CHB0, m_CEU0_CHB0, 0, 0],
            [0, 0, 0, 0],
        ],
        # Now we add the demographic events working backwards in time.
        demographic_events=[
            # Admixed population recoalesces with origin populations (T_ADMIX0)
            msprime.MassMigration(
                time=T_ADMIX0, source=3, destination=0, proportion=mm_AF1
            ),
            msprime.MassMigration(
                time=T_ADMIX0, source=3, destination=1, proportion=mm_CEU0
            ),
            msprime.MassMigration(
                time=T_ADMIX0, source=3, destination=2, proportion=mm_CHB0
            ),
            # Zero out migration rate (desn't matter but added for equality to prod.)
            msprime.MigrationRateChange(time=T_CEU_CHB, rate=0.0),
            # CEU and CHB coalesce and set population to OOA size (T_CEU_CHB)
            msprime.MassMigration(
                time=T_CEU_CHB, source=2, destination=1, proportion=1.0
            ),
            msprime.PopulationParametersChange(
                time=T_CEU_CHB,
                initial_size=N_OOA,
                growth_rate=0.0,
                population_id=1,
            ),
            # Set OOA <--> AF migration rate (T_CEU_CHB)
            msprime.MigrationRateChange(
                time=T_CEU_CHB, rate=m_AF1_OOA, matrix_index=(0, 1)
            ),
            msprime.MigrationRateChange(
                time=T_CEU_CHB, rate=m_AF1_OOA, matrix_index=(1, 0)
            ),
            # Zero out migration rate
            msprime.MigrationRateChange(time=T_AF1_OOA, rate=0.0),
            # OOA and AF1 coalesce (T_OOA)
            msprime.MassMigration(
                time=T_AF1_OOA, source=1, destination=0, proportion=1.0
            ),
            # AF1 -> AF0 population size change (T_AF0_AF1)
            msprime.PopulationParametersChange(
                time=T_AF0_AF1, initial_size=N_AF0, population_id=0
            ),
        ],
        # Mapping of per-epoch population names in the production model to
        # population IDs in this model.
        population_id_map=[
            {"AFR": 0, "EUR": 1, "ASIA": 2, "ADMIX": 3},
            {"AFR": 0, "EUR": 1, "ASIA": 2, "ADMIX": 3},
            {"AFR": 0, "EUR": 1, "ASIA": 2, "ADMIX": 3},
            {"AFR": 0, "EUR": 1, "ASIA": 2, "ADMIX": 3},
            {"AFR": 0, "EUR": 1, "ASIA": 2, "ADMIX": 3},
        ],
        # the model derives from Gravel et al, who used this rate
        mutation_rate=2.36e-8,
    )


_species.get_demographic_model("AmericanAdmixture_4B11").register_qc(BrowningAmerica())


def RagsdaleArchaic():
    id = "QC-OutOfAfricaArchaicAdmixture_5R19"
    populations = [population_sample_0] * 3 + [population_sample_none] * 2

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
    T_AF = 300e3 / generation_time
    T_OOA = 60.7e3 / generation_time
    T_CEU_CHB = 36e3 / generation_time
    T_ARCHAF_split = 499e3 / generation_time
    T_ARCHAF_mig = 125e3 / generation_time
    T_NEAN_split = 559e3 / generation_time
    T_ARCH_ADMIX_end = 18.7e3 / generation_time

    # Calculate population sizes at modern (T=0) time
    N_CEU1 = N_CEU0 * math.exp(r_CEU * T_CEU_CHB)
    N_CHB1 = N_CHB0 * math.exp(r_CHB * T_CEU_CHB)

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=generation_time,
        populations=populations,
        # Set population sizes at T=0
        # pop0 is Africa, pop1 is Europe, pop2 is Asia, pop3 is Neanderthal, pop4 is
        # archaic african
        population_configurations=[
            msprime.PopulationConfiguration(initial_size=N_YRI, growth_rate=0),
            msprime.PopulationConfiguration(initial_size=N_CEU1, growth_rate=r_CEU),
            msprime.PopulationConfiguration(initial_size=N_CHB1, growth_rate=r_CHB),
            msprime.PopulationConfiguration(initial_size=N_0, growth_rate=0),
            msprime.PopulationConfiguration(initial_size=N_0, growth_rate=0),
        ],
        # Setup initial migration matrix
        migration_matrix=[
            [0, m_YRI_CEU, m_YRI_CHB, 0, 0],  # noqa
            [m_YRI_CEU, 0, m_CEU_CHB, 0, 0],  # noqa
            [m_YRI_CHB, m_CEU_CHB, 0, 0, 0],  # noqa
            [0, 0, 0, 0, 0],  # noqa
            [0, 0, 0, 0, 0],  # noqa
        ],
        demographic_events=[
            # Migration between YRI and ARCHAF(E1)
            msprime.MigrationRateChange(
                time=T_ARCH_ADMIX_end, rate=m_AF_ARCHAF, matrix_index=(0, 4)
            ),
            msprime.MigrationRateChange(
                time=T_ARCH_ADMIX_end, rate=m_AF_ARCHAF, matrix_index=(4, 0)
            ),
            # Migration between CEU and NEAN(E1)
            msprime.MigrationRateChange(
                time=T_ARCH_ADMIX_end, rate=m_OOA_NEAN, matrix_index=(1, 3)
            ),
            msprime.MigrationRateChange(
                time=T_ARCH_ADMIX_end, rate=m_OOA_NEAN, matrix_index=(3, 1)
            ),
            # Migration between CHB and NEAN(E1)
            msprime.MigrationRateChange(
                time=T_ARCH_ADMIX_end, rate=m_OOA_NEAN, matrix_index=(2, 3)
            ),
            msprime.MigrationRateChange(
                time=T_ARCH_ADMIX_end, rate=m_OOA_NEAN, matrix_index=(3, 2)
            ),
            # Coalescence of CHB into CEU (E2)
            msprime.MassMigration(time=T_CEU_CHB, source=2, dest=1, proportion=1.0),
            # Reset migration rates (E2)(redundant)*
            msprime.MigrationRateChange(time=T_CEU_CHB, rate=0.0),
            # Migration rate change between OOA(CEU) and AF(YRI)(E2)
            msprime.MigrationRateChange(
                time=T_CEU_CHB, rate=m_AF_B, matrix_index=(0, 1)
            ),
            msprime.MigrationRateChange(
                time=T_CEU_CHB, rate=m_AF_B, matrix_index=(1, 0)
            ),
            # Migration between YRI and ARCHAF (E2)(redundant without mig. rate reset)*
            msprime.MigrationRateChange(
                time=T_CEU_CHB, rate=m_AF_ARCHAF, matrix_index=(0, 4)
            ),
            msprime.MigrationRateChange(
                time=T_CEU_CHB, rate=m_AF_ARCHAF, matrix_index=(4, 0)
            ),
            # Migration between CEU and NEAN (E2)(redundant without mig. rate reset)*
            msprime.MigrationRateChange(
                time=T_CEU_CHB, rate=m_OOA_NEAN, matrix_index=(1, 3)
            ),
            msprime.MigrationRateChange(
                time=T_CEU_CHB, rate=m_OOA_NEAN, matrix_index=(3, 1)
            ),
            # CEU change to fixed population size at the time of the CHB/CEU coal. (E2)
            msprime.PopulationParametersChange(
                time=T_CEU_CHB, initial_size=N_B, growth_rate=0, population_id=1
            ),
            # Coalescence between the OOA and AF pops (E3)
            msprime.MassMigration(time=T_OOA, source=1, destination=0, proportion=1.0),
            # Reset migration rates (E3)
            msprime.MigrationRateChange(time=T_OOA, rate=0.0),
            # Migration between YRI and ARCHAF (E3)
            msprime.MigrationRateChange(
                time=T_OOA, rate=m_AF_ARCHAF, matrix_index=(0, 4)
            ),
            msprime.MigrationRateChange(
                time=T_OOA, rate=m_AF_ARCHAF, matrix_index=(4, 0)
            ),
            # Migration between archaic african and african pop. "ends" (E4)
            msprime.MigrationRateChange(time=T_ARCHAF_mig, rate=0),
            # AF reverts to ancestral population size pre OOA (E5)
            msprime.PopulationParametersChange(
                time=T_AF, initial_size=N_0, population_id=0
            ),
            # Archaic AF population coalesces into AF (E6)
            msprime.MassMigration(
                time=T_ARCHAF_split, source=4, dest=0, proportion=1.0
            ),
            # NEAN pop. coalesces into AF (E7)
            msprime.MassMigration(time=T_NEAN_split, source=3, dest=0, proportion=1.0),
        ],
        # Mapping of per-epoch population names in the production model to
        # population IDs in this model.
        population_id_map=[
            {"YRI": 0, "CEU": 1, "CHB": 2, "Neanderthal": 3, "ArchaicAFR": 4},
            {"YRI": 0, "CEU": 1, "CHB": 2, "Neanderthal": 3, "ArchaicAFR": 4},
            {"YRI": 0, "CEU": 1, "CHB": 2, "Neanderthal": 3, "ArchaicAFR": 4},
            {"YRI": 0, "CEU": 1, "CHB": 2, "Neanderthal": 3, "ArchaicAFR": 4},
            {"YRI": 0, "CEU": 1, "CHB": 2, "Neanderthal": 3, "ArchaicAFR": 4},
            {"YRI": 0, "CEU": 1, "CHB": 2, "Neanderthal": 3, "ArchaicAFR": 4},
            {"YRI": 0, "CEU": 1, "CHB": 2, "Neanderthal": 3, "ArchaicAFR": 4},
            {"YRI": 0, "CEU": 1, "CHB": 2, "Neanderthal": 3, "ArchaicAFR": 4},
        ],
        # the method depends on recombination rate, not mutation rate: "This
        # normalization removes all dependence of the statistics on the overall
        # mutation rate, so that estimates of split times and population sizes
        # are calibrated by the recombination rate per generation instead of
        # the mutation rate"
        mutation_rate=None,
    )


_species.get_demographic_model("OutOfAfricaArchaicAdmixture_5R19").register_qc(
    RagsdaleArchaic()
)


def KammAncientSamples():
    """
    Demographic inferred by momi described in Kamm et al. (2019). The model is
    illustrated in Figure 3, with parameters given in Table 2.
    """
    id = "AncientEurasia_9K19"
    populations = [
        population_sample_0,
        stdpopsim.Population("", "", 8e3 / 25),
        population_sample_0,
        stdpopsim.Population("", "", 7.5e3 / 25),
        stdpopsim.Population("", "", 24e3 / 25),
        population_sample_0,
        stdpopsim.Population("", "", 45e3 / 25),
        stdpopsim.Population("", "", 50e3 / 25),
        population_sample_none,
    ]

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

    # Compute Neanderthal pop size decline rate
    # I'm assuming that the N_Nean is the size of Neanderthal population
    # at the time of sampling the Altai individual
    r_Nean = -math.log(N_Nean_Losch / N_Nean) / (t_Mbu_Losch - t_Altai)

    pop_id_map = {
        "Mbuti": 0,
        "LBK": 1,
        "Sardinian": 2,
        "Loschbour": 3,
        "MA1": 4,
        "Han": 5,
        "UstIshim": 6,
        "Neanderthal": 7,
        "BasalEurasian": 8,
    }

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=generation_time,
        populations=populations,
        # set up populations
        population_configurations=[
            msprime.PopulationConfiguration(  # Mbuti
                initial_size=N_Mbu,
                growth_rate=0,
                metadata={"name": "Mbuti", "sampling_time": t_Mbuti},
            ),
            msprime.PopulationConfiguration(  # LBK
                initial_size=N_LBK,
                growth_rate=0,
                metadata={"name": "LBK", "sampling_time": t_LBK},
            ),
            msprime.PopulationConfiguration(  # Sardinian
                initial_size=N_Sard,
                growth_rate=0,
                metadata={"name": "Sardinian", "sampling_time": t_Sardinian},
            ),
            msprime.PopulationConfiguration(  # Loschbour
                initial_size=N_Losch,
                growth_rate=0,
                metadata={"name": "Loschbour", "sampling_time": t_Loschbour},
            ),
            msprime.PopulationConfiguration(  # MA1
                initial_size=N_MA1,
                growth_rate=0,
                metadata={"name": "MA1", "sampling_time": t_MA1},
            ),
            msprime.PopulationConfiguration(  # Han
                initial_size=N_Han,
                growth_rate=0,
                metadata={"name": "Han", "sampling_time": t_Han},
            ),
            msprime.PopulationConfiguration(  # UstIshim
                initial_size=N_Ust,
                growth_rate=0,
                metadata={"name": "UstIshim", "sampling_time": t_UstIshim},
            ),
            msprime.PopulationConfiguration(  # Neanderthal
                initial_size=N_Nean,
                growth_rate=0,
                metadata={"name": "Neanderthal", "sampling_time": t_Altai},
            ),
            msprime.PopulationConfiguration(  # Basal Eurasian
                initial_size=N_Basal,
                growth_rate=0,
                metadata={"name": "BasalEurasian", "sampling_time": None},
            ),
        ],
        # Using columns in figure in Kamm paper as proxies for pop number
        demographic_events=[
            msprime.MassMigration(
                time=t_GhostWHG_to_Sard,
                source=2,
                destination=3,
                proportion=p_GhostWHG_to_Sard,
            ),
            msprime.MassMigration(
                time=t_Sard_LBK, source=2, destination=1, proportion=1.0
            ),
            msprime.PopulationParametersChange(
                time=t_Sard_LBK, initial_size=N_Sard_LBK, population_id=1
            ),
            msprime.MassMigration(
                time=t_Basal_to_EEF, source=1, destination=8, proportion=p_Basal_to_EEF
            ),
            msprime.MassMigration(
                time=t_LBK_Losch, source=1, destination=3, proportion=1.0
            ),
            msprime.MassMigration(
                time=t_MA1_Losch, source=4, destination=3, proportion=1.0
            ),
            msprime.PopulationParametersChange(
                time=t_Altai, initial_size=N_Nean, growth_rate=r_Nean, population_id=7
            ),
            msprime.MassMigration(
                time=t_Han_Losch, source=5, destination=3, proportion=1.0
            ),
            msprime.PopulationParametersChange(
                time=t_Han_Losch, initial_size=N_Han_Losch, population_id=3
            ),
            msprime.MassMigration(
                time=t_Ust_Losch, source=6, destination=3, proportion=1.0
            ),
            msprime.MassMigration(
                time=t_Nean_to_Eur, source=3, destination=7, proportion=p_Nean_to_Eur
            ),
            msprime.MassMigration(
                time=t_Basal_Losch, source=8, destination=3, proportion=1.0
            ),
            msprime.MassMigration(
                time=t_Mbu_Losch, source=0, destination=3, proportion=1.0
            ),
            msprime.PopulationParametersChange(
                time=t_Mbu_Losch, initial_size=N_Mbu_Losch, population_id=3
            ),
            msprime.PopulationParametersChange(
                time=t_Mbu_Losch,
                initial_size=N_Nean_Losch,
                growth_rate=0,
                population_id=7,
            ),
            msprime.MassMigration(
                time=t_Nean_Losch, source=7, destination=3, proportion=1.0
            ),
            msprime.PopulationParametersChange(
                time=t_Nean_Losch, initial_size=N_Nean_Losch, population_id=3
            ),
        ],
        # Mapping of per-epoch population names in the production model to
        # population IDs in this model.
        # Note: we'll be updating this map as the model is modernised.
        population_id_map=[pop_id_map] * 13,
        # this is a prediction of the paper rather than an assumption
        mutation_rate=1.22e-8,
    )


_species.get_demographic_model("AncientEurasia_9K19").register_qc(KammAncientSamples())


def DenisovanAncestryInPapuans():
    """
    Demographic model from Jacobs et al (2019). The model is
    illustrated on Figure S5, parameters are in Table S5
    """
    id = "QC-PapuansOutOfAfrica_10J19"
    populations = (
        [population_sample_0] * 4
        + [
            stdpopsim.Population("", "", 2058),
            stdpopsim.Population("", "", 2612),
        ]
        + [population_sample_none] * 4
    )

    generation_time = 29  # Just information

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
    t_Den2_to_Papuan = 45.7e3 / generation_time
    p_Den2_to_Papuan = (1 - p1) * 0.04
    t_Den1_to_Papuan = 29.8e3 / generation_time
    p_Den1_to_Papuan = p1 * 0.04
    t_Nean1_to_Papuan = 1412
    p_Nean1_to_Papuan = 0.002
    t_Nean1_to_CEU_CHB = 1566
    p_Nean1_to_CEU_CHB = 0.011
    t_Nean1_to_GhostA = 1853
    p_Nean1_to_GhostA = 0.024

    # set up populations
    population_configurations = [
        msprime.PopulationConfiguration(  # 0 YRI
            initial_size=N_YRI,
            growth_rate=0,
            metadata={"name": "YRI", "sampling_time": 0},
        ),
        msprime.PopulationConfiguration(  # 1 CEU
            initial_size=N_CEU,
            growth_rate=0,
            metadata={"name": "CEU", "sampling_time": 0},
        ),
        msprime.PopulationConfiguration(  # 2 CHB
            initial_size=N_CHB,
            growth_rate=0,
            metadata={"name": "CHB", "sampling_time": 0},
        ),
        msprime.PopulationConfiguration(  # 3 Papuan
            initial_size=N_Papuan,
            growth_rate=0,
            metadata={"name": "Papuan", "sampling_time": 0},
        ),
        msprime.PopulationConfiguration(  # 4 DenA
            initial_size=N_DenA,
            growth_rate=0,
            metadata={"name": "DenA", "sampling_time": t_DenA},
        ),
        msprime.PopulationConfiguration(  # 5 NeanA
            initial_size=N_NeanA,
            growth_rate=0,
            metadata={"name": "NeaA", "sampling_time": t_NeanA},
        ),
        msprime.PopulationConfiguration(  # 6 Den1
            initial_size=N_Den1,
            growth_rate=0,
            metadata={"name": "Den1", "sampling_time": None},
        ),
        msprime.PopulationConfiguration(  # 7 Den2
            initial_size=N_Den2,
            growth_rate=0,
            metadata={"name": "Den2", "sampling_time": None},
        ),
        msprime.PopulationConfiguration(  # 8 Nean1
            initial_size=N_Nean1,
            growth_rate=0,
            metadata={"name": "Nea1", "sampling_time": None},
        ),
        msprime.PopulationConfiguration(  # 9 Ghost
            initial_size=N_Ghost,
            growth_rate=0,
            metadata={"name": "Ghost", "sampling_time": None},
        ),
    ]

    migration_matrix = [[0] * 10 for _ in range(10)]
    migration_matrix[0][9] = m_YRI_Ghost
    migration_matrix[9][0] = m_YRI_Ghost
    migration_matrix[1][9] = m_Ghost_CEU
    migration_matrix[9][1] = m_Ghost_CEU
    migration_matrix[1][2] = m_CEU_CHB
    migration_matrix[2][1] = m_CEU_CHB
    migration_matrix[2][3] = m_CHB_Papuan
    migration_matrix[3][2] = m_CHB_Papuan

    demographic_events = [
        # Coalescence of CEU and CHB into CHB
        msprime.MassMigration(time=t_CEU_CHB, source=1, destination=2, proportion=1.0),
        # Set size of CEU+CHB population
        msprime.PopulationParametersChange(
            time=t_CEU_CHB, initial_size=N_CEU_CHB, population_id=2
        ),
        # Change migration matrix
        msprime.MigrationRateChange(time=t_CEU_CHB, rate=0, matrix_index=(2, 1)),
        msprime.MigrationRateChange(time=t_CEU_CHB, rate=0, matrix_index=(1, 2)),
        msprime.MigrationRateChange(time=t_CEU_CHB, rate=0, matrix_index=(3, 2)),
        msprime.MigrationRateChange(time=t_CEU_CHB, rate=0, matrix_index=(2, 3)),
        msprime.MigrationRateChange(time=t_CEU_CHB, rate=0, matrix_index=(9, 1)),
        msprime.MigrationRateChange(time=t_CEU_CHB, rate=0, matrix_index=(1, 9)),
        msprime.MigrationRateChange(
            time=t_CEU_CHB, rate=m_CEUCHB_Papua, matrix_index=(2, 3)
        ),
        msprime.MigrationRateChange(
            time=t_CEU_CHB, rate=m_CEUCHB_Papua, matrix_index=(3, 2)
        ),
        msprime.MigrationRateChange(
            time=t_CEU_CHB, rate=m_Ghost_CEUCHB, matrix_index=(2, 9)
        ),
        msprime.MigrationRateChange(
            time=t_CEU_CHB, rate=m_Ghost_CEUCHB, matrix_index=(9, 2)
        ),
        # Set bottleneck size of CEU+CHB population
        msprime.PopulationParametersChange(
            time=t_CEU_CHB_bot, initial_size=N_CEU_CHB_bot, population_id=2
        ),
        # Change migration matrix
        msprime.MigrationRateChange(time=t_CEU_CHB_bot, rate=0),
        # Set bottleneck size of Papuan population
        msprime.PopulationParametersChange(
            time=t_Papuan_bot, initial_size=N_Papuan_bot, population_id=3
        ),
        # Coalescence of CEU+CHB and Ghost to Ghost
        msprime.MassMigration(
            time=t_CEU_Ghost, source=2, destination=9, proportion=1.0
        ),
        # Coalescence of Papuan and Ghost to GhostA
        msprime.MassMigration(
            time=t_Papuan_Ghost, source=3, destination=9, proportion=1.0
        ),
        # Set bottleneck size of GhostA population
        msprime.PopulationParametersChange(
            time=t_GhostA_bot, initial_size=N_GhostA_bot, population_id=9
        ),
        # Coalescence of Ghost and YRI into Human (YRI)
        msprime.MassMigration(
            time=t_YRI_GhostA, source=9, destination=0, proportion=1.0
        ),
        # Set size of Human population
        msprime.PopulationParametersChange(
            time=t_YRI_GhostA, initial_size=N_Human, population_id=0
        ),
        # Coalescence of NeanA and Nean1 into NeanAnc
        msprime.MassMigration(
            time=t_Nean1_NeanA, source=8, destination=5, proportion=1.0
        ),
        # Set size of NeanAnc population
        msprime.PopulationParametersChange(
            time=t_Nean1_NeanA, initial_size=N_Nean1, population_id=5
        ),
        # Coalescence of Den1 and DenA into DenAnc (DenA)
        msprime.MassMigration(
            time=t_Den1_DenA, source=6, destination=4, proportion=1.0
        ),
        # Set size of DenA population
        msprime.PopulationParametersChange(
            time=t_Den1_DenA, initial_size=N_DenAnc, population_id=4
        ),
        # Coalescence of DenAnc and Den2 into DenAnc
        msprime.MassMigration(
            time=t_Den1_Den2, source=7, destination=4, proportion=1.0
        ),
        # Set size of DenA population
        msprime.PopulationParametersChange(
            time=t_Den1_Den2, initial_size=N_DenAnc, population_id=4
        ),
        # Coalescence of DenAnc and NeanAnc into Den_Nean (Nean1)
        msprime.MassMigration(time=t_Den_Nean, source=5, destination=4, proportion=1.0),
        # Set size of Den_Nean population
        msprime.PopulationParametersChange(
            time=t_Den_Nean, initial_size=N_Nean1, population_id=4
        ),
        # Coalescence of Den_Nean and Human into Anc (YRI)
        msprime.MassMigration(
            time=t_Human_Den_Nean, source=4, destination=0, proportion=1.0
        ),
        # Set ancestral size of population
        msprime.PopulationParametersChange(
            time=t_Human_Den_Nean, initial_size=N_A, population_id=0
        ),
        # Admixture events
        # Admixture from Den1 to Papuans
        msprime.MassMigration(
            time=t_Den1_to_Papuan, source=3, destination=6, proportion=p_Den1_to_Papuan
        ),
        # Admixture from Den2 to Papuans
        msprime.MassMigration(
            time=t_Den2_to_Papuan, source=3, destination=7, proportion=p_Den2_to_Papuan
        ),
        # Admixture from Nean1 to GhostA
        msprime.MassMigration(
            time=t_Nean1_to_GhostA,
            source=9,
            destination=8,
            proportion=p_Nean1_to_GhostA,
        ),
        # Admixture from Nean1 to CEU+CHB
        msprime.MassMigration(
            time=t_Nean1_to_CEU_CHB,
            source=2,
            destination=8,
            proportion=p_Nean1_to_CEU_CHB,
        ),
        # Admixture from Nean1 to Papuans
        msprime.MassMigration(
            time=t_Nean1_to_Papuan,
            source=3,
            destination=8,
            proportion=p_Nean1_to_Papuan,
        ),
        # Admixture from Neandertal to East Asia population
        msprime.MassMigration(
            time=t_Nean1_to_CHB, proportion=p_Nean1_to_CHB, source=2, destination=8
        ),
    ]
    demographic_events.sort(key=lambda x: x.time)

    pop_id_map = {
        "YRI": 0,
        "CEU": 1,
        "CHB": 2,
        "Papuan": 3,
        "DenA": 4,
        "NeaA": 5,
        "Den1": 6,
        "Den2": 7,
        "Nea1": 8,
        "Ghost": 9,
    }
    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=generation_time,
        populations=populations,
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events,
        # Mapping of per-epoch population names in the production model to
        # population IDs in this model.
        # Note: we'll be updating this map as the model is modernised.
        population_id_map=[pop_id_map] * 19,
        mutation_rate=1.4e-8,
    )


_species.get_demographic_model("PapuansOutOfAfrica_10J19").register_qc(
    DenisovanAncestryInPapuans()
)


def GutenkunstOOA():
    """
    From Gutenkunst et al (2009). Parameters taken from maximum likelihood values
    in Table 1 in that paper.
    """
    id = "QC-OutOfAfrica_3G09"
    populations = [population_sample_0] * 3

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
    T_AF = 220e3 / generation_time
    T_B = 140e3 / generation_time
    T_EU_AS = 21.2e3 / generation_time

    # Calculate population sizes at modern (T=0) time
    N_EUF = N_EU0 * math.exp(r_EU * T_EU_AS)
    N_ASF = N_AS0 * math.exp(r_AS * T_EU_AS)

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=generation_time,
        populations=populations,
        population_configurations=[
            msprime.PopulationConfiguration(initial_size=N_AF, growth_rate=0),
            msprime.PopulationConfiguration(initial_size=N_EUF, growth_rate=r_EU),
            msprime.PopulationConfiguration(initial_size=N_ASF, growth_rate=r_AS),
        ],
        # Setup initial migration matrix
        migration_matrix=[
            [0, m_AF_EU, m_AF_AS],
            [m_AF_EU, 0, m_EU_AS],
            [m_AF_AS, m_EU_AS, 0],
        ],
        demographic_events=[
            # CEU and CHB merge into B, reset migration rates to Af-B, change pop size
            msprime.MassMigration(time=T_EU_AS, source=2, dest=1, proportion=1),
            msprime.MigrationRateChange(time=T_EU_AS, rate=0),
            msprime.MigrationRateChange(time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
            msprime.MigrationRateChange(time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
            msprime.PopulationParametersChange(
                time=T_EU_AS, initial_size=N_B, population_id=1, growth_rate=0
            ),
            # B and AF merge, turn migration off, reset population size
            msprime.MassMigration(time=T_B, source=1, dest=0, proportion=1),
            msprime.MigrationRateChange(time=T_B, rate=0),
            # Ancestral size change, reset population size
            msprime.PopulationParametersChange(
                time=T_AF, initial_size=N_A, population_id=0
            ),
        ],
        # Mapping of per-epoch population names in the production model to
        # population IDs in this model.
        population_id_map=[
            {"YRI": 0, "CEU": 1, "CHB": 2},
            {"YRI": 0, "CEU": 1, "CHB": 2},
            {"YRI": 0, "CEU": 1, "CHB": 2},
            {"YRI": 0, "CEU": 1, "CHB": 2},
        ],
        mutation_rate=2.35e-8,
    )


_species.get_demographic_model("OutOfAfrica_3G09").register_qc(GutenkunstOOA())


def GladsteinAshkSubstructure():
    """
    Parameters largely taken from supplemental table 3 of Gladstein et al (2019).
    """
    id = "QC-AshkSub_7G19"
    # Population indices: YRI, CHB, CEU, M, J, WA, EA = 0, 1, 2, 3, 4, 5, 6
    # M = Middle eastern pop
    # J = Non-Ashk jewish pop
    populations = [population_sample_0] * 7

    generation_time = 25  # parameter estimate section and same as gutenkunst

    # Population sizes from supp tab 3
    N_ANC = 7300
    N_YRI = 10**4.26
    N_CEU = 10**4.52
    N_CHB = 10**3.61
    N_WA = 10**3.82
    N_EA = 10**6.29
    N_Ag = 10**3.04
    N_J = 10**5.55
    N_M = 10**5.64

    # Migration rate from CEU to ASHK ancestral pop
    m = 0.17

    # Epoch times in generations
    T_growth = 220e3 / generation_time  # anc. afr. growth time
    T_AF = 2105  # OOA event
    T_eu_as = 850  # CEU/CHB split
    T_EM = 481  # ME/J anc branch off CEU
    T_MJ = 211  # J/ME split
    T_A = 29  # J branch to anc. EA/WA pop
    T_m = T_A - 1  # from github issue 511 (time of gene flow)
    T_A_EW = 14  # WA/EA split
    T_Ag = T_A_EW - 1  # from github issue 511 (time of growth in WA/EA)

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=generation_time,
        populations=populations,
        population_configurations=[
            msprime.PopulationConfiguration(initial_size=N_YRI, growth_rate=0),
            msprime.PopulationConfiguration(initial_size=N_CHB, growth_rate=0),
            msprime.PopulationConfiguration(initial_size=N_CEU, growth_rate=0),
            msprime.PopulationConfiguration(initial_size=N_M, growth_rate=0),
            msprime.PopulationConfiguration(initial_size=N_J, growth_rate=0),
            msprime.PopulationConfiguration(initial_size=N_WA, growth_rate=0),
            msprime.PopulationConfiguration(initial_size=N_EA, growth_rate=0),
        ],
        # Population indices: YRI, CHB, CEU, M, J, WA, EA = 0, 1, 2, 3, 4, 5, 6
        demographic_events=[
            # EA & WA pop size change
            msprime.PopulationParametersChange(
                time=T_Ag, initial_size=N_Ag, population_id=5
            ),
            msprime.PopulationParametersChange(
                time=T_Ag, initial_size=N_Ag, population_id=6
            ),
            # EA & WA merge
            msprime.MassMigration(time=T_A_EW, source=6, dest=5, proportion=1),
            # Mass migration from europe into EA/WA anc.
            msprime.MassMigration(time=T_m, source=5, dest=2, proportion=m),
            # EA/WA anc merge with J
            msprime.MassMigration(time=T_A, source=5, dest=4, proportion=1),
            # J/ME merge
            msprime.MassMigration(time=T_MJ, source=4, dest=3, proportion=1),
            # J/ME anc merge with CEU
            msprime.MassMigration(time=T_EM, source=3, dest=2, proportion=1),
            # CEU/CHB merge
            msprime.MassMigration(time=T_eu_as, source=2, dest=1, proportion=1),
            # OOA
            msprime.MassMigration(time=T_AF, source=1, dest=0, proportion=1),
            # Ancestral size change, reset population size
            msprime.PopulationParametersChange(
                time=T_growth, initial_size=N_ANC, population_id=0
            ),
        ],
        # Mapping of per-epoch population names in the production model to
        # population IDs in this model.
        population_id_map=[
            {"YRI": 0, "CHB": 1, "CEU": 2, "ME": 3, "J": 4, "WAJ": 5, "EAJ": 6},
        ]
        * 10,
        mutation_rate=2.5e-8,
    )


_species.get_demographic_model("AshkSub_7G19").register_qc(GladsteinAshkSubstructure())


def ZigZag():
    """
    Model from Schiffels et al (2014) used to test inference methods on a "zigzag"
    demography. Single population, repeated growth and decay of pop size.
    """
    id = "QC-Zigzag_1S14"
    populations = [population_sample_0]

    # from section 7 of the supplement
    generation_time = 30
    mu = 1.25e-8
    L = 10000000
    theta = 7156
    Ne = theta / 4 / mu / L  # ancestral size for rescaling times/rates
    N0 = 5 * Ne  # multiply pop sizes to match 'ms -eN 0 5'

    population_configurations = [
        msprime.PopulationConfiguration(
            initial_size=N0, growth_rate=0, metadata={"sampling_time": 0}
        )
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
                time=time * 4 * Ne, growth_rate=rate / 4 / Ne, population_id=0
            )
        )

    de.append(
        msprime.PopulationParametersChange(
            time=0.596236 * 4 * Ne,
            growth_rate=0,
            population_id=0,
            initial_size=0.1 * N0,
        )
    )

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=generation_time,
        populations=populations,
        population_configurations=population_configurations,
        demographic_events=de,
        population_id_map=[{"generic": 0}] * 7,
        # This is a test model, not inferred from actual data.
        mutation_rate=None,
    )


_species.get_demographic_model("Zigzag_1S14").register_qc(ZigZag())


def JouganousOOA2017():
    """
    Four population OOA model from Jouganous et al. 2017, fit to the
    SFS using moments. The populations are YRI, CEU, CHB, and JPT.
    The model is shown in Figure 3 and the parameters are taken from
    table 4.
    The model includes population splits, directional migration and
    exponential growth in the various populations. The model was fit
    using data from 20 samples from each population and assumes a
    mutation rate of 1.44e-8 per generation and a generation time
    of 29 years.
    """
    id = "QC-OutOfAfrica_4J17"

    N_A = 11293
    N_AF = 23721
    N_B = 2831
    N_EU0 = 2512
    r_EU = 0.0016
    N_AS0 = 1019
    r_AS = 0.0026
    N_JP0 = 4384
    r_JP = 0.0129
    m_AF_B = 16.8e-5
    m_AF_EU = 1.14e-5
    m_AF_AS = 0.56e-5
    m_EU_AS = 4.75e-5
    m_CH_JP = 3.3e-5
    T_AF_y = 357e3
    T_AF_g = int(T_AF_y / 29)
    T_B_y = 119e3
    T_B_g = int(T_B_y / 29)
    T_EU_AS_y = 46e3
    T_EU_AS_g = int(T_EU_AS_y / 29)
    T_CH_JP_y = 9e3
    T_CH_JP_g = int(T_CH_JP_y / 29)

    N_curr_CEU = N_EU0 * math.exp(T_EU_AS_g * r_EU)
    N_curr_CHB = N_AS0 * math.exp(T_EU_AS_g * r_AS)
    N_curr_JPT = N_JP0 * math.exp(T_CH_JP_g * r_JP)

    de = msprime.Demography()
    de.add_population(name="YRI", initial_size=N_AF, initially_active=True)
    de.add_population(
        name="CEU", initial_size=N_curr_CEU, growth_rate=r_EU, initially_active=True
    )
    de.add_population(
        name="CHB", initial_size=N_curr_CHB, growth_rate=r_AS, initially_active=True
    )
    de.add_population(
        name="JPT", initial_size=N_curr_JPT, growth_rate=r_JP, initially_active=True
    )

    de.set_symmetric_migration_rate(["YRI", "CEU"], m_AF_EU)
    de.set_symmetric_migration_rate(["YRI", "CHB"], m_AF_AS)
    de.set_symmetric_migration_rate(["CEU", "CHB"], m_EU_AS)
    de.set_symmetric_migration_rate(["CHB", "JPT"], m_CH_JP)

    de.add_mass_migration(time=T_CH_JP_g, source="JPT", dest="CHB", proportion=1)
    de.add_symmetric_migration_rate_change(
        time=T_CH_JP_g, populations=["CHB", "JPT"], rate=0
    )
    de.add_mass_migration(time=T_EU_AS_g, source="CHB", dest="CEU", proportion=1)
    de.add_symmetric_migration_rate_change(
        time=T_EU_AS_g, populations=["CEU", "CHB"], rate=0
    )
    de.add_symmetric_migration_rate_change(
        time=T_EU_AS_g, populations=["YRI", "CHB"], rate=0
    )

    de.add_population_parameters_change(
        time=T_EU_AS_g, initial_size=N_B, growth_rate=0, population="CEU"
    )
    de.add_symmetric_migration_rate_change(
        time=T_EU_AS_g, populations=["CEU", "YRI"], rate=m_AF_B
    )
    de.add_mass_migration(time=T_B_g, source="CEU", dest="YRI", proportion=1)
    de.add_symmetric_migration_rate_change(
        time=T_B_g, populations=["CEU", "YRI"], rate=0
    )
    de.add_population_parameters_change(time=T_AF_g, initial_size=N_A, population="YRI")

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=29,
        model=de,
        mutation_rate=1.44e-8,
        population_id_map=None,
    )


_species.get_demographic_model("OutOfAfrica_4J17").register_qc(JouganousOOA2017())


def Boyko2008():
    """
    African-American two-epoch instantaneous growth model from Boyko et al
    2008, fit to the synonymous SFS for the 11 of 15 African Americans showing
    the least European ancestry, using coalescent simulations with
    recombination with the maximum likelihood method of Williamson et al 2005;
    times were calibrated assuming 3e5 generations since human-chimp divergence
    and fitting the number of synonymous human-chimp differences.
    """
    id = "QC-Africa_1B08"
    Nanc = 7778
    Ncurr = 25636
    Texp = 6809  # generations ago

    de = msprime.Demography()
    de.add_population(
        initial_size=Ncurr,
        name="African_Americans",
        description="African-Americans from Boyko et al",
    )

    de.add_population_parameters_change(
        time=Texp,
        initial_size=Nanc,
    )

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=None,
        model=de,
        mutation_rate=1.8e-8,
    )


_species.get_demographic_model("Africa_1B08").register_qc(Boyko2008())


def Iasi2021():
    """
    Extended adxmixture pulse model from Neandertals into modern humans from
    Iasi et al 2021. The implemented model corresponds to the simple model described
    in Supplemental Figure 1A.
    """
    id = "QC-NeaAdmixPulse_3I21"

    generation_time = 29

    # Population sizes
    N_AFR = 10000
    N_OOA = 10000
    N_NEA = 10000

    # Split times
    T_OOA = 2550  # in generations
    T_NEA = 10000

    # Admixture pulse times
    T_NEA_ADMIX_start = 1834  # 53.2 kya
    T_NEA_ADMIX_end = 1034  # 30 kya

    # Migration rates
    M_NEA_admix = 0.03

    de = msprime.Demography()
    de.add_population(
        initial_size=N_AFR,
        name="YRI",
    )
    de.add_population(
        initial_size=N_OOA,
        name="CEU",
    )
    de.add_population(
        initial_size=N_NEA,
        name="NEA",
    )

    # Copying `extended_pulse()` method here to get the specific migration rates
    # Note: migration_start and migration_end don't actually specify the start and
    # end of migration, but rather are used to create the gamma distribution of
    # migration rates. Values smaller than the migration_cutoff are then set to zero
    def extended_pulse(
        split_time,
        migration_start,
        migration_stop,
        total_migration_rate,
        source,
        dest,
        migration_cutoff=1e-5,
    ):

        """
        This function creates a dataframe of migration rate changes to simulate
        an extended pulse of unidirectional gene flow from a dest to a source
        population (forward in time) in msprime.
        The extended pulse models the migration rate m(t) as a rescaled
        Gamma distribution with a total contribution of migrant alleles entering
        between the start (backwards in time) (migration_start) and end
        (migration_stop) time of gene flow.
        The total migration rate is defined by the total_migration_rate.
        The start and end are not hard start and endpoints of the gene flow.
        The split_time time gives the maximum range (0 to split time)
        for the extended pulse.
        The migration_cutoff gives the lower cutoff for the per generation
        migration rate. The function returns a dataframe of the migration rates.
        """

        """
        The shape and scale parameters are calculated from the input.
        """
        tm = (migration_stop - migration_start) / 2 + migration_start
        k = 1 / ((migration_stop - migration_start) / (4 * tm)) ** 2
        m = stdpopsim.utils.gamma_pdf(x=range(split_time + 1), a=k, loc=0, scale=tm / k)

        """
        Scaling the distribution by the total migration rate.
        """
        m[abs(m) < migration_cutoff] = 0
        m_scaled = m * total_migration_rate / sum(m)

        """
        Filling the table of migration rate for each generation.
        """
        D_extended_pulse = []
        for x in range(split_time + 1):

            """
            Writing gene flow events which are inside the set time boarders and
            over the set migration cutoff. They will be included in the m(t)
            distribution.
            """
            rate = m_scaled[x]
            if rate > 0:
                D_extended_pulse.append(
                    dict(time=x, rate=rate, source=source, dest=dest)
                )

        return D_extended_pulse

    # Get a dataframe specifying the rates in the migration pulse
    pulse_df = extended_pulse(
        split_time=T_OOA,
        migration_start=T_NEA_ADMIX_start,
        migration_stop=T_NEA_ADMIX_end,
        total_migration_rate=M_NEA_admix,
        source=1,
        dest=2,
        migration_cutoff=1e-5,
    )

    # Need to insert zero migration rates before and after the pulse
    start = {"time": pulse_df[0]["time"] - 1, "rate": 0, "source": 1, "dest": 2}
    end = {"time": pulse_df[-1]["time"] + 1, "rate": 0, "source": 1, "dest": 2}
    pulse_df = [start] + pulse_df + [end]

    # Don't love iterating over rows here, but the DataFrame is small so works fine
    for record in pulse_df:
        de.add_migration_rate_change(
            time=record["time"],
            rate=record["rate"],
            source=record["source"],
            dest=record["dest"],
        )

    # Merge CEU into YRI and NEA into YRI
    de.add_mass_migration(time=T_OOA, source="CEU", dest="YRI", proportion=1)
    de.add_mass_migration(time=T_NEA, source="NEA", dest="YRI", proportion=1)

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=generation_time,
        model=de,
        mutation_rate=2e-8,
    )


_species.get_demographic_model(
    "OutOfAfricaExtendedNeandertalAdmixturePulse_3I21"
).register_qc(Iasi2021())


# Currently this is not in use, but we kept it here in case we need it when the
# catalog is updated to use demes.
_ancient_europe_demes_str = """
time_units: generations

demes:
  - name: OOA
    epochs:
      - start_size: 5000
        end_time: 1500
  - name: NE
    ancestors: [OOA]
    epochs:
      - start_size: 5000
        end_time: 600
  - name: WA
    ancestors: [OOA]
    epochs:
      - start_size: 5000
        end_time: 800
  - name: CHG
    ancestors: [WA]
    epochs:
      - start_size: 10000
        end_time: 180
  - name: ANA
    ancestors: [WA]
    epochs:
      - start_size: 50000
        end_time: 200
  - name: WHG
    ancestors: [NE]
    epochs:
      - start_size: 10000
        end_time: 200
  - name: EHG
    ancestors: [NE]
    epochs:
      - start_size: 10000
        end_time: 180
  - name: YAM
    ancestors: [EHG, CHG]
    proportions: [0.5, 0.5]
    start_time: 180
    epochs:
      - start_size: 5000
        end_time: 140
  - name: NEO
    ancestors: [WHG, ANA]
    proportions: [0.25, 0.75]
    start_time: 200
    epochs:
      - start_size: 50000
        end_time: 140
  - name: Bronze
    ancestors: [YAM, NEO]
    proportions: [0.5, 0.5]
    start_time: 140
    epochs:
      - start_size: 50000
        end_size: 592450737.709282
        end_time: 0
"""


def PearsonAncientEurope():
    """
    Demographic history of ancient Europe from Allentoft et al. 2022. It
    describes population splits and migration pulses for four main branches
    in Fig. S3i.1 in Supplementary Information (Part I) found at
    https://www.biorxiv.org/content/10.1101/2022.05.04.490594v2.supplementary-material
    The branches are colored (black, red, yellow and purple).
    """
    id = "QC-AncientEurope_4A21"
    r_EU = 0.067  # TODO: not sure where this came from
    demog = msprime.Demography()
    # TODO: not clear from the picture that the OOA pop started with 5000 people
    demog.add_population(
        name="OOA", description="Out of Africa population", initial_size=5_000
    )
    demog.add_population(
        name="NE", description="Northern Europeans", initial_size=5_000
    )
    demog.add_population(name="WA", description="West Asians", initial_size=5_000)
    demog.add_population_split(time=1_500, derived=["NE", "WA"], ancestral="OOA")
    demog.add_population(
        name="CHG",
        description="Caucus Hunter-gatherers",
        initial_size=10_000,
        default_sampling_time=300,
    )
    demog.add_population(
        name="ANA",
        description="Anatolian Farmers",
        initial_size=50_000,
        default_sampling_time=260,
    )
    demog.add_population_split(time=800, derived=["ANA", "CHG"], ancestral="WA")
    demog.add_population(
        name="WHG",
        description="Western hunter-gatherers",
        initial_size=10_000,
        default_sampling_time=250,
    )
    demog.add_population(
        name="EHG",
        description="Eastern hunter-gatherers",
        initial_size=10_000,
        default_sampling_time=250,
    )
    demog.add_population_split(time=600, derived=["EHG", "WHG"], ancestral="NE")
    demog.add_population(
        name="YAM",
        description="Yamnayans",
        initial_size=5_000,
        default_sampling_time=160,
    )
    demog.add_population(
        name="NEO",
        description="Neolithic farmers",
        initial_size=50_000,
        default_sampling_time=180,
    )
    demog.add_admixture(
        time=200, ancestral=["WHG", "ANA"], proportions=[0.25, 0.75], derived="NEO"
    )
    demog.add_admixture(
        time=180, ancestral=["EHG", "CHG"], proportions=[0.5, 0.5], derived="YAM"
    )
    demog.add_population(
        name="Bronze",
        description="Bronze age",
        initial_size=592450737.7092822,
        growth_rate=r_EU,
        default_sampling_time=135,
    )
    demog.add_admixture(
        time=140, ancestral=["YAM", "NEO"], proportions=[0.5, 0.5], derived="Bronze"
    )
    demog.sort_events()
    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=29,  # TODO: not clear where this came from,
        # but it is not necessary, given the times are in generations already
        model=demog,
        mutation_rate=1.25e-8,
    )


_species.get_demographic_model("AncientEurope_4A21").register_qc(PearsonAncientEurope())


def Huber2017():
    """
    Gamma DFE from Huber et al. 2017 PNAS.
    """

    id = "Huber2017_gamma_dfe"
    neutral = stdpopsim.MutationType()
    negative = stdpopsim.MutationType(
        dominance_coeff=0.5,
        distribution_type="g",  # gamma distribution
        distribution_args=[-0.014, 0.19],
    )

    return stdpopsim.DFE(
        id=id,
        description=id,
        long_description=id,
        mutation_types=[neutral, negative],
        proportions=[0.3, 0.7],
    )


_species.get_dfe("Gamma_H17").register_qc(Huber2017())
