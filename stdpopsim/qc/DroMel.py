import msprime

import stdpopsim


_species = stdpopsim.get_species("DroMel")


def LiStephanTwoPopulation():
    id = "QC-OutOfAfrica_2L06"
    populations = [
        stdpopsim.Population("AFR", ""),
        stdpopsim.Population("EUR", ""),
    ]

    # Parameters for the African population are taken from the section Demographic
    # History of the African Population
    generation_time = 0.1  # 10  generations per year
    N_A0 = 8.603e6  # modern African pop. size
    N_A1 = N_A0 / 5.0  # African pop. size before expansion

    # Parameters for the European population are taken from the section Demographic
    # History of the European Population
    N_E0 = 1.075e6  # modern European pop. size
    N_E1 = 2.2e3  # European founder pop. size

    # Times from from the section Demographic History of the * Population
    T_A0 = 6e4 / generation_time  # time of 1st expansion in African pop.
    T_E_A = 15.8e3 / generation_time  # European/African divergence time
    T_EE = T_E_A - 340 / generation_time  # Time of European pop. re-expansion

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=generation_time,
        populations=populations,
        # Set population sizes at T=0
        # pop0 is Africa, pop1 is Europe
        population_configurations=[
            msprime.PopulationConfiguration(initial_size=N_A0, growth_rate=0),
            msprime.PopulationConfiguration(initial_size=N_E0, growth_rate=0),
        ],
        # Now we add the demographic events working backwards in time.
        demographic_events=[
            # OOA bottleneck
            msprime.PopulationParametersChange(
                time=T_EE, initial_size=N_E1, population_id=1
            ),
            # E and A coalesce
            msprime.MassMigration(time=T_E_A, source=1, destination=0, proportion=1.0),
            # Pre-expansion Africa
            msprime.PopulationParametersChange(
                time=T_A0, initial_size=N_A1, population_id=0
            ),
        ],
        population_id_map=[
            {"AFR": 0, "EUR": 1},
            {"AFR": 0, "EUR": 1},
            {"AFR": 0, "EUR": 1},
            {"AFR": 0, "EUR": 1},
        ],
    )


_species.get_demographic_model("OutOfAfrica_2L06").register_qc(LiStephanTwoPopulation())


def SheehanSongThreeEpic():
    id = "QC-African3Epoch_1S16"
    populations = [stdpopsim.Population("AFR", "")]

    # Model from paper https://doi.org/10.1371/journal.pcbi.1004845

    # Parameters are taken from table 7 using the average stat prediction values
    # as those were generally stated to be the best
    N_1 = 544.2e3  # recent
    N_2 = 145.3e3  # bottleneck
    N_3 = 652.7e3  # ancestral

    # Times taken from simulating data section based on PSMC and converted to
    # number of generations from coalescent units using the baseline effective
    # population size. Note that the coalescent values are calculated by
    N_ref = 1e5
    t_1_coal = 0.5
    t_2_coal = 5
    T_1 = t_1_coal * 4 * N_ref
    T_2 = (t_1_coal + t_2_coal) * 4 * N_ref

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=_species.generation_time,
        populations=populations,
        # Set population sizes at T=0
        population_configurations=[
            msprime.PopulationConfiguration(initial_size=N_1, growth_rate=0),
        ],
        # Now we add the demographic events working backwards in time.
        demographic_events=[
            # Bottleneck
            msprime.PopulationParametersChange(
                time=T_1, initial_size=N_2, population_id=0
            ),
            # Ancestral population size
            msprime.PopulationParametersChange(
                time=T_2, initial_size=N_3, population_id=0
            ),
        ],
        population_id_map=[{"AFR": 0}] * 3,
    )


_species.get_demographic_model("African3Epoch_1S16").register_qc(SheehanSongThreeEpic())
