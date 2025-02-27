import msprime
import math
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
        mutation_rate=1.450e-9,
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
        mutation_rate=8.4e-9,
    )


_species.get_demographic_model("African3Epoch_1S16").register_qc(SheehanSongThreeEpic())


def Gamma_H17():
    id = "Gamma_H17"
    description = "Deleterious Gamma DFE"
    long_description = "QC version"
    neutral = stdpopsim.MutationType()
    gamma_shape = 0.33
    gamma_scale = 1.2e-3
    gamma_mean = gamma_shape * gamma_scale
    h = 0.5  # dominance coefficient
    negative = stdpopsim.MutationType(
        dominance_coeff=h,
        distribution_type="g",  # gamma distribution
        # (1+s for homozygote in SLiM versus 1+2s in dadi)
        distribution_args=[-2 * gamma_mean, gamma_shape],
    )
    # LNS = 2.85 * LS
    # prop_synonymous = 1/(1+2.85) = 0.26
    prop_synonymous = 0.26
    return stdpopsim.DFE(
        id=id,
        description=description,
        long_description=long_description,
        mutation_types=[neutral, negative],
        proportions=[prop_synonymous, 1 - prop_synonymous],
    )


_species.get_dfe("Gamma_H17").register_qc(Gamma_H17())


def ZhenPos():
    id = "GammaPos_H17"
    description = "Deleterious Gamma DFE with fixed-s beneficials"
    # Same model as Huber 2017
    neutral = stdpopsim.MutationType()
    gamma_shape = 0.33
    gamma_scale = 6.01e-4
    gamma_mean = gamma_shape * gamma_scale
    h = 0.5  # dominance coefficient
    negative = stdpopsim.MutationType(
        dominance_coeff=h,
        distribution_type="g",  # gamma distribution
        distribution_args=[round(-2 * gamma_mean, 7), gamma_shape],
    )
    # LNS = 2.85 * LS
    # prop_synonymous = 1/(1+2.85) = 0.26
    prop_synonymous = 0.26
    prop_beneficial = (1 - prop_synonymous) * 6.75e-4
    selection_coefficient = 10 ** (-4.801)
    prop_deleterious = 1 - (prop_synonymous + prop_beneficial)
    positive = stdpopsim.MutationType(
        dominance_coeff=0.5,
        distribution_type="f",
        distribution_args=[selection_coefficient],
    )
    return stdpopsim.DFE(
        id=id,
        description=description,
        long_description=description,
        mutation_types=[neutral, negative, positive],
        proportions=[prop_synonymous, prop_deleterious, prop_beneficial],
    )


_species.get_dfe("GammaPos_H17").register_qc(ZhenPos())

def RagsdalePos():
    id = "RagsdalePos"
    description = "lognormal DFE based on Ragsdale et al. 2021"
    neutral = stdpopsim.MutationType()
    # selection coefficients are scaled by 2xeffective population size
    Ne = 2.8e6
    # NS = 2.5 * S
    prop_synonymous = 1 / (1 + 2.5)
    prop_nonsynonymous = 1 - prop_synonymous
    prop_beneficial = 0.0079 * prop_nonsynonymous
    prop_nonsynonymous = prop_nonsynonymous - prop_beneficial
    # beneficial selection coefficient
    sval_pos = 39.9 / Ne / 2
    # scale for DaDi
    sval_pos = 2 * sval_pos
    positive = stdpopsim.MutationType(
        dominance_coeff=0.5,
        distribution_type="f",
        distribution_args=[sval_pos],
    )
    # lognormal DFE parameters for deleterious mutations,
    # scale by 2Ne
    muval = 5.42 - math.log(2 * Ne)
    sigmaval = 3.36
    # adjust mu so that mean is 2x former mean
    # (to match DaDi's scaling)
    expectedmean = math.exp(muval + (sigmaval**2 / 2))
    targetmean = expectedmean * 2
    muval = math.log(targetmean) - (sigmaval**2 / 2)
    negative = stdpopsim.MutationType(
        dominance_coeff=0.5,
        distribution_type="ln",
        distribution_args=[muval, sigmaval],
    )
    return stdpopsim.DFE(
        id=id,
        description=description,
        long_description=description,
        mutation_types=[neutral, negative, positive],
        proportions=[prop_synonymous, prop_nonsynonymous, prop_beneficial],
    )


_species.get_dfe("LognormalPlusPositive_R16").register_qc(RagsdalePos())
