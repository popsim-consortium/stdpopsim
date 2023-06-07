import msprime
import numpy as np

import stdpopsim


_species = stdpopsim.get_species("AraTha")


def Durvasula2017MSMC():
    id = "QC-SouthMiddleAtlas_1D17"
    populations = [stdpopsim.Population("SouthMiddleAtlas", "")]

    # Both of the following are directly
    # converted from MSMC output scaled by A.Thaliana
    # mutation rate 7e-9 and 1 generation
    # per year.

    times = np.array(
        [
            6.990000e02,
            2.796000e03,
            6.068000e03,
            9.894000e03,
            1.437000e04,
            1.960600e04,
            2.573000e04,
            3.289400e04,
            4.127500e04,
            5.107700e04,
            6.254400e04,
            7.595800e04,
            9.164800e04,
            1.100010e05,
            1.314710e05,
            1.565840e05,
            1.859600e05,
            2.203240e05,
            2.605200e05,
            3.075400e05,
            3.625410e05,
            4.268790e05,
            5.021390e05,
            5.901730e05,
            6.931510e05,
            8.136100e05,
            9.545170e05,
            1.119341e06,
            1.312147e06,
            1.537686e06,
            1.801500e06,
            2.110100e06,
        ]
    )

    sizes = np.array(
        [
            4.2252426e07,
            4.2252426e07,
            6.0323000e04,
            7.2174000e04,
            4.0591000e04,
            2.1158000e04,
            2.1442000e04,
            3.9942000e04,
            7.8908000e04,
            1.1113200e05,
            1.1074500e05,
            9.6283000e04,
            8.7661000e04,
            8.3932000e04,
            8.3829000e04,
            9.1813000e04,
            1.1164400e05,
            1.4345600e05,
            1.8157100e05,
            2.1733100e05,
            2.4140000e05,
            2.4698400e05,
            2.3859300e05,
            2.2822200e05,
            2.1775200e05,
            1.9801900e05,
            1.6521000e05,
            1.2179600e05,
            1.2179600e05,
            7.3989000e04,
            7.3989000e04,
            7.3989000e04,
        ]
    )

    # The first 8 epochs are "masked" to
    # the last Ne at 40kya due to
    # the limitations of MSMC to infer
    # population size in this range.
    #
    # Similarly, the last 2 entries
    # are set to equal the third last.
    #
    # Durvasula et al 2017 shows that
    # MSMC has power in A.Thaliana
    # between 40kya and 1.6Mya.
    sizes[:8] = sizes[8]
    sizes[30:32] = sizes[30]
    demographic_events = []
    population_configurations = [msprime.PopulationConfiguration(initial_size=sizes[0])]

    for i, t in enumerate(times):
        demographic_events.append(
            msprime.PopulationParametersChange(
                time=t, initial_size=sizes[i], population_id=0
            )
        )

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=_species.generation_time,
        populations=populations,
        population_configurations=population_configurations,
        demographic_events=demographic_events,
        population_id_map=[{"SouthMiddleAtlas": 0}] * 33,
        mutation_rate=7.1e-9,
    )


_species.get_demographic_model("SouthMiddleAtlas_1D17").register_qc(Durvasula2017MSMC())


def HuberTwoEpoch():
    id = "QC-African2Epoch_1H18"
    populations = [
        stdpopsim.Population(id="SouthMiddleAtlas", description="A. thalina"),
    ]

    # Time of second epoch
    T_2 = 568344
    # population sizes
    N_ANC = 746148
    N_2 = 100218

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=_species.generation_time,
        populations=populations,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_2, metadata=populations[0].asdict()
            ),
        ],
        demographic_events=[
            msprime.PopulationParametersChange(
                time=T_2, initial_size=N_ANC, population_id=0
            ),
        ],
        population_id_map=[{"SouthMiddleAtlas": 0}] * 2,
        # Huber et al say "7e-9" but then refer to Ossowski,
        # which Durvasula reported as giving 7.1e-9
        mutation_rate=7e-9,
    )


_species.get_demographic_model("African2Epoch_1H18").register_qc(HuberTwoEpoch())


def HuberThreeEpoch():
    id = "QC-African3Epoch_1H18"
    populations = [
        stdpopsim.Population(id="SouthMiddleAtlas", description="A. thalina"),
    ]

    # Time of second epoch
    T_2 = 7420
    T_3 = 14534
    # population sizes
    N_ANC = 161744
    N_2 = 24076
    N_3 = 203077

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=_species.generation_time,
        populations=populations,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_3, metadata=populations[0].asdict()
            ),
        ],
        demographic_events=[
            msprime.PopulationParametersChange(
                time=T_3, initial_size=N_2, population_id=0
            ),
            msprime.PopulationParametersChange(
                time=T_2 + T_3, initial_size=N_ANC, population_id=0
            ),
        ],
        population_id_map=[{"SouthMiddleAtlas": 0}] * 3,
        # Huber et al say "7e-9" but then refer to Ossowski,
        # which Durvasula reported as giving 7.1e-9
        mutation_rate=7e-9,
    )


_species.get_demographic_model("African3Epoch_1H18").register_qc(HuberThreeEpoch())


def Huber2018DFE():
    # From Supp Table 4, line "Genome-Wide; Only A.thaliana; Additive"
    id = "Gamma_H18"
    gamma_shape = 0.155
    gamma_scale = 0.00612
    neutral = stdpopsim.MutationType()
    negative = stdpopsim.MutationType(
        dominance_coeff=0.5,
        distribution_type="g",  # gamma distribution
        # extra factor of 2 is to convert dadi to SLiM
        # (1+s for homozygote in SLiM versus 1+2s in dadi)
        distribution_args=[-2 * gamma_shape * gamma_scale, gamma_shape],  # mean, shape
    )
    # proportion of single-nuc changes that are synonymous,
    # from counting the codon table
    prop_neutral = 330 / 1728

    return stdpopsim.DFE(
        id=id,
        description=id,
        long_description=id,
        mutation_types=[neutral, negative],
        proportions=[prop_neutral, 1 - prop_neutral],
    )


_species.get_dfe("Gamma_H18").register_qc(Huber2018DFE())
