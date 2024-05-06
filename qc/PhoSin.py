import msprime

import stdpopsim


_species = stdpopsim.get_species("PhoSin")


def Robinson2022_TwoEpoch():
    id = "QC-Vaquita2Epoch_1R22"

    populations = [
        stdpopsim.Population(id="Vaquita", description=""),
    ]

    # Time of second epoch
    T_2 = 2162
    # population sizes
    N_ANC = 4485
    N_2 = 2807

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=_species.generation_time,
        mutation_rate=5.83e-9,
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
        population_id_map=[{"Vaquita": 0}] * 2,
    )


_species.get_demographic_model("Vaquita2Epoch_1R22").register_qc(
    Robinson2022_TwoEpoch()
)


def Robinson2022():
    """
    Gamma DFE from Robinson et al. 2022 Science.
    """

    id = "Robinson2022_gamma_dfe"
    neutral = stdpopsim.MutationType()

    negative = stdpopsim.MutationType(
        dominance_coeff_list=[0.0, 0.01, 0.1, 0.4],
        dominance_coeff_breaks=[-0.1, -0.01, -0.001],
        distribution_type="g",  # gamma distribution
        distribution_args=[-0.0257, 0.131],
    )

    return stdpopsim.DFE(
        id=id,
        description=id,
        long_description=id,
        mutation_types=[neutral, negative],
        proportions=[0.3, 0.7],
    )


_species.get_dfe("Gamma_R22").register_qc(Robinson2022())
