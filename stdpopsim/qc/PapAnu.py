import msprime
import numpy as np

import stdpopsim


_species = stdpopsim.get_species("PapAnu")


def WallOnePopPapio():
    id = "QC-SinglePopSMCpp_1W22"
    populations = [stdpopsim.Population("PAnubis_SNPRC", "")]

    # This is a demographic model demographic model produced by SMC++
    # generation time of 11 years and a mutation rate of 5.7e-9 per bp per gen
    generation_time = 11
    mutation_rate = 5.7 * 10e-9

    # the time and population sizes are on Sup. Fig2
    # the exact values were provided by the authours
    # here years are divided by its generation time
    times = np.array(
        [
            0.0000,
            221.2812,
            463.8092,
            15027.7847,
            39328.5751,
            71092.5441,
            110830.9203,
            186053.2682,
        ]
    )
    sizes = np.array(
        [
            335505.4808,
            120758.1302,
            51822.58297,
            41841.54229,
            30714.33863,
            72998.86202,
            55968.42221,
            93362.02606,
        ]
    )

    demographic_events = []
    for sz, t in zip(sizes, times):
        demographic_events.append(
            msprime.PopulationParametersChange(time=t, initial_size=sz, population_id=0)
        )

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=generation_time,
        populations=populations,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=335505, metadata=populations[0].asdict()
            ),
        ],
        demographic_events=demographic_events,
        mutation_rate=mutation_rate,
    )


_species.get_demographic_model("SinglePopSMCpp_1W22").register_qc(WallOnePopPapio())
