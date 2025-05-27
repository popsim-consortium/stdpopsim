import msprime
import stdpopsim
import numpy as np

_species = stdpopsim.get_species("PruDul")


def _pop1_almond():
    # the size during the interval times[k] to times[k+1] = sizes[k]
    times = np.array(
        [
            14440,
            22074,
            26655,
            29949,
            32601,
            34771,
            36700,
            38307,
            41120,
            43611,
            45861,
            47950,
            49959,
            51888,
        ]
    )
    sizes = np.array(
        [
            216523,
            219758,
            110605,
            65380,
            49918,
            45483,
            44353,
            47616,
            57496,
            69583,
            83878,
            98724,
            116879,
            133409,
        ]
    )

    demographic_events = []
    for sz, t in zip(sizes, times):
        demographic_events.append(
            msprime.PopulationParametersChange(time=t, initial_size=sz, population_id=0)
        )

    populations = [
        stdpopsim.Population(
            id="Prunus_dulcis",
            description="Prunus dulcis population",
        )
    ]

    return stdpopsim.DemographicModel(
        id="AlmondHist_1V16",
        description="Prunus_dulcis_population",
        long_description="""
            This model comes from MSMC using five
            individuals (PD03, PD04, PD05, PD06, and PD07) from Turkey,
              Iran, Pakistan and Italy.
            The model is estimated with 14 time periods.
            depiction of the model can be found in Figure S8a of the paper.
        """,
        populations=populations,
        citations=[
            stdpopsim.Citation(
                author="Velasco et al.",
                year=2016,
                doi="https://doi.org/10.1534/g3.116.032672",
                reasons={stdpopsim.CiteReason.DEM_MODEL},
            )
        ],
        generation_time=10,
        mutation_rate=10e-8,
        demographic_events=demographic_events,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=sizes[0], metadata=populations[0].asdict()
            )
        ],
    )


_species.add_demographic_model(_pop1_almond())
