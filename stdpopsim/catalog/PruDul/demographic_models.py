import msprime
import stdpopsim
import numpy as np

_species = stdpopsim.get_species("PruDul")


def _pop1_almond():
    # the size during the interval times[k] to times[k+1] = sizes[k]
    times = np.array(
    [
        1821,
        3726,
        5685,
        7754,
        9879,
        12047,
        14409,
        16819,
        21943,
        27673,
        34063,
        41325,
        49651,
        59366
    ]
    )
    sizes = np.array(
            [
        216627,
        220902,
        111163,
        66270,
        49881,
        44180,
        43467,
        49168,
        58432,
        70546,
        84798,
        100475,
        117577,
        134679
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
        mutation_rate=1e-8,
        demographic_events=demographic_events,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=sizes[0], metadata=populations[0].asdict()
            )
        ],
    )


_species.add_demographic_model(_pop1_almond())
