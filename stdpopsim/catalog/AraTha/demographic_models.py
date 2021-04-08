import msprime
import numpy as np

import stdpopsim

_species = stdpopsim.get_species("AraTha")


def _sma_1pop():
    # the size during the interval times[k] to times[k+1] = sizes[k]
    times = np.array(
        [
            699,
            2796,
            6068,
            9894,
            14370,
            19606,
            25730,
            32894,
            41275,
            51077,
            62544,
            75958,
            91648,
            110001,
            131471,
            156584,
            185960,
            220324,
            260520,
            307540,
            362541,
            426879,
            502139,
            590173,
            693151,
            813610,
            954517,
            1119341,
            1312147,
            1537686,
            1801500,
            2110100,
        ]
    )
    sizes = np.array(
        [
            42252426,
            42252426,
            60323,
            72174,
            40591,
            21158,
            21442,
            39942,
            78908,
            111132,
            110745,
            96283,
            87661,
            83932,
            83829,
            91813,
            111644,
            143456,
            181571,
            217331,
            241400,
            246984,
            238593,
            228222,
            217752,
            198019,
            165210,
            121796,
            121796,
            73989,
            73989,
            73989,
        ]
    )

    # MSMC is accurate from 40Kya-1.6Mya for A.thaliana (Durvasula et al 2017)
    # set the first 7 sizes
    # equal to the size at 8 (~40Kya)
    sizes[:8] = sizes[8]
    # set the last 2 entries equal
    # to the size at 30 (~1.6Mya)
    sizes[30:32] = sizes[30]

    demographic_events = []
    for sz, t in zip(sizes, times):
        demographic_events.append(
            msprime.PopulationParametersChange(time=t, initial_size=sz, population_id=0)
        )

    populations = [
        stdpopsim.Population(
            id="SouthMiddleAtlas",
            description="Arabidopsis Thaliana South Middle Atlas population",
        )
    ]

    return stdpopsim.DemographicModel(
        id="SouthMiddleAtlas_1D17",
        description="South Middle Atlas piecewise constant size",
        long_description="""
            This model comes from MSMC using two randomly sampled homozygous
            individuals (Khe32 and Ifr4) from the South Middle Atlas region
            from the Middle Atlas Mountains in Morocco. The model is estimated
            with 32 time periods. Because estimates from the recent and ancient
            past are less accurate, we set the population size in the first 7
            time periods equal to the size at the 8th time period and the size
            during last 2 time periods equal to the size in the 30th time
            period.
        """,
        populations=populations,
        citations=[
            stdpopsim.Citation(
                author="Durvasula et al.",
                year=2017,
                doi="https://doi.org/10.1073/pnas.1616736114",
                reasons={stdpopsim.CiteReason.DEM_MODEL},
            )
        ],
        generation_time=1,
        demographic_events=demographic_events,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=sizes[0], metadata=populations[0].asdict()
            )
        ],
    )


_species.add_demographic_model(_sma_1pop())


def _afr_2epoch():
    N_A = 746148
    N_0 = 100218
    t_1 = 568344
    populations = [
        stdpopsim.Population(
            id="SouthMiddleAtlas",
            description="Arabidopsis Thaliana South Middle Atlas population",
        )
    ]
    return stdpopsim.DemographicModel(
        id="African2Epoch_1H18",
        description="South Middle Atlas African two epoch model",
        long_description="""
            Model estimated from site frequency spectrum of synonymous
            SNPs from African South Middle Atlas samples using
            Williamson et al. 2005 methodology. Values come from supplementary
            table 1 of Huber et al 2018. Sizes change from N_A -> N_0 and t_1 is
            time of the second epoch.
        """,
        populations=populations,
        citations=[
            stdpopsim.Citation(
                author="Huber et al.",
                year=2018,
                doi="https://doi.org/10.1038/s41467-018-05281-7",
                reasons={stdpopsim.CiteReason.DEM_MODEL},
            )
        ],
        generation_time=1,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_0, metadata=populations[0].asdict()
            )
        ],
        demographic_events=[
            msprime.PopulationParametersChange(
                time=t_1, initial_size=N_A, population_id=0
            )
        ],
    )


_species.add_demographic_model(_afr_2epoch())


def _afr_3epoch():
    # values from the supplementary Table 1.
    # the size changed as N_A -> N_2 -> N_3.
    # t_2 is time of 2nd epoch and t_3 of the third epoch
    N_A = 161744
    N_2 = 24076
    N_3 = 203077
    t_2 = 7420
    t_3 = 14534
    populations = [
        stdpopsim.Population(
            id="SouthMiddleAtlas",
            description="Arabidopsis Thaliana South Middle Atlas population",
        )
    ]
    return stdpopsim.DemographicModel(
        id="African3Epoch_1H18",
        description="South Middle Atlas African three epoch model",
        long_description="""
            Model estimated from site frequency spectrum of synonymous
            SNPs from African (South Middle Atlas) samples using Williamson et
            al. 2005 methodology. Values come from supplementary table 1 of
            Huber et al 2018. Sizes change from N_A -> N_2 -> N_3 and t_2 is
            the time of the second epoch and t_3 is the time of the 3rd epoch.
        """,
        populations=populations,
        citations=[
            stdpopsim.Citation(
                author="Huber et al.",
                year=2018,
                doi="https://doi.org/10.1038/s41467-018-05281-7",
                reasons={stdpopsim.CiteReason.DEM_MODEL},
            )
        ],
        generation_time=1,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_3, metadata=populations[0].asdict()
            )
        ],
        demographic_events=[
            msprime.PopulationParametersChange(
                time=t_3, initial_size=N_2, population_id=0
            ),
            msprime.PopulationParametersChange(
                time=t_2 + t_3, initial_size=N_A, population_id=0
            ),
        ],
    )


_species.add_demographic_model(_afr_3epoch())
