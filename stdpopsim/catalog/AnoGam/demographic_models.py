import msprime
import numpy as np
import stdpopsim

_species = stdpopsim.get_species("AnoGam")


def _GAS_sp():
    # the size during the interval times[k] to times[k+1] = sizes[k]
    times = np.array(
        [
            1.08221472e01,
            1.77815418e01,
            2.43957877e01,
            2.62194838e01,
            2.80715527e01,
            2.99059882e01,
            3.17590591e01,
            3.36213534e01,
            3.74218265e01,
            6.87139222e01,
            7.13468519e01,
            7.40546632e01,
            7.69788992e01,
            8.00932536e01,
            8.33883683e01,
            9.33624771e01,
            1.05867180e02,
            1.19778311e02,
            1.34317538e02,
            1.49486416e02,
            1.65081438e02,
            1.81069660e02,
            1.96884420e02,
            2.13010501e02,
            2.28659684e02,
            2.43883726e02,
            2.58403062e02,
            2.70573820e02,
            2.75664382e02,
            2.80730912e02,
            2.85888999e02,
            2.91187697e02,
            2.96632655e02,
            3.02218199e02,
            3.07951917e02,
            3.13860240e02,
            3.19955507e02,
            3.26244790e02,
            3.32745865e02,
            4.77774801e02,
            5.64379009e02,
            6.19045020e02,
            6.39454764e02,
            6.61129552e02,
            6.84129309e02,
            7.08609246e02,
            7.34719738e02,
            7.62657505e02,
            7.92661043e02,
            8.24963147e02,
            8.59801339e02,
            9.38516578e02,
            1.14545402e03,
            1.21187638e03,
            1.42580467e03,
            1.66762718e03,
            1.94206662e03,
            2.61982155e03,
            4.14416552e03,
            4.86878040e03,
            5.74277111e03,
            6.06376420e03,
            6.48115044e03,
            7.16974182e03,
            9.94668728e03,
            1.08341379e04,
            1.12602248e04,
        ]
    )
    sizes = np.array(
        [
            2570469.78880419,
            2398970.46335773,
            649423.64089016,
            647425.87204344,
            629386.55333025,
            623896.55036331,
            615171.6634657,
            609770.71698915,
            606896.78477367,
            609770.71698915,
            612860.29930411,
            646628.97921259,
            672653.84034255,
            694949.57011815,
            2053484.48904133,
            2512451.15337624,
            2726859.95574373,
            2779609.28386339,
            2827486.59107239,
            2833326.90480091,
            2830283.03259496,
            2726859.95574373,
            2707367.42825967,
            2557240.92209016,
            2420531.19687943,
            2245241.1128414,
            1829782.22978925,
            743769.25832605,
            719107.79112942,
            710882.47527273,
            708782.9061541,
            706605.83914068,
            702884.83678823,
            699330.14877435,
            698106.79217618,
            697332.03584512,
            696317.93467906,
            696167.64460966,
            695633.55941906,
            694677.60526595,
            693429.24269629,
            692716.93526192,
            692379.41939462,
            690176.6983718,
            688681.83478372,
            687162.86783629,
            686235.36369922,
            686150.90910499,
            685951.34145025,
            685005.95081601,
            684576.59330447,
            683889.43627665,
            687172.58638815,
            1980227.21562606,
            1989715.58702617,
            1992430.35210815,
            1998953.44050118,
            1998160.88128216,
            1972764.03439268,
            1903551.76041575,
            543762.2932594,
            530289.21416049,
            624896.67937685,
            1680052.006021,
            322144.56539297,
            77334.77862527,
            409527.12801416,
        ]
    )

    demographic_events = []
    for sz, t in zip(sizes, times):
        demographic_events.append(
            msprime.PopulationParametersChange(time=t, initial_size=sz, population_id=0)
        )
    populations = [
        stdpopsim.Population(
            id="GAS",
            description="Gabon gambiae population",
        )
    ]

    return stdpopsim.DemographicModel(
        id="GabonAg1000G_1A17",
        description="Stairwayplot estimates of N(t) for Gabon sample",
        long_description="""
        These estimates were done as part of the Ag1000G 2017 Consortium paper.
        Stairwayplot was run with the addition of a
        misorientation parameter using SFS information from each population.
        The model contains 110 distinct epochs, so only some summaries
        are reported in the population table.
        """,
        populations=populations,
        citations=[
            stdpopsim.Citation(
                doi="https://doi.org/10.1038/nature24995",
                year=2017,
                author="Ag1000G Consortium",
                reasons={stdpopsim.CiteReason.DEM_MODEL},
            )
        ],
        # generation times and mutation rate given at bottom of page 32 in supp info,
        # and confirmed with Andy Kern that 3.5e-9 was indeed used for these stairwayplot
        # estimates as well as dadi
        generation_time=1 / 11,
        mutation_rate=3.5e-9,
        demographic_events=demographic_events,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=4069863, metadata=populations[0].asdict()
            )
        ],
    )


_species.add_demographic_model(_GAS_sp())
