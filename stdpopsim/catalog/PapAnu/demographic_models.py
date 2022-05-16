import msprime
import numpy as np
import stdpopsim

_species = stdpopsim.get_species("PapAnu")


def _pap_anu():
    # the size during the interval times[k] to times[k+1] = sizes[k]
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
    populations = [
        stdpopsim.Population(
            id="PAnubis_SNPRC",
            description="Papio Anubis population from SNPRC",
        )
    ]

    return stdpopsim.DemographicModel(
        id="SinglePopSMCpp_1W22",
        description="SMC++ estimates of N(t) for Papio Anubis individuals",
        long_description="""
        These estimates were obtained from a sample of Papio Anubis
        individuals from the colony housed at the Southwest National
        Primate Research Center (SNPRC). SMC++ was run with a subset of
        36 individuals from the population.
        """,
        populations=populations,
        citations=[
            stdpopsim.Citation(
                doi="https://doi.org/10.1093/gbe/evac040",
                year=2022,
                author="Wall et. al.",
                reasons={stdpopsim.CiteReason.DEM_MODEL},
            )
        ],
        # citations for generation time and mutation rate can be
        # found in species.py
        generation_time=11,
        mutation_rate=5.7e-9,
        demographic_events=demographic_events,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=335505, metadata=populations[0].asdict()
            )
        ],
    )


_species.add_demographic_model(_pap_anu())
