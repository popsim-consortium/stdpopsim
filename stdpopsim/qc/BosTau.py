import msprime
import numpy as np

import stdpopsim


_species = stdpopsim.get_species("BosTau")


def McLeod2013_1Pop():
    id = "QC-HolsteinFriesian_1M13"
    populations = [stdpopsim.Population("Holstein_Friesian", "Holstein_Friesian")]

    # mutation rate taken from description of simulations
    # used to generate data mimicing the real genomic
    # data. See page 5-6 in SI of McLeod et al. (2013)
    mutation_rate = 9.4e-9

    # value used by original model (QC-HolsteinFriesian_1M13)
    # mutation_rate = 1.0e-8

    # Arrays of time durations and Ne's taken from
    # Supp Tabls S1 of McLeod et al. (2013).
    # These are values estimated in original analysis
    # (Figure 4A in paper) and used in validation simulations

    times = np.array(
        [
            # 4, Original model (QC-HolsteinFriesian_1M13) adds a 1 generation epoch
            3,
            3,
            6,
            6,
            6,
            130,
            300,
            200,
            1100,
            600,
            1000,
            29800
            # , 900,000 LAST TIME PERIOD LAST FOEVER IN COAL MODELS
        ]
    )

    sizes = np.array(
        [90, 120, 250, 350, 1000, 1500, 2000, 2500, 3500, 7000, 10000, 17000, 62000]
    )

    demographic_events = []
    population_configurations = [msprime.PopulationConfiguration(initial_size=sizes[0])]
    curr_time = 0
    for i, t in enumerate(times):
        curr_time += t
        demographic_events.append(
            msprime.PopulationParametersChange(
                time=curr_time, initial_size=sizes[i + 1], population_id=0
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
        population_id_map=[{"Holstein_Friesian": 0}] * 13,
        mutation_rate=mutation_rate,
    )


_species.get_demographic_model("HolsteinFriesian_1M13").register_qc(McLeod2013_1Pop)
