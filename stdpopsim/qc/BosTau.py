import msprime
import numpy as np

import stdpopsim


_species = stdpopsim.get_species("BosTau")

population_HolsteinFriesian = [
    stdpopsim.Population("Holstein_Friesian", "Holstein_Friesian")
]


def McLeod2013_1Pop():
    id = "QC-HolsteinFriesian_1M13"
    populations = population_HolsteinFriesian

    # values used by original model (QC-HolsteinFriesian_1M13)
    # mutation_rate = 1.0e-8
    recombination_rate = 1.0e-8

    # mutation rate taken from description of simulations
    # used to generate data mimicing the real genomic
    # data. See page 5-6 in SI of McLeod et al. (2013)
    mutation_rate = 9.4e-9

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
            29800,
            # , 900,000 LAST TIME PERIOD LAST FOREVER IN COAL MODELS
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
        recombination_rate=recombination_rate,
    )


_species.get_demographic_model("HolsteinFriesian_1M13").register_qc(McLeod2013_1Pop())


# Same for all four models in Boitard et al. (2016)
boitard_mut_rate = 1.0e-8


def Boitard2016_Friesian():
    id = "QC-HolsteinFriesian_1B16"
    populations = population_HolsteinFriesian

    # Boitard et al. (2016), page 15
    recombination_rate = 3.66e-9

    mutation_rate = boitard_mut_rate

    # arrays of times and sizes taken from personal communication with Simon Biotard

    times = np.array(
        [
            9,
            24,
            47,
            83,
            139,
            228,
            367,
            584,
            923,
            1455,
            2287,
            3590,
            5629,
            8821,
            13817,
            21638,
            33881,
            53045,
            83043,
            129999,
        ]
    )

    sizes = np.array(
        [
            793,
            1076,
            1320,
            1815,
            3892,
            6908,
            11367,
            15602,
            20159,
            23783,
            27155,
            36512,
            50726,
            52882,
            63092,
            72841,
            71209,
            79079,
            86206,
            65485,
            31652,
        ]
    )

    demographic_events = []
    population_configurations = [msprime.PopulationConfiguration(initial_size=sizes[0])]

    for i, t in enumerate(times):

        demographic_events.append(
            msprime.PopulationParametersChange(
                time=t, initial_size=sizes[i + 1], population_id=0
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
        population_id_map=[{"Holstein_Friesian": 0}] * 21,
        mutation_rate=mutation_rate,
        recombination_rate=recombination_rate,
    )


_species.get_demographic_model("HolsteinFriesian_1B16").register_qc(
    Boitard2016_Friesian()
)


def Boitard2016_Jersey():
    id = "QC-Jersey_1B16"
    populations = [stdpopsim.Population("Jersey", "Jersey")]

    # Boitard et al. (2016), page 15
    recombination_rate = 4.58e-9

    # Same for all four models in Boitard et al. (2016)
    mutation_rate = boitard_mut_rate

    # arrays of times and sizes taken from personal communication with Simon Biotard

    times = np.array(
        [
            9,
            24,
            47,
            83,
            139,
            228,
            367,
            584,
            923,
            1455,
            2287,
            3590,
            5629,
            8821,
            13817,
            21638,
            33881,
            53045,
            83043,
            129999,
        ]
    )

    sizes = np.array(
        [
            388,
            561,
            947,
            1640,
            3704,
            5766,
            8106,
            10474,
            13460,
            13960,
            22549,
            33713,
            51777,
            63277,
            78608,
            77971,
            75739,
            67505,
            76905,
            53856,
            30275,
        ]
    )

    demographic_events = []
    population_configurations = [msprime.PopulationConfiguration(initial_size=sizes[0])]

    for i, t in enumerate(times):

        demographic_events.append(
            msprime.PopulationParametersChange(
                time=t, initial_size=sizes[i + 1], population_id=0
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
        population_id_map=[{"Jersey": 0}] * 21,
        mutation_rate=mutation_rate,
        recombination_rate=recombination_rate,
    )


_species.get_demographic_model("Jersey_1B16").register_qc(Boitard2016_Jersey())


def Boitard2016_Fleckvieh():
    id = "QC-Fleckvieh_1B16"
    populations = [stdpopsim.Population("Fleckvieh", "Fleckvieh")]

    # Boitard et al. (2016), page 15
    recombination_rate = 3.89e-9

    # Same for all four models in Boitard et al. (2016)
    mutation_rate = boitard_mut_rate

    # arrays of times and sizes taken from personal communication with Simon Biotard

    times = np.array(
        [
            9,
            24,
            47,
            83,
            139,
            228,
            367,
            584,
            923,
            1455,
            2287,
            3590,
            5629,
            8821,
            13817,
            21638,
            33881,
            53045,
            83043,
            129999,
        ]
    )

    sizes = np.array(
        [
            2227,
            2812,
            3395,
            4184,
            5689,
            7414,
            10546,
            13604,
            18176,
            20963,
            28547,
            38069,
            54358,
            61479,
            72920,
            82844,
            82448,
            85607,
            91228,
            70789,
            31131,
        ]
    )

    demographic_events = []
    population_configurations = [msprime.PopulationConfiguration(initial_size=sizes[0])]

    for i, t in enumerate(times):
        demographic_events.append(
            msprime.PopulationParametersChange(
                time=t, initial_size=sizes[i + 1], population_id=0
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
        population_id_map=[{"Fleckvieh": 0}] * 21,
        mutation_rate=mutation_rate,
        recombination_rate=recombination_rate,
    )


_species.get_demographic_model("Fleckvieh_1B16").register_qc(Boitard2016_Fleckvieh())

# Angus_1B16


def Boitard2016_Angus():
    id = "QC-Angus_1B16"
    populations = [stdpopsim.Population("Angus", "Angus")]

    # Boitard et al. (2016), page 15
    recombination_rate = 5.00e-9

    # Same for all four models in Boitard et al. (2016)
    mutation_rate = boitard_mut_rate

    # arrays of times and sizes taken from personal communication with Simon Biotard

    times = np.array(
        [
            9,
            24,
            47,
            83,
            139,
            228,
            367,
            584,
            923,
            1455,
            2287,
            3590,
            5629,
            8821,
            13817,
            21638,
            33881,
            53045,
            83043,
            129999,
        ]
    )

    sizes = np.array(
        [
            291,
            412,
            709,
            1224,
            2558,
            4836,
            6837,
            10827,
            12042,
            16441,
            24300,
            35998,
            47892,
            56775,
            69328,
            80651,
            83576,
            87426,
            85228,
            56981,
            33810,
        ]
    )

    demographic_events = []
    population_configurations = [msprime.PopulationConfiguration(initial_size=sizes[0])]
    for i, t in enumerate(times):
        demographic_events.append(
            msprime.PopulationParametersChange(
                time=t, initial_size=sizes[i + 1], population_id=0
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
        population_id_map=[{"Angus": 0}] * 21,
        mutation_rate=mutation_rate,
        recombination_rate=recombination_rate,
    )


_species.get_demographic_model("Angus_1B16").register_qc(Boitard2016_Angus())
