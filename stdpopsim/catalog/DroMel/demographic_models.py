import msprime

import stdpopsim

_species = stdpopsim.get_species("DroMel")


_LiAndStephan = stdpopsim.Citation(
    author="Li et al.", year=2006, doi="https://doi.org/10.1371/journal.pgen.0020166"
)

# population definitions that are reused.
_afr_population = stdpopsim.Population(
    id="AFR", description="African D. melanogaster population"
)
_eur_population = stdpopsim.Population(
    id="EUR", description="European D. melanogaster population"
)


def _afr_3epoch():
    id = "African3Epoch_1S16"
    description = "Three epoch African population"
    long_description = """
        The three epoch (modern, bottleneck, ancestral) model estimated for a
        single African Drosophila Melanogaster population from Sheehan and Song (2016).
        Population sizes are estimated by a
        deep learning model trained on simulation data. NOTE: Due to differences in
        coalescence units between PSMC (2N) and msms (4N) the number of generations were
        doubled from PSMC estimates when simulating data from msms in the original
        publication. We have faithfully represented the published model here.
    """
    populations = [_afr_population]
    citations = [
        stdpopsim.Citation(
            author="Sheehan and Song",
            year=2016,
            doi="https://doi.org/10.1371/journal.pcbi.1004845",
            reasons={stdpopsim.CiteReason.DEM_MODEL},
        )
    ]
    generation_time = _species.generation_time

    # Parameter values from "Simulating Data" section
    # these are assumptions, not estimates
    N_ref = 100000
    t_1_coal = 0.5
    t_2_coal = 5.0
    # estimates from the ANN
    N_R = 544200
    N_B = 145300
    N_A = 652700
    # Times are provided in 4N_ref generations, so we convert into generations.
    # generation_time = 10 / year
    t_1 = t_1_coal * 4 * N_ref
    t_2 = (t_1_coal + t_2_coal) * 4 * N_ref

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_R, metadata=populations[0].asdict()
            )
        ],
        demographic_events=[
            # Size change at bottleneck (back in time; BIT)
            msprime.PopulationParametersChange(
                time=t_1, initial_size=N_B, population_id=0
            ),
            # Size change at recovery (BIT)
            msprime.PopulationParametersChange(
                time=t_2, initial_size=N_A, population_id=0
            ),
        ],
    )


_species.add_demographic_model(_afr_3epoch())


def _ooa_2():
    id = "OutOfAfrica_2L06"
    description = "Three epoch model for African and European populations"
    long_description = """
        The three epoch (modern, bottleneck, ancestral) model estimated for two
        Drosophila Melanogaster populations: African (ancestral) and European (derived)
        from Li and Stephan (2006).
    """
    populations = [_afr_population, _eur_population]
    citations = [_LiAndStephan.because(stdpopsim.CiteReason.DEM_MODEL)]
    generation_time = _species.generation_time

    # African Parameter values from "Demographic History of the African
    # Population" section
    N_A0 = 8.603e06
    t_A0 = 600000  # assuming 10 generations / year
    N_A1 = N_A0 / 5.0

    # European Parameter values from "Demo History of Euro Population"
    N_E0 = 1.075e06
    N_E1 = 2200
    t_AE = 158000  # generations
    t_E1 = t_AE - 3400

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_A0, metadata=populations[0].asdict()
            ),
            msprime.PopulationConfiguration(
                initial_size=N_E0, metadata=populations[1].asdict()
            ),
        ],
        demographic_events=[
            # Size change at Euro bottleneck
            msprime.PopulationParametersChange(
                time=t_E1, initial_size=N_E1, population_id=1
            ),
            # Split
            msprime.MassMigration(time=t_AE, source=1, destination=0, proportion=1.0),
            # African bottleneck
            msprime.PopulationParametersChange(
                time=t_A0, initial_size=N_A1, population_id=0
            ),
        ],
    )


_species.add_demographic_model(_ooa_2())
