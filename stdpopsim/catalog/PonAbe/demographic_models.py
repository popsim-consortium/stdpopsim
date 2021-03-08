import math

import msprime

import stdpopsim


_species = stdpopsim.get_species("PonAbe")


_locke2011 = stdpopsim.Citation(
    author="Locke et al.", year=2011, doi="http://doi.org/10.1038/nature09687"
)


def _orangutan():
    id = "TwoSpecies_2L11"
    description = "Two population orangutan model"
    long_description = """
        The two orang-utan species, Sumatran (Pongo abelii) and Bornean (Pongo
        pygmaeus) inferred from the joint-site frequency spectrum with ten
        individuals from each population. This model is an isolation-with-
        migration model, with exponential growth or decay in each population
        after the split. The Sumatran population grows in size, while the
        Bornean population slightly declines.
    """

    citations = [_locke2011.because(stdpopsim.CiteReason.DEM_MODEL)]

    populations = [
        stdpopsim.Population("Bornean", "Pongo pygmaeus (Bornean) population"),
        stdpopsim.Population("Sumatran", "Pongo abelii (Sumatran) population"),
    ]

    # Parameters from paper:
    # ancestral size, before split
    Na = 17934

    # time of split
    T_split_years = 403149
    # get split time in units of generations
    generation_time = _species.generation_time
    T_split = T_split_years / generation_time

    # proportion of ancestral pop to branch B
    s = 0.592

    # sizes at split
    Na_B = 17934 * s
    Na_S = 17934 * (1 - s)

    # present sizes
    N_B = 8805
    N_S = 37661

    # get growth rates
    r_B = -1 * math.log(Na_B / N_B) / T_split
    r_S = -1 * math.log(Na_S / N_S) / T_split

    # migration rates
    m_S_B = 0.395 / 2 / Na
    m_B_S = 0.239 / 2 / Na

    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=B and 1=S
    # initially.

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        citations=citations,
        populations=populations,
        generation_time=generation_time,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_B, growth_rate=r_B, metadata=populations[0].asdict()
            ),  # NOQA
            msprime.PopulationConfiguration(
                initial_size=N_S, growth_rate=r_S, metadata=populations[1].asdict()
            ),  # NOQA
        ],
        migration_matrix=[
            [0, m_B_S],  # NOQA
            [m_S_B, 0],  # NOQA
        ],
        demographic_events=[
            # Merge populations and turn off migration, change to size Na
            msprime.MassMigration(
                time=T_split, source=1, destination=0, proportion=1.0
            ),
            msprime.MigrationRateChange(time=T_split, rate=0),
            msprime.PopulationParametersChange(
                time=T_split, initial_size=Na, growth_rate=0, population_id=0
            ),
        ],
    )


_species.add_demographic_model(_orangutan())
