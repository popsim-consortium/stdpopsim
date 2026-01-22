import msprime

import stdpopsim

_species = stdpopsim.get_species("AnoCar")


def _GA_vs_EF_IMex():
    id = "Anoli_5B19"
    description = "Anoli GA-EK split"
    long_description = """
        This is a model fit to whole genomes of Anolis carolinesis
        using five populations and fitting a model of
        isolation with constant migration and no interruption
        of gene flow.
        See Figure 3 and table 2 Bourgeois et al 2019.
    """
    T = 1  # in generations
    N_GulfAtlantic = 364329
    N_EasternFlorida = 6755470
    N_Anc = 2156641
    # the migration rate
    m_ge = 2.36e-07
    m_eg = 1.93e-07
    # Tsc, time during which stable populations stay connected
    # Tscg, time since population size change (with gene flow)
    # Tsc = 61217
    # Tscg = 1962652
    # after inital split initial populations have a population
    # size of N_GulfAtlantic_a and N_EasternFlorida_a
    # N_GulfAtlantic_a = 5971380
    # N_EasternFlorida_a = 2137291

    model = msprime.Demography()
    model.add_population(
        initial_size=N_GulfAtlantic,
        name="Gulf_Atlantic",
        description="Anoles from Gulf Atlantic",
    )
    model.add_population(
        initial_size=N_EasternFlorida,
        name="Eastern_Florida",
        description="Anoles from Eastern Florida",
    )
    model.add_population(
        initial_size=N_Anc,
        name="Ancestral",
        description="Ancestral population",
    )
    model.set_migration_rate(rate=m_ge, source="Gulf_Atlantic", dest="Eastern_Florida")
    model.set_migration_rate(rate=m_eg, source="Eastern_Florida", dest="Gulf_Atlantic")

    model.add_population_split(
        time=T, derived=["Gulf_Atlantic", "Eastern_Florida"], ancestral="Ancestral"
    )

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        citations=[
            stdpopsim.Citation(
                author="Bourgeois et al.",
                year=2019,
                doi="https://doi.org/10.1093/gbe/evz110",
                reasons={stdpopsim.CiteReason.DEM_MODEL},
            )
        ],
        generation_time=1.5,
        model=model,
        # mutation_rate=2.1e-10,
    )


_species.add_demographic_model(_GA_vs_EF_IMex)
