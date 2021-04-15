import msprime

import stdpopsim

_species = stdpopsim.get_species("AnaPla")


def _mallard_black_split():
    id = "MallardBlackDuck_2L19"
    description = "North American Mallard/Black Duck split"
    long_description = """
        This is a model fit to contemporary samples of wild North American
        mallard and black duck, using the "split-migration" model of dadi.
        See Figure 6 of Lavretsky et al 2019.
    """
    T = 632305 / 4  # in generations, not years
    N_BlackDuck = 1.57e6
    N_Mallard = 1.37e6
    # personal communication from Joshua Brown 13 Apr 2021:
    # "Based on the contemporary dataset, the ancestral population size
    # "for the Black duck/Mallard dadi model was 819535."
    N_Anc = 819535
    # the migration rate is reported as 2.82 in each direction.  From the dadi
    # manual, m12 is "the fraction of individuals each generation in pop 1 that
    # are new migrants from pop 2, times the 2Nref". To convert back to real
    # time units (fraction replaced per generation) we divide by 2 * N_anc.
    m = 2.82 / (2 * N_Anc)

    model = msprime.Demography()
    model.add_population(
        initial_size=N_Mallard,
        name="Mallard",
        description="Wild North American mallards",
    )
    model.add_population(
        initial_size=N_BlackDuck,
        name="Black_duck",
        description="Wild black ducks",
    )
    model.add_population(
        initial_size=N_Anc,
        name="Ancestral",
        description="Ancestral population",
    )
    model.set_symmetric_migration_rate(populations=["Mallard", "Black_duck"], rate=m)

    model.add_population_split(
        time=T, derived=["Mallard", "Black_duck"], ancestral="Ancestral"
    )

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        citations=[
            stdpopsim.Citation(
                author="Lavretsky et al.",
                year=2019,
                doi="https://doi.org/10.1111/mec.15343",
                reasons={stdpopsim.CiteReason.DEM_MODEL},
            )
        ],
        generation_time=4,
        model=model,
        mutation_rate=4.83e-9,
    )


_species.add_demographic_model(_mallard_black_split())
