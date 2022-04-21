import msprime
import stdpopsim

_species = stdpopsim.get_species("BosTau")


def _HolsteinFriesian_1M13():
    id = "HolsteinFriesian_1M13"
    description = "Piecewise size changes in Holstein-Friesian cattle."
    long_description = """
    The piecewise-constant population size model of Holstein-Friesian cattle
    from MacLeod et al. (2013). Population sizes were estimated from inferred
    runs of homozygosity, with parameter values taken from Figure 4A and Table S1.
    """
    populations = [
        stdpopsim.Population(id="Holstein_Friesian", description="Holstein-Friesian"),
    ]
    citations = [
        stdpopsim.Citation(
            author="MacLeod et al.",
            year=2013,
            doi="https://doi.org/10.1093/molbev/mst125",
            reasons={stdpopsim.CiteReason.DEM_MODEL},
        )
    ]
    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=_species.generation_time,
        mutation_rate=0.94e-8,
        population_configurations=[
            #      3 generations at    90,     1-     3
            msprime.PopulationConfiguration(
                initial_size=90, metadata=populations[0].asdict()
            )
        ],
        # Here 'time' should be in generation notation ie. how many
        # generations ago were that Ne (effective population size)
        demographic_events=[
            #      3 generations at   120,     4-     6
            msprime.PopulationParametersChange(
                time=4, initial_size=120, population_id=0
            ),
            #      6 generations at   250,     7-    12
            msprime.PopulationParametersChange(
                time=7, initial_size=250, population_id=0
            ),
            #      6 generations at   350,    13-    18
            msprime.PopulationParametersChange(
                time=13, initial_size=350, population_id=0
            ),
            #      6 generations at  1000,    19-    24
            msprime.PopulationParametersChange(
                time=19, initial_size=1000, population_id=0
            ),
            #    130 generations at  1500,    25-   154
            msprime.PopulationParametersChange(
                time=25, initial_size=1500, population_id=0
            ),
            #    200 generations at  2000,   155-   454
            msprime.PopulationParametersChange(
                time=155, initial_size=2000, population_id=0
            ),
            #    200 generations at  2500,   455-   654
            msprime.PopulationParametersChange(
                time=455, initial_size=2500, population_id=0
            ),
            #   1100 generations at  3500,   655-  1754
            msprime.PopulationParametersChange(
                time=655, initial_size=3500, population_id=0
            ),
            #    600 generations at  7000,  1755-  2354
            msprime.PopulationParametersChange(
                time=1755, initial_size=7000, population_id=0
            ),
            #   1000 generations at 10000,  2355-  3354
            msprime.PopulationParametersChange(
                time=2355, initial_size=10000, population_id=0
            ),
            #  29800 generations at 17000,  3355- 33154
            msprime.PopulationParametersChange(
                time=3355, initial_size=17000, population_id=0
            ),
            # 900000 generations at 62000, 33155-933154
            msprime.PopulationParametersChange(
                time=33155, initial_size=62000, population_id=0
            ),
        ],
    )


_species.add_demographic_model(_HolsteinFriesian_1M13())
