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
            #       3 generations at     90,      [0-      3)
            msprime.PopulationConfiguration(
                initial_size=90, metadata=populations[0].asdict()
            )
        ],
        # Here 'time' should be in generation notation ie. how many
        # generations ago were that Ne (effective population size)
        # Note the zero-based counting in the code
        demographic_events=[
            #       3 generations at    120,      [3-      6)
            msprime.PopulationParametersChange(
                time=3, initial_size=120, population_id=0
            ),
            #       6 generations at    250,      [6-     12)
            msprime.PopulationParametersChange(
                time=6, initial_size=250, population_id=0
            ),
            #       6 generations at    350,     [12-     18)
            msprime.PopulationParametersChange(
                time=12, initial_size=350, population_id=0
            ),
            #       6 generations at  1,000,     [18-     24)
            msprime.PopulationParametersChange(
                time=18, initial_size=1000, population_id=0
            ),
            #     130 generations at  1,500,     [24-    154)
            msprime.PopulationParametersChange(
                time=24, initial_size=1500, population_id=0
            ),
            #     300 generations at  2,000,    [154-    454)
            msprime.PopulationParametersChange(
                time=154, initial_size=2000, population_id=0
            ),
            #     200 generations at  2,500,    [454-    654)
            msprime.PopulationParametersChange(
                time=454, initial_size=2500, population_id=0
            ),
            #   1,100 generations at  3,500,    [654-  1,754)
            msprime.PopulationParametersChange(
                time=654, initial_size=3500, population_id=0
            ),
            #     600 generations at  7,000,  [1,754-  2,354)
            msprime.PopulationParametersChange(
                time=1754, initial_size=7000, population_id=0
            ),
            #   1,000 generations at 10,000,  [2,354-  3,354)
            msprime.PopulationParametersChange(
                time=2354, initial_size=10000, population_id=0
            ),
            #  29,800 generations at 17 000,  [3,354- 33,154)
            msprime.PopulationParametersChange(
                time=3354, initial_size=17000, population_id=0
            ),
            # 900,000 generations at 62,000, [33,154-93,3154)
            msprime.PopulationParametersChange(
                time=33154, initial_size=62000, population_id=0
            ),
        ],
    )


_species.add_demographic_model(_HolsteinFriesian_1M13())
