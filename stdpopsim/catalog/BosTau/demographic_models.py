import msprime
import stdpopsim

_species = stdpopsim.get_species("BosTau")


def _HolsteinFriesian_1M13():
    id = "HolsteinFriesian_1M13"
    description = (
        "Piecewise size model for Holstein-Friesian cattle (MacLeod et al., 2013)."
    )
    long_description = """
    The piecewise-constant population size model of Holstein-Friesian cattle
    from MacLeod et al. (2013). Effective population sizes were estimated based
    on runs of homozygosity observed in a single individual,
    using the following assumptions:
    a generation interval of 5 years (page 2213),
    a mutation rate of 0.94e-8 (pages 2221 and 2215), and
    a recombination rate of 1e-8 (pages 2215 and 2221).
    Effective population sizes are given in Figure 4A and Table S1.
    The single individual is the bull Walkway Chief Mark.
    Mark and his father are two of the most influentiual sires in this breed - see
    Larkin et al. (2012) http://www.pnas.org/cgi/doi/10.1073/pnas.1114546109.
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
        recombination_rate=1e-8,
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
            #    many generations at 62,000, [33,154-93,3154)
            msprime.PopulationParametersChange(
                time=33154, initial_size=62000, population_id=0
            ),
        ],
    )


_species.add_demographic_model(_HolsteinFriesian_1M13())


def _HolsteinFriesian_1B16():
    id = "HolsteinFriesian_1B16"
    description = (
        "Piecewise size model for Holstein-Friesian cattle (Boitard et al., 2016)."
    )
    long_description = """
    The piecewise-constant population size model of Holstein-Friesian cattle
    from Boitard et al. (2016). Effective population sizes were estimated using
    Approximate Bayesian Computation with allele frequency spectrum and
    linkage-disequlibrium statistics (both on SNPs with minor allele frequency
    above 0.2) observed in 25 individuals,
    using the following assumptions:
    a generation interval of 5 years (page 10),
    a mutation rate of 1e-8 (pages 8 and 22), and
    a recombination rate of 3.66e-9 (page 15).
    Effective population sizes are given in Figure 6,
    with exact values provided by the lead author in personal communication.
    The 25 individuals' genomes were obtained from the 1000 bull genomes Run 2
    (of which 129 were Holstein-Friesian;
    the authors chose 52 from Australia to get a homogeneous sub-population,
    kept the least inbred and unrelated samples, and from these randomly selected 25).
    """
    populations = [
        stdpopsim.Population(id="Holstein_Friesian", description="Holstein-Friesian"),
    ]
    citations = [
        stdpopsim.Citation(
            author="Boitard et al.",
            year=2016,
            doi="https://doi.org/10.1371/journal.pgen.1005877",
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
        mutation_rate=1e-8,
        recombination_rate=3.66e-9,
        population_configurations=[
            #     9 generations at    793,     [0-     9)      9-    0=    9
            msprime.PopulationConfiguration(
                initial_size=793, metadata=populations[0].asdict()
            )
        ],
        # Here 'time' should be in generation notation ie. how many
        # generations ago were that Ne (effective population size)
        # Note the zero-based counting in the code
        demographic_events=[
            #    15 generations at   1076,     [9-    24)     24-    9=   15
            msprime.PopulationParametersChange(
                time=9, initial_size=1076, population_id=0
            ),
            #    23 generations at   1320,    [24-    47)     47-   24=   23
            msprime.PopulationParametersChange(
                time=24, initial_size=1320, population_id=0
            ),
            #    36 generations at  1815,    [47-    83)     83-   47=   36
            msprime.PopulationParametersChange(
                time=47, initial_size=1815, population_id=0
            ),
            #    56 generations at  3892,    [83-   139)    139-   83=   56
            msprime.PopulationParametersChange(
                time=83, initial_size=3892, population_id=0
            ),
            #    89 generations at  6908,   [139-   228)    228-  139=   89
            msprime.PopulationParametersChange(
                time=139, initial_size=6908, population_id=0
            ),
            #   139 generations at 11367,   [228-   367)    367-  228=  139
            msprime.PopulationParametersChange(
                time=228, initial_size=11367, population_id=0
            ),
            #   217 generations at 15602,   [367-   584)    584-  367=  217
            msprime.PopulationParametersChange(
                time=367, initial_size=15602, population_id=0
            ),
            #   339 generations at 20159,   [584-   923)    923-  584=  339
            msprime.PopulationParametersChange(
                time=584, initial_size=20159, population_id=0
            ),
            #   532 generations at 23783,   [923-  1455)   1455-  923=  532
            msprime.PopulationParametersChange(
                time=923, initial_size=23783, population_id=0
            ),
            #   832 generations at 27155,  [1455-  2287)   2287- 1455=  832
            msprime.PopulationParametersChange(
                time=1455, initial_size=27155, population_id=0
            ),
            #  1303 generations at 36512,  [2287-  3590)   3590- 2287= 1303
            msprime.PopulationParametersChange(
                time=2287, initial_size=36512, population_id=0
            ),
            #  2039 generations at 50726,  [3590-  5629)   5629- 3590= 2039
            msprime.PopulationParametersChange(
                time=3590, initial_size=50726, population_id=0
            ),
            #  3192 generations at 52882,  [5629-  8821)   8821- 5629= 3192
            msprime.PopulationParametersChange(
                time=5629, initial_size=52882, population_id=0
            ),
            #  4996 generations at 63092,  [8821- 13817)  13817- 8821= 4996
            msprime.PopulationParametersChange(
                time=8821, initial_size=63092, population_id=0
            ),
            #  7821 generations at 72841, [13817- 21638)  21638-13817= 7821
            msprime.PopulationParametersChange(
                time=13817, initial_size=72841, population_id=0
            ),
            # 12243 generations at 71209, [21638- 33881)  33881-21638=12243
            msprime.PopulationParametersChange(
                time=21638, initial_size=71209, population_id=0
            ),
            # 19164 generations at 79079, [33881- 53045)  53045-33881=19164
            msprime.PopulationParametersChange(
                time=33881, initial_size=79079, population_id=0
            ),
            # 29998 generations at 86206, [53045- 83043)  83043-53045=29998
            msprime.PopulationParametersChange(
                time=53045, initial_size=86206, population_id=0
            ),
            # 46956 generations at 65485, [83043-129999) 129999-83043=46956
            msprime.PopulationParametersChange(
                time=83043, initial_size=65485, population_id=0
            ),
            # "many" generations at 31652, [129999-inf) inf-129999="many"
            msprime.PopulationParametersChange(
                time=129999, initial_size=31652, population_id=0
            ),
        ],
    )


_species.add_demographic_model(_HolsteinFriesian_1B16())


def _Fleckvieh_1B16():
    id = "Fleckvieh_1B16"
    description = "Piecewise size model for Fleckvieh cattle (Boitard et al., 2016)."
    long_description = """
    The piecewise-constant population size model of Fleckvieh cattle
    from Boitard et al. (2016). Effective population sizes were estimated using
    Approximate Bayesian Computation with allele frequency spectrum and
    linkage-disequlibrium statistics (both on SNPs with minor allele frequency
    above 0.2) observed in 25 individuals,
    using the following assumptions:
    a generation interval of 5 years (page 10),
    a mutation rate of 1e-8 (pages 8 and 22), and
    a recombination rate of 3.89e-9 (page 15).
    Effective population sizes are given in Figure 6,
    with exact values provided by the lead author in personal communication.
    The 25 individuals' genomes were obtained from the 1000 bull genomes Run 2
    (of which 43 were Fleckvieh;
    the authors kept the least inbred and unrelated samples, and
    from these randomly selected 25).
    """
    populations = [
        stdpopsim.Population(id="Fleckvieh", description="Fleckvieh"),
    ]
    citations = [
        stdpopsim.Citation(
            author="Boitard et al.",
            year=2016,
            doi="https://doi.org/10.1371/journal.pgen.1005877",
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
        mutation_rate=1e-8,
        recombination_rate=3.89e-9,
        population_configurations=[
            #     9 generations at  2227,     [0-     9)      9-    0=    9
            msprime.PopulationConfiguration(
                initial_size=2227, metadata=populations[0].asdict()
            )
        ],
        # Here 'time' should be in generation notation ie. how many
        # generations ago were that Ne (effective population size)
        # Note the zero-based counting in the code
        demographic_events=[
            #    15 generations at  2812,     [9-    24)     24-    9=   15
            msprime.PopulationParametersChange(
                time=9, initial_size=2812, population_id=0
            ),
            #    23 generations at  3395,    [24-    47)     47-   24=   23
            msprime.PopulationParametersChange(
                time=24, initial_size=3395, population_id=0
            ),
            #    36 generations at  4184,    [47-    83)     83-   47=   36
            msprime.PopulationParametersChange(
                time=47, initial_size=4184, population_id=0
            ),
            #    56 generations at  5689,    [83-   139)    139-   83=   56
            msprime.PopulationParametersChange(
                time=83, initial_size=5689, population_id=0
            ),
            #    89 generations at  7414,   [139-   228)    228-  139=   89
            msprime.PopulationParametersChange(
                time=139, initial_size=7414, population_id=0
            ),
            #   139 generations at 10546,   [228-   367)    367-  228=  139
            msprime.PopulationParametersChange(
                time=228, initial_size=10546, population_id=0
            ),
            #   217 generations at 13604,   [367-   584)    584-  367=  217
            msprime.PopulationParametersChange(
                time=367, initial_size=13604, population_id=0
            ),
            #   339 generations at 18176,   [584-   923)    923-  584=  339
            msprime.PopulationParametersChange(
                time=584, initial_size=18176, population_id=0
            ),
            #   532 generations at 20963,   [923-  1455)   1455-  923=  532
            msprime.PopulationParametersChange(
                time=923, initial_size=20963, population_id=0
            ),
            #   832 generations at 28547,  [1455-  2287)   2287- 1455=  832
            msprime.PopulationParametersChange(
                time=1455, initial_size=28547, population_id=0
            ),
            #  1303 generations at 38069,  [2287-  3590)   3590- 2287= 1303
            msprime.PopulationParametersChange(
                time=2287, initial_size=38069, population_id=0
            ),
            #  2039 generations at 54358,  [3590-  5629)   5629- 3590= 2039
            msprime.PopulationParametersChange(
                time=3590, initial_size=54358, population_id=0
            ),
            #  3192 generations at 61479,  [5629-  8821)   8821- 5629= 3192
            msprime.PopulationParametersChange(
                time=5629, initial_size=61479, population_id=0
            ),
            #  4996 generations at 72920,  [8821- 13817)  13817- 8821= 4996
            msprime.PopulationParametersChange(
                time=8821, initial_size=72920, population_id=0
            ),
            #  7821 generations at 82844, [13817- 21638)  21638-13817= 7821
            msprime.PopulationParametersChange(
                time=13817, initial_size=82844, population_id=0
            ),
            # 12243 generations at 82448, [21638- 33881)  33881-21638=12243
            msprime.PopulationParametersChange(
                time=21638, initial_size=82448, population_id=0
            ),
            # 19164 generations at 85607, [33881- 53045)  53045-33881=19164
            msprime.PopulationParametersChange(
                time=33881, initial_size=85607, population_id=0
            ),
            # 29998 generations at 91228, [53045- 83043)  83043-53045=29998
            msprime.PopulationParametersChange(
                time=53045, initial_size=91228, population_id=0
            ),
            # 46956 generations at 70789, [83043-129999) 129999-83043=46956
            msprime.PopulationParametersChange(
                time=83043, initial_size=70789, population_id=0
            ),
            #  many generations at 31131, [129999-inf) inf-129999="many"
            msprime.PopulationParametersChange(
                time=129999, initial_size=31131, population_id=0
            ),
        ],
    )


_species.add_demographic_model(_Fleckvieh_1B16())


def _Jersey_1B16():
    id = "Jersey_1B16"
    description = "Piecewise size model for Jersey cattle (Boitard et al., 2016)."
    long_description = """
    The piecewise-constant population size model of Jersey cattle
    from Boitard et al. (2016). Effective population sizes were estimated using
    Approximate Bayesian Computation with allele frequency spectrum and
    linkage-disequlibrium statistics (both on SNPs with minor allele frequency
    above 0.2) observed in 15 individuals,
    using the following assumptions:
    a generation interval of 5 years (page 10),
    a mutation rate of 1e-8 (pages 8 and 22), and
    a recombination rate of 4.58e-9 (page 15).
    Effective population sizes are given in Figure 6,
    with exact values provided by the lead author in personal communication.
    The 15 individuals' genomes were obtained from the 1000 bull genomes Run 2.
    """
    populations = [
        stdpopsim.Population(id="Jersey", description="Jersey"),
    ]
    citations = [
        stdpopsim.Citation(
            author="Boitard et al.",
            year=2016,
            doi="https://doi.org/10.1371/journal.pgen.1005877",
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
        mutation_rate=1e-8,
        recombination_rate=4.58e-9,
        population_configurations=[
            #     9 generations at   388,     [0-     9)      9-    0=    9
            msprime.PopulationConfiguration(
                initial_size=388, metadata=populations[0].asdict()
            )
        ],
        # Here 'time' should be in generation notation ie. how many
        # generations ago were that Ne (effective population size)
        # Note the zero-based counting in the code
        demographic_events=[
            #    15 generations at   561,     [9-    24)     24-    9=   15
            msprime.PopulationParametersChange(
                time=9, initial_size=561, population_id=0
            ),
            #    23 generations at   947,    [24-    47)     47-   24=   23
            msprime.PopulationParametersChange(
                time=24, initial_size=947, population_id=0
            ),
            #    36 generations at  1640,    [47-    83)     83-   47=   36
            msprime.PopulationParametersChange(
                time=47, initial_size=1640, population_id=0
            ),
            #    56 generations at  3704,    [83-   139)    139-   83=   56
            msprime.PopulationParametersChange(
                time=83, initial_size=3704, population_id=0
            ),
            #    89 generations at  5766,   [139-   228)    228-  139=   89
            msprime.PopulationParametersChange(
                time=139, initial_size=5766, population_id=0
            ),
            #   139 generations at  8106,   [228-   367)    367-  228=  139
            msprime.PopulationParametersChange(
                time=228, initial_size=8106, population_id=0
            ),
            #   217 generations at 10474,   [367-   584)    584-  367=  217
            msprime.PopulationParametersChange(
                time=367, initial_size=10474, population_id=0
            ),
            #   339 generations at 13460,   [584-   923)    923-  584=  339
            msprime.PopulationParametersChange(
                time=584, initial_size=13460, population_id=0
            ),
            #   532 generations at 13960,   [923-  1455)   1455-  923=  532
            msprime.PopulationParametersChange(
                time=923, initial_size=13960, population_id=0
            ),
            #   832 generations at 22549,  [1455-  2287)   2287- 1455=  832
            msprime.PopulationParametersChange(
                time=1455, initial_size=22549, population_id=0
            ),
            #  1303 generations at 33713,  [2287-  3590)   3590- 2287= 1303
            msprime.PopulationParametersChange(
                time=2287, initial_size=33713, population_id=0
            ),
            #  2039 generations at 51777,  [3590-  5629)   5629- 3590= 2039
            msprime.PopulationParametersChange(
                time=3590, initial_size=51777, population_id=0
            ),
            #  3192 generations at 63277,  [5629-  8821)   8821- 5629= 3192
            msprime.PopulationParametersChange(
                time=5629, initial_size=63277, population_id=0
            ),
            #  4996 generations at 78608,  [8821- 13817)  13817- 8821= 4996
            msprime.PopulationParametersChange(
                time=8821, initial_size=78608, population_id=0
            ),
            #  7821 generations at 77971, [13817- 21638)  21638-13817= 7821
            msprime.PopulationParametersChange(
                time=13817, initial_size=77971, population_id=0
            ),
            # 12243 generations at 75739, [21638- 33881)  33881-21638=12243
            msprime.PopulationParametersChange(
                time=21638, initial_size=75739, population_id=0
            ),
            # 19164 generations at 67505, [33881- 53045)  53045-33881=19164
            msprime.PopulationParametersChange(
                time=33881, initial_size=67505, population_id=0
            ),
            # 29998 generations at 76905, [53045- 83043)  83043-53045=29998
            msprime.PopulationParametersChange(
                time=53045, initial_size=76905, population_id=0
            ),
            # 46956 generations at 53856, [83043-129999) 129999-83043=46956
            msprime.PopulationParametersChange(
                time=83043, initial_size=53856, population_id=0
            ),
            #  many generations at 30275, [129999-inf) inf-129999="many"
            msprime.PopulationParametersChange(
                time=129999, initial_size=30275, population_id=0
            ),
        ],
    )


_species.add_demographic_model(_Jersey_1B16())


def _Angus_1B16():
    id = "Angus_1B16"
    description = "Piecewise size model for Angus cattle (Boitard et al., 2016)."
    long_description = """
    The piecewise-constant population size model of Angus cattle
    from Boitard et al. (2016). Effective population sizes were estimated using
    Approximate Bayesian Computation with allele frequency spectrum and
    linkage-disequlibrium statistics (both on SNPs with minor allele frequency
    above 0.2) observed in 15 individuals,
    using the following assumptions:
    a generation interval of 5 years (page 10),
    a mutation rate of 1e-8 (pages 8 and 22), and
    a recombination rate of 5.00e-9 (page 15).
    Effective population sizes are given in Figure 6,
    with exact values provided by the lead author in personal communication.
    The 25 individuals' genomes were obtained from the 1000 bull genomes Run 2
    (of which 47 were Angus;
    the authors kept the least inbred and unrelated samples, and
    from these randomly selected 25).
    """
    populations = [
        stdpopsim.Population(id="Angus", description="Angus"),
    ]
    citations = [
        stdpopsim.Citation(
            author="Boitard et al.",
            year=2016,
            doi="https://doi.org/10.1371/journal.pgen.1005877",
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
        mutation_rate=1e-8,
        recombination_rate=5.00e-9,
        population_configurations=[
            #     9 generations at   291,     [0-     9)      9-    0=    9
            msprime.PopulationConfiguration(
                initial_size=291, metadata=populations[0].asdict()
            )
        ],
        # Here 'time' should be in generation notation ie. how many
        # generations ago were that Ne (effective population size)
        # Note the zero-based counting in the code
        demographic_events=[
            #    15 generations at   412,     [9-    24)     24-    9=   15
            msprime.PopulationParametersChange(
                time=9, initial_size=412, population_id=0
            ),
            #    23 generations at   709,    [24-    47)     47-   24=   23
            msprime.PopulationParametersChange(
                time=24, initial_size=709, population_id=0
            ),
            #    36 generations at  1224,    [47-    83)     83-   47=   36
            msprime.PopulationParametersChange(
                time=47, initial_size=1224, population_id=0
            ),
            #    56 generations at  2558,    [83-   139)    139-   83=   56
            msprime.PopulationParametersChange(
                time=83, initial_size=2558, population_id=0
            ),
            #    89 generations at  4836,   [139-   228)    228-  139=   89
            msprime.PopulationParametersChange(
                time=139, initial_size=4836, population_id=0
            ),
            #   139 generations at  6837,   [228-   367)    367-  228=  139
            msprime.PopulationParametersChange(
                time=228, initial_size=6837, population_id=0
            ),
            #   217 generations at 10827,   [367-   584)    584-  367=  217
            msprime.PopulationParametersChange(
                time=367, initial_size=10827, population_id=0
            ),
            #   339 generations at 12042,   [584-   923)    923-  584=  339
            msprime.PopulationParametersChange(
                time=584, initial_size=12042, population_id=0
            ),
            #   532 generations at 16441,   [923-  1455)   1455-  923=  532
            msprime.PopulationParametersChange(
                time=923, initial_size=16441, population_id=0
            ),
            #   832 generations at 24300,  [1455-  2287)   2287- 1455=  832
            msprime.PopulationParametersChange(
                time=1455, initial_size=24300, population_id=0
            ),
            #  1303 generations at 35998,  [2287-  3590)   3590- 2287= 1303
            msprime.PopulationParametersChange(
                time=2287, initial_size=35998, population_id=0
            ),
            #  2039 generations at 47892,  [3590-  5629)   5629- 3590= 2039
            msprime.PopulationParametersChange(
                time=3590, initial_size=47892, population_id=0
            ),
            #  3192 generations at 56775,  [5629-  8821)   8821- 5629= 3192
            msprime.PopulationParametersChange(
                time=5629, initial_size=56775, population_id=0
            ),
            #  4996 generations at 69328,  [8821- 13817)  13817- 8821= 4996
            msprime.PopulationParametersChange(
                time=8821, initial_size=69328, population_id=0
            ),
            #  7821 generations at 80651, [13817- 21638)  21638-13817= 7821
            msprime.PopulationParametersChange(
                time=13817, initial_size=80651, population_id=0
            ),
            # 12243 generations at 83576, [21638- 33881)  33881-21638=12243
            msprime.PopulationParametersChange(
                time=21638, initial_size=83576, population_id=0
            ),
            # 19164 generations at 87426, [33881- 53045)  53045-33881=19164
            msprime.PopulationParametersChange(
                time=33881, initial_size=87426, population_id=0
            ),
            # 29998 generations at 85228, [53045- 83043)  83043-53045=29998
            msprime.PopulationParametersChange(
                time=53045, initial_size=85228, population_id=0
            ),
            # 46956 generations at 56981, [83043-129999) 129999-83043=46956
            msprime.PopulationParametersChange(
                time=83043, initial_size=56981, population_id=0
            ),
            #  many generations at 33810, [129999-inf) inf-129999="many"
            msprime.PopulationParametersChange(
                time=129999, initial_size=33810, population_id=0
            ),
        ],
    )


_species.add_demographic_model(_Angus_1B16())
