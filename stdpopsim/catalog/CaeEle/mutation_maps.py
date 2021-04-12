import stdpopsim

_species = stdpopsim.get_species("CaeEle")

_gm = stdpopsim.GeneticMap(
    species=_species,
    id="KonradMutationMap",
    description="Mutation map from mutation accomulation lines",
    long_description="""
        The authors conducted a mutation accumulation (MA) experiment
        at varying population sizes in the nematode Caenorhabditis elegans,
        evolving 35 lines in parallel for 409 generations at three population
        sizes (N = 1, 10, and 100 individuals).
        """,
    url=("http://sesame.uoregon.edu/~ateterina/" "konrad_mut_maps.tgz"),
    sha256="1a571531bca5c63ac0f9715954ddf48c1fc17e1ca238e2ac631996561ceaa202",
    file_pattern="C.elegans.Konrad.2019.{id}.hapmap.txt",
    citations=[
        stdpopsim.Citation(
            author="Konrad et al.",
            doi="https://doi.org/10.1534/genetics.119.302054",
            year=2019,
            reasons={stdpopsim.CiteReason.MUT_MAP},
        ),
        stdpopsim.Citation(
            author="Konrad et al.",
            doi="https://doi.org/10.1093/molbev/msx051",
            year=2017,
            reasons={stdpopsim.CiteReason.MUT_MAP},
        ),
    ],
)
_species.add_mutation_map(_gm)
