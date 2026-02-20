import stdpopsim

_species = stdpopsim.get_species("MusMus")

_genetic_map_citation = stdpopsim.Citation(
    doi="https://doi.org/10.1534/genetics.109.105486",
    author="Cox et al.",
    year=2009,
    reasons={stdpopsim.CiteReason.GEN_MAP},
)

_gm = stdpopsim.GeneticMap(
    species=_species,
    id="Cox_etal_2009_v3_GRCm39",  # ID for genetic map, see naming conventions
    description="Revised Cox genetic map, for mouse genome build 39 coordinates",
    long_description="""
        This map results from the Cox et al. (2009) genetic map which has been
        lifted over to the GRCm39 genome which was released in 2020. As described in
        the abstract of Cox et al. (2009), the genetic map resulted from a large
        heterogeneous mouse population and incorporates 10,195 SNPs using a set
        of 47 families comprising 3546 meioses. The methods by which this original
        map was lifted to the GRCm39 genome are described on a github repo:
        https://github.com/kbroman/CoxMapV3
        """,
    url=("https://us-west-2.console.aws.amazon.com/s3/object/stdpopsim/genetic_maps/MusMus/cox.tar.gz"),
    sha256="d589a7adf334f31343aff722aeebe96b1663cab573856d80a883f07bd3cb9b11",
    file_pattern="Cox_etal_2009_v3_GRCm39_{id}.txt",
    citations=[_genetic_map_citation],
)

_species.add_genetic_map(_gm)
