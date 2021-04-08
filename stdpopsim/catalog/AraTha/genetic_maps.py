import stdpopsim

_species = stdpopsim.get_species("AraTha")


_gm = stdpopsim.GeneticMap(
    species=_species,
    id="SalomeAveraged_TAIR7",
    description="Crossover frequency map averaged over 17 populations",
    long_description="""
        This map is based on the study of crossover frequencies in over 7000
        plants in 17 F2 populations derived from crosses between 18 A. thaliana
        accessions. Salomé et al provide genetic maps for each of these
        populations. To get a single map for each chromosome, the Haldane map
        function distances were converted to recombination rates (cM/Mb) for
        each cross and then averaged across the 17 populations using loess.
        """,
    url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/genetic_maps/"
        "AraTha/salome2012_maps.tar.gz"
    ),
    sha256="49745e1cab87d59e33eacfdf66303839632d3b07883dd55a99fe1dc27b336ac6",
    file_pattern="arab_chr{id}_map_loess.txt",
    citations=[
        stdpopsim.Citation(
            doi="https://doi.org/10.1038/hdy.2011.95",
            author="Salomé et al.",
            year=2012,
            reasons={stdpopsim.CiteReason.GEN_MAP},
        )
    ],
)
_species.add_genetic_map(_gm)
