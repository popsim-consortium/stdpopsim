import stdpopsim

_species = stdpopsim.get_species("ZeaMay")

_gm = stdpopsim.GeneticMap(
    species=_species,
    id="FineMap_B73v5",
    description="Composite crossover map calibrated to the Ogut NAM map",
    long_description="""
        Piecewise-constant crossover map on Zm-B73-REFERENCE-NAM-5.0
        coordinates (the finemap_v5.bed map of FineMap). About 402,000
        crossover intervals from four published datasets (US and Chinese NAM
        RILs from Rodgers-Melnick et al. 2015; European populations from
        Bauer et al. 2013; a landrace and a teosinte population from Samayoa
        et al. 2021) were lifted over to B73 v5. Each crossover was spread
        uniformly over its interval and the densities summed. Each
        chromosome's profile was then scaled so that its total length equals
        the chromosome length of the Ogut et al. (2015) NAM consensus map.
        Rates within a chromosome are relative rates; absolute totals come
        from the Ogut map. Resolution is limited by marker spacing in the
        source data, so kb-scale peaks should not be read as hotspots.
        Regions beyond the first and last mapped crossover interval on each
        chromosome are given a rate of zero.
        """,
    url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/"
        "genetic_maps/ZeaMay/FineMap_B73v5.tar.gz"
    ),
    sha256="5bd69f16bdc48da60195ed21ba8a97195b24d5c2fe746601e67759cc89af9117",
    file_pattern="FineMap_B73v5_chr{id}.txt",
    citations=[
        stdpopsim.Citation(
            author="Ross-Ibarra",
            year=2026,
            doi="https://doi.org/10.5281/zenodo.19639077",
            reasons={stdpopsim.CiteReason.GEN_MAP},
        ),
        stdpopsim.Citation(
            author="Ogut et al.",
            year=2015,
            doi="https://doi.org/10.1038/hdy.2014.123",
            reasons={stdpopsim.CiteReason.GEN_MAP},
        ),
    ],
)
_species.add_genetic_map(_gm)
