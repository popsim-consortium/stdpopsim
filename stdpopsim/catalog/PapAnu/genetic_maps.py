import stdpopsim


_species = stdpopsim.get_species("PapAnu")

_gm = stdpopsim.GeneticMap(
    species=_species,
    id="Pyrho_PAnubis1.0",
    description="Pyrho inferred genetic map for Papio Anubis",
    long_description="""
        These estimates were obtained from a sample Papio Anubis
        individuals from the colony housed at the Southwest National
        Primate Research Center (SNPRC).
        """,
    url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/genetic_maps/"
        "PapAnu/papio_anubis_genetic_map.tar.gz"
    ),
    sha256="",
    file_pattern="Pyrho_PAnubis1.0_chr{id}.txt",
    citations=[
        stdpopsim.Citation(
            year=2022,
            author="Wall et. al.",
            doi="https://doi.org/10.1093/gbe/evac040",
            reasons={stdpopsim.CiteReason.GEN_MAP},
        ),
    ],
)
_species.add_genetic_map(_gm)
