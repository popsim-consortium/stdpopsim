import stdpopsim


_species = stdpopsim.get_species("PapAnu")

_gm = stdpopsim.GeneticMap(
    species=_species,
    id="Pyrho_PAnubis1_0",
    description="Pyrho inferred genetic map for Papio Anubis",
    long_description="""
        These estimates were obtained from a sample of Papio Anubis
        individuals from the colony housed at the Southwest National
        Primate Research Center (SNPRC).
        """,
    url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/genetic_maps/"
        "PapAnu/papio_anubis_genetic_map.gz"
    ),
    sha256="bf33f975dda8582cacbb7794ba057c2c4f926e9a11dfc9f0978860c33f1fbe9b",
    file_pattern="pyrho_chr{id}_PAnubis1.0.txt",
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
