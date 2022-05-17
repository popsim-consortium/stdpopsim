import stdpopsim

_species = stdpopsim.get_species("PonAbe")

_nater2017 = stdpopsim.Citation(
    author="Nater et al.", year=2017, doi="https://doi.org/10.1016/j.cub.2017.09.047"
)

# There are two genetic maps available for Orangutan species: one for Pongo
# abelii (Sumatran orangutan) and one for Pongo pygmaeus (Bornean orangutan).
# Both recombination maps were inferred using LDhat in Nater et al. (2017),
# doi: 10.1016/j.cub.2017.09.047. Both recombination maps are mapped to PonAbe2.
# Recombination maps from Nater et al. were converted from rho/kbp to cM using
# Watterson's estimator of theta to estimate Ne = 41,000 (Sumatra) and
# Ne = 27,000 (Borneo). See supporting information in Nater et al. for details.
# Maps have been lifted over to PonAbe3 (aka Susie_PABv2).


_gm_pa = stdpopsim.GeneticMap(
    species=_species,
    id="NaterPA_PonAbe3",
    description="From Nater et al. (2017) for Pongo abelii",
    long_description="""
        This genetic map is from the Nater et al. (2017) study, inferred using
        LDhat from n=15 whole-genome sequenced Sumatran orangutan individuals.
        See https://doi.org/10.1016/j.cub.2017.09.047 for more details.
        Lifted over from assembly PonAbe2 (as used in Nater et al.) to PonAbe3.
        """,
    url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/genetic_maps/PonAbe/"
        "NaterPA_PonAbe3.tar.gz"
    ),
    sha256="80470e8d0a237f9cc938cac4595b13a6511b6988706a86db905c43e4619dfe37",
    file_pattern="genetic_map_PonAbe3_{id}.txt",
    citations=[_nater2017.because(stdpopsim.CiteReason.GEN_MAP)],
)
_species.add_genetic_map(_gm_pa)

_gm_pp = stdpopsim.GeneticMap(
    species=_species,
    id="NaterPP_PonAbe3",
    description="From Nater et al. (2017) for Pongo pygmaeus",
    long_description="""
        This genetic map is from the Nater et al. (2017) study, inferred using
        LDhat from n=20 whole-genome sequenced Bornean orangutan individuals.
        See https://doi.org/10.1016/j.cub.2017.09.047 for more details.
        Lifted over from assembly PonAbe2 (as used in Nater et al.) to PonAbe3.
        """,
    url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/genetic_maps/PonAbe/"
        "NaterPP_PonAbe3.tar.gz"
    ),
    sha256="ccdaf29c3a7e515e137bf16305da6da3f48f10d70998588335c650b01a622e4f",
    file_pattern="genetic_map_PonAbe3_{id}.txt",
    citations=[_nater2017.because(stdpopsim.CiteReason.GEN_MAP)],
)
_species.add_genetic_map(_gm_pp)
