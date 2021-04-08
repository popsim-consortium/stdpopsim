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


_gm_pa = stdpopsim.GeneticMap(
    species=_species,
    id="NaterPA_PonAbe2",
    description="From Nater et al. (2017) for Pongo abelii",
    long_description="""
        This genetic map is from the Nater et al. (2017) study, inferred using
        LDhat from n=15 whole-genome sequenced Sumatran orangutan individuals.
        See https://doi.org/10.1016/j.cub.2017.09.047 for more details.
        """,
    url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/genetic_maps/PonAbe/"
        "NaterPA_PonAbe2.tar.gz"
    ),
    sha256="33b0162a6d945dd341bd1086d213ad1bc16949c9210d3b49d92692fd7f831ace",
    file_pattern="Nater_et_al_PA_chr{id}_PonAbe2.txt",
    citations=[_nater2017.because(stdpopsim.CiteReason.GEN_MAP)],
)
_species.add_genetic_map(_gm_pa)

_gm_pp = stdpopsim.GeneticMap(
    species=_species,
    id="NaterPP_PonAbe2",
    description="From Nater et al. (2017) for Pongo pygmaeus",
    long_description="""
        This genetic map is from the Nater et al. (2017) study, inferred using
        LDhat from n=20 whole-genome sequenced Bornean orangutan individuals.
        See https://doi.org/10.1016/j.cub.2017.09.047 for more details.
        """,
    url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/genetic_maps/PonPyg/"
        "NaterPP_PonAbe2.tar.gz"
    ),
    sha256="f4858b7efe15abe28b9367e7e9dc16614dc614df9326f9eedbcd943ca075bed7",
    file_pattern="Nater_et_al_PP_chr{id}_PonAbe2.txt",
    citations=[_nater2017.because(stdpopsim.CiteReason.GEN_MAP)],
)
_species.add_genetic_map(_gm_pp)
