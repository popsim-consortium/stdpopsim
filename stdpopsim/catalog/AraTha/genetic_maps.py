import stdpopsim
import msprime
from .species import _selfing_correction

_species = stdpopsim.get_species("AraTha")

_genetic_map_citation = stdpopsim.Citation(
    doi="https://doi.org/10.1038/hdy.2011.95",
    author="Salomé et al.",
    year=2012,
    reasons={stdpopsim.CiteReason.GEN_MAP},
)


# define a subclass to correct recombination rates for selfing
# The original rate map in "SalomeAveraged_TAIR10" is calculated for
# a set of outcrossing F2 populations, but A. thaliana is highly selfing.
# The selfing correction is derived in Nordborg 2000 Genetics
# and is equal to 1 - σ/(2 - σ), where σ is the selfing rate.
# We use σ = 0.97 from Platt et al. 2010 PLOS Genetics.
class _AraThaGeneticMap(stdpopsim.GeneticMap):
    def get_chromosome_map(self, id):
        rate_map = super().get_chromosome_map(id)
        return msprime.RateMap(
            position=rate_map.position,
            rate=rate_map.rate * _selfing_correction,
        )


# Here, we use `_AraThaGeneticMap` in place of `stdpopsim.GeneticMap`
_gm = _AraThaGeneticMap(
    species=_species,
    id="SalomeAveraged_TAIR10",  # ID for genetic map, see naming conventions
    description="Crossover frequency map averaged over 17 populations",
    long_description="""
        This map is based on the study of crossover frequencies in over 7000
        plants in 17 F2 populations derived from crosses between 18 A. thaliana
        accessions. Salomé et al provide genetic maps for each of these
        populations. To get a single map for each chromosome, the Haldane map
        function distances were converted to recombination rates (cM/Mb) for
        each cross and then averaged across the 17 populations using loess.
        To account for the high selfing rate in A. thaliana, the recombination
        rates were adjusted following Nordborg 2000 Genetics, section
        'Recombination and selfing' (bottom left of page 925).
        A selfing rate of 0.97 was used from Platt et al. 2010 PLOS Genetics.
        The map was constructed on genome version TAIR7, and lifted to TAIR10.
        """,
    url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/"
        "genetic_maps/AraTha/salome2012_liftover_maps.tar.gz"
    ),
    sha256="6eae8d22eff46eca27c288d0431f25b51cf22617be8eeb2fa1aaa0bfe00c0f76",
    file_pattern="salome2012_liftover_maps/TAIR10_lifted_chr{id}.txt",
    citations=[_genetic_map_citation],
)

_species.add_genetic_map(_gm)
