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


# define a new subclass to handle the selfing rate correction
# The new class inherits all attributes from the original class
# with the exception of the `get_chromosome_map` method
# which is redefined to account for selfing.
# Below, the `get_chromosome_map` method is redefined for the new class.
# For a given chromosome `id` we retrieve the original map with the superclass method.
# Then, we return a new map with the same position array but with the rate array
# multiplied by the selfing correction factor.
# Details of the selfing correction from Nordborg 2000 Genetics
# are outlined in `species.py`.
class _AraThaGeneticMap(stdpopsim.GeneticMap):
    def get_chromosome_map(self, id):  # redefine this method for the new class
        rate_map = super().get_chromosome_map(id)  # retrieve the original map
        return msprime.RateMap(
            position=rate_map.position,  # copy the position array
            rate=rate_map.rate * _selfing_correction,  # apply selfing correction
        )


# Here, we use the `_AraThaGeneticMap` subclass in place of `stdpopsim.GeneticMap`
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
