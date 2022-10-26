import stdpopsim

_species = stdpopsim.get_species("ApiMel")

_genetic_map_citation = stdpopsim.Citation(
    doi="https://doi.org/10.1186/s13059-014-0566-0",
    author="Liu et al.",
    year=2015,
    reasons={stdpopsim.CiteReason.GEN_MAP},
)

_gm = stdpopsim.GeneticMap(
    species=_species,
    id="Liu_Amel_HAv3.1",  # ID for genetic map, see naming conventions
    description="Crossover events mapped in 55 individuals",
    long_description="""
        This map is based on the study of crossover events in 55 individuals
        from three colonies of Apis mellifera ligustica from China. The study
        included three diploid queens, six diploid workers, and 46 haploid
        drones. The map was first lifted from the Amel_4.5 to Amel_HAv3.1
        assembly. To convert the crossover events to genetic position, the
        recombination rate was assumed stable between crossovers and the
        genetic position was inferred from the physical positions of SNPs
        within these crossover ranges.
        """,
    url=(
        "https://stdpopsim.s3.us-west-2.amazonaws.com/"
        "genetic_maps/ApiMel/Liu2015_litfover_maps.tar.gz"
    ),
    sha256="551a1819fb007573b6fa00a088493964dacd089a5d5863cf03b7b73eac0bd405",
    file_pattern="Liu2015_litfover_maps/Amel_HAv3.1_lifted_chr{id}.txt",
    citations=[_genetic_map_citation],
)

_species.add_genetic_map(_gm)
