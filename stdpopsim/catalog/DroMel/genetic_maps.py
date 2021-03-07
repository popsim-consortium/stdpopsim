import stdpopsim

_species = stdpopsim.get_species("DroMel")

_gm = stdpopsim.GeneticMap(
    species=_species,
    id="ComeronCrossover_dm6",
    description="Crossover map from meioses products of 8 lab crosses",
    long_description="""
        The crossover map from a study of 8 crosses of 12 highly
        inbred lines of D. melanogaster. This is based on the
        products of 5,860 female meioses from whole genome sequencing data.
        Recombination rates were calculated from the density of individual
        recombination events that were detected in crosses. This map was
        subsequently lifted over to the dm6 assembly.
        """,
    url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/genetic_maps/"
        "DroMel/comeron2012_maps.tar.gz"
    ),
    sha256="08185a0e3b0ad26eefe69fc6bdb8f3f599a760e11e87dd343335b33d1563f62a",
    file_pattern="genetic_map_comeron2012_dm6_chr{id}.txt",
    citations=[
        stdpopsim.Citation(
            author="Comeron et al",
            doi="https://doi.org/10.1371/journal.pgen.1002905",
            year=2012,
            reasons={stdpopsim.CiteReason.GEN_MAP},
        )
    ],
)
_species.add_genetic_map(_gm)

_gm = stdpopsim.GeneticMap(
    species=_species,
    id="ComeronCrossoverV2_dm6",
    description="Crossover map from meioses products of 8 lab crosses",
    long_description="""
        The crossover map from a study of 8 crosses of 12 highly
        inbred lines of D. melanogaster. This is based on the
        products of 5,860 female meioses from whole genome sequencing data.
        Recombination rates were calculated from the density of individual
        recombination events that were detected in crosses. This map was
        subsequently lifted over to the dm6 assembly using the available
        maintenance code command:

        python liftOver_comeron2012.py \
            --winLen 1000 \
            --gapThresh 1000000 \
            --useAdjacentAvg \
            --retainIntermediates
        """,
    url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/genetic_maps/"
        "DroMel/comeron2012v2_maps.tar.gz"
    ),
    sha256="79484e6042c36b06946af672cabf2ec3012b6bfc6e0018517c8e7be4769cd334",
    file_pattern="genetic_map_comeron2012v2_dm6_chr{id}.txt",
    citations=[
        stdpopsim.Citation(
            author="Comeron et al",
            doi="https://doi.org/10.1371/journal.pgen.1002905",
            year=2012,
            reasons={stdpopsim.CiteReason.GEN_MAP},
        )
    ],
)
_species.add_genetic_map(_gm)
