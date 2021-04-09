import stdpopsim

from . import genome_data

_BeyeEtAl = stdpopsim.Citation(
    doi="https://dx.doi.org/10.1101%2Fgr.5680406",
    year=2006,
    author="Beye et al.",
    reasons={stdpopsim.CiteReason.REC_RATE},
)
_recombination_rate = {
    "CM009931.2": 23.9e-8,
    "CM009932.2": 24.6e-8,
    "CM009933.2": 24.1e-8,
    "CM009934.2": 27.6e-8,
    "CM009935.2": 21.4e-8,
    "CM009936.2": 21.2e-8,
    "CM009937.2": 23.4e-8,
    "CM009938.2": 20.9e-8,
    "CM009939.2": 24.6e-8,
    "CM009940.2": 24.8e-8,
    "CM009941.2": 20.3e-8,
    "CM009942.2": 21.2e-8,
    "CM009943.2": 23.4e-8,
    "CM009944.2": 24.6e-8,
    "CM009945.2": 22.1e-8,
    "CM009946.2": 22.8e-8,
    "CM009947.2": 0,
}


_YangEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1038/nature14649",
    year=2015,
    author="Yang et al.",
    reasons={stdpopsim.CiteReason.MUT_RATE},
)

_overall_mutation_rate = 9 * 10e-9
_mutation_rate = {
    "CM009931.2": _overall_mutation_rate,
    "CM009932.2": _overall_mutation_rate,
    "CM009933.2": _overall_mutation_rate,
    "CM009934.2": _overall_mutation_rate,
    "CM009935.2": _overall_mutation_rate,
    "CM009936.2": _overall_mutation_rate,
    "CM009937.2": _overall_mutation_rate,
    "CM009938.2": _overall_mutation_rate,
    "CM009939.2": _overall_mutation_rate,
    "CM009940.2": _overall_mutation_rate,
    "CM009941.2": _overall_mutation_rate,
    "CM009942.2": _overall_mutation_rate,
    "CM009943.2": _overall_mutation_rate,
    "CM009944.2": _overall_mutation_rate,
    "CM009945.2": _overall_mutation_rate,
    "CM009946.2": _overall_mutation_rate,
    "CM009947.2": _overall_mutation_rate,
}

_TheHGSC = stdpopsim.Citation(
    doi="http://dx.doi.org/10.1038%2Fnature05260",
    year=2006,
    author="The Honeybee Genome Sequencing Consortium",
    reasons={stdpopsim.CiteReason.ASSEMBLY},
)

_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
    citations=[
        _YangEtAl,
        _BeyeEtAl,
        _TheHGSC,
    ],
)


_NelsonEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1111/mec.14122",
    year=2017,
    author="Nelson et al.",
    reasons={stdpopsim.CiteReason.GEN_TIME},
)

_WallbergEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1038/ng.3077",
    year=2014,
    author="Wallberg et al.",
    reasons={stdpopsim.CiteReason.POP_SIZE},
)

_species = stdpopsim.Species(
    id="ApiMel",
    ensembl_id="apis_mellifera",
    name="Apis mellifera",
    common_name="Apis mellifera (DH4)",
    genome=_genome,
    generation_time=2,
    population_size=2e05,
    citations=[
        _WallbergEtAl,
        _NelsonEtAl,
    ],
)

stdpopsim.register_species(_species)
