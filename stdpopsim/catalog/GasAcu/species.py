import stdpopsim

from . import genome_data

_recombination_rate = {
    "1": 3.11e-08,
    "2": 3.11e-08,
    "3": 3.11e-08,
    "4": 3.11e-08,
    "5": 3.11e-08,
    "6": 3.11e-08,
    "7": 3.11e-08,
    "8": 3.11e-08,
    "9": 3.11e-08,
    "10": 3.11e-08,
    "11": 3.11e-08,
    "12": 3.11e-08,
    "13": 3.11e-08,
    "14": 3.11e-08,
    "15": 3.11e-08,
    "16": 3.11e-08,
    "17": 3.11e-08,
    "18": 3.11e-08,
    "19": 3.11e-08,
    "20": 3.11e-08,
    "21": 3.11e-08,
    "Y": 0,
    "MT": 0,
}

_overall_rate = 3.7e-8
_mutation_rate = {
    "1": _overall_rate,
    "2": _overall_rate,
    "3": _overall_rate,
    "4": _overall_rate,
    "5": _overall_rate,
    "6": _overall_rate,
    "7": _overall_rate,
    "8": _overall_rate,
    "9": _overall_rate,
    "10": _overall_rate,
    "11": _overall_rate,
    "12": _overall_rate,
    "13": _overall_rate,
    "14": _overall_rate,
    "15": _overall_rate,
    "16": _overall_rate,
    "17": _overall_rate,
    "18": _overall_rate,
    "19": _overall_rate,
    "20": _overall_rate,
    "21": _overall_rate,
    "Y": _overall_rate,
    "MT": _overall_rate,
}

_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
    citations=[
        stdpopsim.Citation(
            author="Peichel et al.",
            year=2017,
            doi="https://doi.org/10.1093/jhered/esx058",
            reasons={stdpopsim.CiteReason.ASSEMBLY},
        ),
        stdpopsim.Citation(
            author="Roesti et al.",
            year=2013,
            doi="https://10.1111/mec.12322",
            reasons={stdpopsim.CiteReason.REC_RATE},
        ),
        stdpopsim.Citation(
            author="Liu et al.",
            year=2016,
            doi="https://10.1111/mec.13827",
            reasons={stdpopsim.CiteReason.MUT_RATE},
        ),
    ],
)

_species = stdpopsim.Species(
    id="GasAcu",
    ensembl_id="9307941",
    name="Gasterosteus aculeatus",
    common_name="Three-spined stickleback",
    genome=_genome,
    generation_time=1,
    population_size=1e4,
    citations=[
        stdpopsim.Citation(
            author="Liu et al.",
            year=2016,
            doi="https://10.1111/mec.13827",
            reasons={stdpopsim.CiteReason.POP_SIZE, stdpopsim.CiteReason.GEN_TIME},
        ),
    ],
)

stdpopsim.register_species(_species)
