import stdpopsim

from . import genome_data

_recombination_rate = {
    "1": 3.17e-8,
    "2": 5.61e-8,
    "3": 5.10e-8,
    "4": 4.97e-8,
    "5": 5.15e-8,
    "6": 3.40e-8,
    "7": 3.76e-8,
    "8": 5.28e-8,
    "9": 5.31e-8,
    "10": 3.16e-8,
    "11": 4.47e-8,
    "12": 3.13e-8,
    "13": 3.08e-8,
    "14": 5.47e-8,
    "15": 4.78e-8,
    "16": 4.71e-8,
    "17": 3.94e-8,
    "18": 3.16e-8,
    "19": 3.11e-8,
    "20": 3.45e-8,
    "21": 3.71e-8,
}

_overall_rate = 2.9e-9
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
}

_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
    citations=[
        stdpopsim.Citation(
            author="Davey et al",
            year=2017,
            doi="https://doi.org/10.1002/evl3.12",
            reasons={stdpopsim.CiteReason.ASSEMBLY, stdpopsim.CiteReason.REC_RATE},
        ),
        stdpopsim.Citation(
            author="Keightley et al",
            year=2015,
            doi="https://doi.org/10.1093/molbev/msu302",
            reasons={stdpopsim.CiteReason.MUT_RATE},
        ),
    ],
)
stdpopsim.utils.append_common_synonyms(_genome)

_species = stdpopsim.Species(
    id="HelMel",
    ensembl_id="heliconius_melpomene",
    name="Heliconius melpomene",
    common_name="Heliconius melpomene",
    genome=_genome,
    generation_time=1 / 10,
    population_size=2.1e6,
    citations=[
        stdpopsim.Citation(
            author="Pardo-Diaz et al",
            year=2012,
            doi="https://doi.org/10.1371/journal.pgen.1002752",
            reasons={stdpopsim.CiteReason.POP_SIZE, stdpopsim.CiteReason.GEN_TIME},
        ),
    ],
)

stdpopsim.register_species(_species)
