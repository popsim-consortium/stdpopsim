import stdpopsim

from . import genome_data

_ChakrabortyEtAl = stdpopsim.Citation(
        doi="https://doi.org/10.1101/gr.263442.120",
        year=2021,
        author="Chakraborty et al.",
        reasons={stdpopsim.CiteReason.ASSEMBLY},
)

_ComeronEtAl = stdpopsim.Citation(
        doi="https://doi.org/10.1371/journal.pgen.1002905",
        year=2012,
        author="Comeron et al.",
        reasons={stdpopsim.CiteReason.REC_RATE},
)

_LegrandEtAl = stdpopsim.Citation(
        doi="https://doi.org/10.1534/genetics.108.092080",
        year=2009,
        author="Legrand et al.",
        reasons={stdpopsim.CiteReason.GEN_TIME,
            stdpopsim.CiteReason.MUT_RATE,
            stdpopsim.CiteReason.POP_SIZE,
            },
)

_recombination_rate = {
        "2L": 2.4125016027908946e-08,
        "2R": 2.2366522822806982e-08,
        "3L": 1.7985794693631893e-08,
        "3R": 1.7165556232922828e-08,
        "X": 2.9151053903465754e-08,
        "4": 2.0085234464525437e-08,
}

_overall_rate = 1.5e-9
_mutation_rate = {
        "2L": _overall_rate,
        "2R": _overall_rate,
        "3L": _overall_rate,
        "3R": _overall_rate,
        "X": _overall_rate,
        "4": _overall_rate
}

_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
    citations=[
        _ChakrabortyEtAl,
        _ComeronEtAl,
        _LegrandEtAl,
    ],
)

_species = stdpopsim.Species(
    id="DroSec",
    ensembl_id="drosophila_sechellia",
    name="Drosophila sechellia",
    common_name="Drosophila sechellia (Rob3c)",
    genome=_genome,
    generation_time=0.05,
    population_size=100000,
    citations=[_LegrandEtAl],
)

stdpopsim.register_species(_species)
