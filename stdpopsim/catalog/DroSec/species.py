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
    reasons={
        stdpopsim.CiteReason.GEN_TIME,
        stdpopsim.CiteReason.MUT_RATE,
        stdpopsim.CiteReason.POP_SIZE,
    },
)

# Mean chromosomal rates, calculated by taking the
# average of all rates in the
# "ComeronCrossoverV2_dm6" genetic map, weighted by
# their respective interval lengths.
# Chromosome 4 isn't in this map, so the average of
# 2L, 2R, 3L and 3R weighted by their respective
# chromosome lengths was used instead.
_recombination_rate = {
    "2L": 2.40462600791e-08,
    "2R": 2.23458641776e-08,
    "3L": 1.79660308862e-08,
    "3R": 1.71642045777e-08,
    "4": 2.00579550709e-08,
    "X": 2.89650687913e-08,
}

# Mutation rate is set to that used by
# Legrand et al. in an ABC selection of
# demographic scenarios (page 1200).
# They state this is an estimate taken
# from D. simulans (Wall et al. 2002).
_overall_rate = 1.5e-9
_mutation_rate = {
    "2L": _overall_rate,
    "2R": _overall_rate,
    "3L": _overall_rate,
    "3R": _overall_rate,
    "X": _overall_rate,
    "4": _overall_rate,
}

# We could not auto-pull the genome data from ensemble
# so instead we used the most up-to-date assembly
# currently available from NCBI.
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
stdpopsim.utils.append_common_synonyms(_genome)

# Generation time was set to that used by
# by Legrand et al. in an ABC selection of demographic
# scenarios (page 1200).
# Population size was estimated in the same paper (page 1202).
_species = stdpopsim.Species(
    id="DroSec",
    ensembl_id="drosophila_sechellia",
    name="Drosophila sechellia",
    common_name="Drosophila sechellia",
    genome=_genome,
    generation_time=0.05,
    population_size=100000,
    citations=[_LegrandEtAl],
)

stdpopsim.register_species(_species)
