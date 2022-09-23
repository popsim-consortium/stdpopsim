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

# Legrand et al. used the mutation rate from this paper.
_WallEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1093/genetics/162.1.203",
    year=2002,
    author="Wall et al.",
    reasons={
        stdpopsim.CiteReason.MUT_RATE,
    },
)

_LegrandEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1534/genetics.108.092080",
    year=2009,
    author="Legrand et al.",
    reasons={
        stdpopsim.CiteReason.GEN_TIME,
        stdpopsim.CiteReason.POP_SIZE,
    },
)

# Mean recombination rates, calculated by averaging across
# the crosses the values of each chromosome's genetic map lengths.
# These values were then divided by
# the chromosome size in Mbp (values represent cM/Mbp).
# The recombinatio rate of chromosome 4 was set to 0 as in
# https://petrov.stanford.edu/cgi-bin/recombination-rates_updateR5.pl
_recombination_rate = {
    "2L": 2.28e-8,
    "2R": 2.51e-8,
    "3L": 1.88e-8,
    "3R": 1.91e-8,
    "4": 0,
    "X": 2.85e-8,
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

# Generic and chromosome-specific ploidy
_species_ploidy = 2
_ploidy = {
    "2L": _species_ploidy,
    "2R": _species_ploidy,
    "3L": _species_ploidy,
    "3R": _species_ploidy,
    "X": _species_ploidy,
    "4": _species_ploidy,
}

# We could not auto-pull the genome data from ensemble
# so instead we used the most up-to-date assembly
# currently available from NCBI.
_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
    ploidy=_ploidy,
    citations=[_ChakrabortyEtAl, _ComeronEtAl, _WallEtAl],
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
    ploidy=_species_ploidy,
    citations=[_LegrandEtAl],
)

stdpopsim.register_species(_species)
