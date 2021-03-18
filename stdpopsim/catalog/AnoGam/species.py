import stdpopsim

from . import genome_data

# [The following are notes for implementers and should be deleted
#  once the recombination rates have been inserted]
# This is the per-chromosome recombination rate, typically the mean
# rate along the chromosome.
# Values in this dictionary are set to -1 by default, so you have
# to update each one. These should be derived from the most reliable
# data and how they were obtained should be documented here.
# The appropriate citation must be added to the list of
# recombination_rate_citations in the Genome object.

_recombination_rate = {
    "2L": 0,  # setting to zero because of inversion
    "2R": 1.3e-8,
    "3L": 1.6e-8,
    "3R": 1.3e-8,
    "X": 1e-8,
    "Mt": 0,
}

_PombiEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.4269/ajtmh.2006.75.901",
    year=2006,
    author="Pombi et al.",
    reasons={stdpopsim.CiteReason.REC_RATE},
)

# [The following are notes for implementers and should be deleted
#  once the mutation rates have been inserted]
# This is the per-chromosome mutation rate, typically the mean
# rate along the chromosome. If per chromosome rates are not available,
# the same value should be used for each chromosome. In this case,
# please use a variable to store this value, rather than repeating
# the same numerical constant, e.g.
# _mutation_rate = {
#    1: _overall_rate,
#    2: _overall_rate,
#    ...
# Values in this dictionary are set to -1 by default, so you have
# to update each one. These should be derived from the most reliable
# data and how they were obtained should be documented here.
# The appropriate citation must be added to the list of
# mutation_rate_citations in the Genome object.

_Ag1000G = stdpopsim.Citation(
    doi="https://doi.org/10.1038/nature24995",
    year=2017,
    author="Ag1000G Consortium",
    reasons={
        stdpopsim.CiteReason.MUT_RATE,
        stdpopsim.CiteReason.GEN_TIME,
        stdpopsim.CiteReason.POP_SIZE,
    },
)

_overall_rate = 5.49e-9
_mutation_rate = {
    "2L": _overall_rate,
    "2R": _overall_rate,
    "3L": _overall_rate,
    "3R": _overall_rate,
    "X": _overall_rate,
    "Mt": _overall_rate,
}

_SharakhovaEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1186/gb-2007-8-1-r5",
    year=2006,
    author="Sharakhova et al.",
    reasons={stdpopsim.CiteReason.ASSEMBLY},
)

_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
    citations=[
        _SharakhovaEtAl,
        _Ag1000G,
        _PombiEtAl,
    ],
)

_species = stdpopsim.Species(
    id="AnoGam",
    ensembl_id="anopheles_gambiae",
    name="Anopheles gambiae",
    common_name="Anopheles gambiae",
    genome=_genome,
    generation_time=1 / 11,
    population_size=6e6,  # Ghana population
    citations=[_Ag1000G],
)

stdpopsim.register_species(_species)
