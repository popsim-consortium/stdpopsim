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
    "A01": -1,
    "A02": -1,
    "A03": -1,
    "A04": -1,
    "A05": -1,
    "A06": -1,
    "A07": -1,
    "A08": -1,
    "A09": -1,
    "A10": -1,
}

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

_mutation_rate = {
    "A01": -1,
    "A02": -1,
    "A03": -1,
    "A04": -1,
    "A05": -1,
    "A06": -1,
    "A07": -1,
    "A08": -1,
    "A09": -1,
    "A10": -1,
}

_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
    # [ Implementers: please insert citations for the papers you are basing
    # the estimates for recombination and mutation rates. The assembly
    # citation is optional and can be deleted if not needed.]
    citations=[
        stdpopsim.Citation(
            author="", year=-1, doi="", reasons={stdpopsim.CiteReason.ASSEMBLY}
        ),
        stdpopsim.Citation(
            author="", year=-1, doi="", reasons={stdpopsim.CiteReason.REC_RATE}
        ),
        stdpopsim.Citation(
            author="", year=-1, doi="", reasons={stdpopsim.CiteReason.MUT_RATE}
        ),
    ],
)

_species = stdpopsim.Species(
    id="BraRap",
    ensembl_id="brassica_rapa",
    name="Brassica rapa",
    common_name="Brassica rapa",
    genome=_genome,
    # [Implementers: you must provide an estimate of the generation_time.
    # Please also add a citation for this.]
    generation_time=0,
    # [Implementers: you must provide an estimate of the population size.
    # TODO: give a definition of what this should be.
    # Please also add a citation for this below..]
    population_size=0,
    citations=[
        stdpopsim.Citation(
            author="", year=-1, doi="", reasons={stdpopsim.CiteReason.POP_SIZE}
        ),
        stdpopsim.Citation(
            author="", year=-1, doi="", reasons={stdpopsim.CiteReason.GEN_TIME}
        ),
    ],
)

stdpopsim.register_species(_species)
