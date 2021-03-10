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
    "1": -1,
    "2": -1,
    "3": -1,
    "4": -1,
    "5": -1,
    "6": -1,
    "7": -1,
    "8": -1,
    "9": -1,
    "10": -1,
    "11": -1,
    "12": -1,
    "13": -1,
    "14": -1,
    "15": -1,
    "16": -1,
    "17": -1,
    "18": -1,
    "19": -1,
    "X": -1,
    "Y": -1,
    "MT": -1,
}

# [The following are notes for implementers and should be deleted
#  once the recombination rates have been inserted]
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
    "1": -1,
    "2": -1,
    "3": -1,
    "4": -1,
    "5": -1,
    "6": -1,
    "7": -1,
    "8": -1,
    "9": -1,
    "10": -1,
    "11": -1,
    "12": -1,
    "13": -1,
    "14": -1,
    "15": -1,
    "16": -1,
    "17": -1,
    "18": -1,
    "19": -1,
    "X": -1,
    "Y": -1,
    "MT": -1,
}

_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
)

# [Implementers: you must add citations for the values that are
# provided. Do this like:
# _genome.recombination_rate_citations.append(
#    stdpopsim.Citation(author=x, date=y, doi=z))
# _genome.mutation_rate_citations.append(
#    stdpopsim.Citation(author=x, date=y, doi=z))

_species = stdpopsim.Species(
    id="MusMus",
    ensembl_id="mus_musculus",
    name="Mus musculus",
    common_name="Mouse",
    genome=_genome,
    # [Implementers: you must provide an estimate of the generation_time.
    # Please also add a citation for this.]
    generation_time=0,
    # [Implementers: you must provide an estimate of the population size.
    # TODO: give a definition of what this should be.
    # Please also add a citation for this.]
    population_size=0,
    population_size_citations=[],
    generation_time_citations=[],
)

stdpopsim.register_species(_species)
