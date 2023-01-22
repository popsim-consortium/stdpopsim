import stdpopsim

from . import genome_data

<<<<<<< HEAD
# Chromosome-specific recombination rate
_recombination_rate = {
    "1": 2.75145215530419e-09,
    "2": 3.31510028178352e-09,
    "3": 3.51437699680511e-09,
    "4": 3.56171480378189e-09,
    "5": 3.86074685482261e-09,
    "6": 3.43642611683849e-09,
    "7": 3.15498882608124e-09,
    "8": 4.33723938517729e-09,
    "9": 4.87408610885459e-09,
    "10": 4.28382272762676e-09,
    "11": 5.70613409415121e-09,
    "12": 4.32252701579385e-09,
    "13": 5.07783234828935e-09,
    "14": 4.19354838709677e-09,
    "15": 5.82580833090591e-09,
    "16": 7.61865541027489e-09,
    "17": 6.3034188034188e-09,
    "18": 7.13171383998217e-09,
    "19": 1.45096731154103e-08,
    "X": 2.72149984880556e-09,
    "Y": 0,
    "MT": 0,
}

# Generic and chromosome-specific ploidy
_species_ploidy = 2
_ploidy = {
    "1": _species_ploidy,
    "2": _species_ploidy,
    "3": _species_ploidy,
    "4": _species_ploidy,
    "5": _species_ploidy,
    "6": _species_ploidy,
    "7": _species_ploidy,
    "8": _species_ploidy,
    "9": _species_ploidy,
    "10": _species_ploidy,
    "11": _species_ploidy,
    "12": _species_ploidy,
    "13": _species_ploidy,
    "14": _species_ploidy,
    "15": _species_ploidy,
    "16": _species_ploidy,
    "17": _species_ploidy,
    "18": _species_ploidy,
    "19": _species_ploidy,
    "X": _species_ploidy,
    "Y": 1,
    "MT": 1,
}

# Generic and chromosome-specific mutation rate
_overall_rate = 5.4e-09
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
    "X": _overall_rate,
    "Y": _overall_rate,
    "MT": 3.7e-08,
=======
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
>>>>>>> 52774ae (ran maintenance script to pull ensembl info)
}

_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
<<<<<<< HEAD
    ploidy=_ploidy,
    citations=[
        stdpopsim.Citation(
            author="Cox et al.",
            year=2009,
            doi="https://doi.org/10.1534/genetics.109.105486",
            reasons={stdpopsim.CiteReason.REC_RATE},
        ),
        stdpopsim.Citation(
            author="Uchimura et al.",
            year=2015,
            doi="http://www.genome.org/cgi/doi/10.1101/gr.186148.114",
            reasons={stdpopsim.CiteReason.MUT_RATE},
        ),
        stdpopsim.Citation(
            author="Goios et al.",
            year=2007,
            doi="http://www.genome.org/cgi/doi/10.1101/gr.5941007",
            reasons={stdpopsim.CiteReason.MUT_RATE},
=======
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
>>>>>>> 52774ae (ran maintenance script to pull ensembl info)
        ),
    ],
)

_species = stdpopsim.Species(
    id="MusMus",
    ensembl_id="mus_musculus",
    name="Mus musculus",
    common_name="Mouse",
    genome=_genome,
<<<<<<< HEAD
    generation_time=0.75,
    population_size=500000,
    citations=[
        stdpopsim.Citation(
            author="Fujiwara et al.",
            year=2022,
            doi="https://doi.org/10.1093/gbe/evac068",
            reasons={stdpopsim.CiteReason.GEN_TIME, stdpopsim.CiteReason.POP_SIZE},
        ),
        stdpopsim.Citation(
            author="Phifer-Rixey, Harr et al.",
            year=2020,
            doi="https://doi.org/10.1186/s12862-020-01666-9",
            reasons={stdpopsim.CiteReason.GEN_TIME},
        ),
        stdpopsim.Citation(
            author="Phifer-Rixey, Bonhomme, et al.",
            year=2012,
            doi="https://doi.org/10.1093/molbev/mss105",
            reasons={stdpopsim.CiteReason.POP_SIZE},
=======
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
>>>>>>> 52774ae (ran maintenance script to pull ensembl info)
        ),
    ],
)

stdpopsim.register_species(_species)
