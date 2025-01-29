import stdpopsim

from . import genome_data

_ZhangEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1016/j.gpb.2022.02.001",
    year=2022,
    author="Zhang et al.",
    reasons={
        stdpopsim.CiteReason.MUT_RATE,
        stdpopsim.CiteReason.GEN_TIME,
        stdpopsim.CiteReason.POP_SIZE,
    },
)
_JohnssonEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1186/s12711-021-00643-0",
    year=2021,
    author="Johnsson et al.",
    reasons={stdpopsim.CiteReason.REC_RATE},
)
_WarrEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1093/gigascience/giaa051",
    year=2020,
    author="Warr et al.",
    reasons={stdpopsim.CiteReason.ASSEMBLY},
)
# This is the per-chromosome recombination rate, typically the mean
# rate along the chromosome.
# Values in this dictionary are set to -1 by default, so you have
# to update each one. These should be derived from the most reliable
# data and how they were obtained should be documented here.
# The appropriate citation must be added to the list of
# recombination_rate_citations in the Genome object.
# X recombination rate is mean of autosomes because it is not in paper.
_recombination_rate = {
    "1": 5.335247608983185e-09,
    "2": 8.6057260999254e-09,
    "3": 9.751141303348547e-09,
    "4": 9.79460049114614e-09,
    "5": 1.1992632296111e-08,
    "6": 8.786173521918499e-09,
    "7": 1.0887772285689678e-08,
    "8": 8.923031414594611e-09,
    "9": 9.337506174170301e-09,
    "10": 1.654233205589542e-08,
    "11": 1.2145120574910172e-08,
    "12": 1.6905062453467094e-08,
    "13": 6.162796630328349e-09,
    "14": 8.721569375548237e-09,
    "15": 8.114468762954893e-09,
    "16": 1.1407464686635653e-08,
    "17": 1.3779949516098884e-08,
    "18": 1.3679923888648167e-08,
    "X": 9.41539636725104e-09,
    "Y": 0,
    "MT": 0,
}

# the de novo mutation rate is inferred by Zhang et al 2022,
#  which was used for the PSMC demographic inference
_overall_rate = 3.6e-9
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
    "X": _overall_rate,
    "Y": _overall_rate,
    "MT": _overall_rate,
}

# This is the per-chromosome ploidy.

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
    "X": _species_ploidy,
    "Y": 1,
    "MT": 1,
}

_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
    ploidy=_ploidy,
    citations=[
        _WarrEtAl,
        _JohnssonEtAl,
        _ZhangEtAl,
    ],
)

_species = stdpopsim.Species(
    id="SusScr",
    ensembl_id="sus_scrofa",
    name="Sus scrofa",
    common_name="Pig",
    genome=_genome,
    ploidy=2,
    # "Pigs are typically social animals, breeding in the form of polygamy [25]. 
    # The ages of estrus in sows and boars in the wild are different. 
    # The age at first pregnancy varies in the wild from about 10 to 20 months [26],
    # while boars begin rut when they are three to five years old. 
    # The first rut age of 4–5 years was documented in Russian wild boars [27], 
    # and 3–4 years was recorded in Chinese wild boars [28]. 
    # Pigs are multiparous animals, and the gestation period lasts about 114–130 days [29]. 
    # Comparing to the animals with single birth, such as cattle and yaks, 
    # we believe that the generation transmission of pigs is of good continuity. 
    # Therefore, we set the generation interval of pigs as 3 years, 
    # which is roughly equivalent to the average age of the first pregnancy in sows and the beginning of rut in boars 
    # plus the pregnant gestation period of sows." (Zhang et al., p. 1042 Genomics Proteomics Bioinformatics).
    generation_time=3,
    # [Implementers: you must provide an estimate of the population size.
    # TODO: give a definition of what this should be.
    # Please also add a citation for this below..]
    population_size=0,
    citations=[
        stdpopsim.Citation(
            author="", year=-1, doi="", reasons={stdpopsim.CiteReason.POP_SIZE}
        ),
        _ZhangEtAl,
    ],
)

stdpopsim.register_species(_species)
