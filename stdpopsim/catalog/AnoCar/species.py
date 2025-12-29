import stdpopsim

from . import genome_data

_LovernEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1093/ilar.45.1.54",
    year=2004,
    author="Lovern et al.",
    reasons={stdpopsim.CiteReason.GEN_TIME},
)

_BourgeoisEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1093/gbe/evz110",
    year=2019,
    author="Pombi et al.",
    reasons={
        stdpopsim.CiteReason.POP_SIZE,
        stdpopsim.CiteReason.MUT_RATE,
        stdpopsim.CiteReason.REC_RATE,
    },
)

_YurchenkoEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1093/gbe/evaa161",
    year=2020,
    author="Yurchenko et al.",
    reasons={stdpopsim.CiteReason.REC_RATE},
)

# Recombination rate from Yurchenko et al.
# which was calculated from a linkage map
# from a different lizard species, Zootoca vivipara
# the "common lizard"
# they estimated male and female recombination rates
# at 1.49 cM/Mb and 1.69 cM/Mb, respectively
# we can use the average of the two
_recombo_rate = 1.59e-8

_recombination_rate = {
    "1": _recombo_rate,
    "2": _recombo_rate,
    "3": _recombo_rate,
    "4": _recombo_rate,
    "5": _recombo_rate,
    "6": _recombo_rate,
    "LGa": _recombo_rate,
    "LGb": _recombo_rate,
    "LGc": _recombo_rate,
    "LGd": _recombo_rate,
    "LGf": _recombo_rate,
    "LGg": _recombo_rate,
    "LGh": _recombo_rate,
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
    "LGa": _species_ploidy,
    "LGb": _species_ploidy,
    "LGc": _species_ploidy,
    "LGd": _species_ploidy,
    "LGf": _species_ploidy,
    "LGg": _species_ploidy,
    "LGh": _species_ploidy,
    "MT": 1,
}

# Mutation rate
_overall_rate = 2.1e-10

_mutation_rate = {
    "1": _overall_rate,
    "2": _overall_rate,
    "3": _overall_rate,
    "4": _overall_rate,
    "5": _overall_rate,
    "6": _overall_rate,
    "LGa": _overall_rate,
    "LGb": _overall_rate,
    "LGc": _overall_rate,
    "LGd": _overall_rate,
    "LGf": _overall_rate,
    "LGg": _overall_rate,
    "LGh": _overall_rate,
    "MT": _overall_rate,
}

_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
    ploidy=_ploidy,
    citations=[_BourgeoisEtAl],
)
stdpopsim.utils.append_common_synonyms(_genome)

_species = stdpopsim.Species(
    id="AnoCar",
    ensembl_id="anolis_carolinensis",
    name="Anolis carolinensis",
    common_name="Anole lizard",
    separate_sexes=True,
    genome=_genome,
    generation_time=1.5,
    # they live between 1-2 years after they are able to mate
    # they mature 8 to 9 months after they are born
    # can live up to 8 years in captivity
    population_size=3.05e6,
    # poulation size caculated from theta caculations
    # theta = 4Neu, theta from table 1
    # Ne averaged across the 5 populations from BourgeoisEtAl
    ploidy=_species_ploidy,
    citations=[_LovernEtAl, _BourgeoisEtAl],
)

stdpopsim.register_species(_species)
