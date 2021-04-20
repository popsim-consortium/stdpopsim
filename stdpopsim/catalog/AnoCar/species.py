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

# No recombination rate yet for this species.
# Author of BourgeoisEtAl is sending the recombination map
# Placeholder rate of 1cM/Mb used
_recombo_rate = 1e-8

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
    citations=[_BourgeoisEtAl],
)
stdpopsim.utils.append_common_synonyms(_genome)

_species = stdpopsim.Species(
    id="AnoCar",
    ensembl_id="anolis_carolinensis",
    name="Anolis carolinensis",
    common_name="Anole lizard",
    genome=_genome,
    generation_time=1.5,
    # they live between 1-2 years after they are able to mate
    # they mature 8 to 9 months after they are born
    # can live up to 8 years in captivity
    population_size=3.05e6,
    # poulation size caculated from theta caculations
    # theta = 4Neu, theta from table 1
    # Ne averaged across the 5 populations from BourgeoisEtAl
    citations=[_LovernEtAl, _BourgeoisEtAl],
)

stdpopsim.register_species(_species)
