import stdpopsim

from . import genome_data

# These are in Table 1 of Juneja et al:
_recombination_rate = {"1": 1.97036e-6, "2": 1.51361e-6, "3": 1.886e-6,
                       "4": 2.22347e-6,
                        "5":2.54212e-6,
                         "6": 1.9829e-6,
                          "7": 2.61388e-6,
                           "8": 2.11591e-6,
                           "MT": 0,
                           "Pt": 0,
                           "Test":1.97036e-8}

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
    "MT": 1,
    "Pt":1,
    "Test":_species_ploidy
}


_VelascoEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1534/g3.116.032672",
    year=2016,
    author="Velasco et al.",
    reasons={
        stdpopsim.CiteReason.GEN_TIME,
        stdpopsim.CiteReason.POP_SIZE,
        stdpopsim.CiteReason.MUT_RATE,
    },
)
_CastaneraEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1093/hr/uhae106",
    year=2024,
    author="Castanera et al.",
    reasons={
        stdpopsim.CiteReason.ASSEMBLY
    },
)
# _MasGomezEtAl = stdpopsim.Citation(
#     doi="https://doi.org/10.3389/fpls.2023.1165847", #Real doi is pending
#     year=2023,
#     author="Mas-Gomez et al.",
#     reasons={
#         stdpopsim.CiteReason.REC_RATE
#     }),
_DAmico_WillmanEtal = stdpopsim.Citation(
    doi="https://doi.org/10.1093/g3journal/jkac065",
    year=2022,
    author="D'Amico-Willman et al.",
    reasons={
        stdpopsim.CiteReason.ASSEMBLY #mitochondrial
    },
)
_WangEtal = stdpopsim.Citation(
    doi="https://doi.org/10.1038/s41598-020-67264-3",
    year=2020,
    author="Wang et al.",
    reasons={
        stdpopsim.CiteReason.ASSEMBLY #chloroplast
    },
)


_overall_rate = 10e-8  # per generation, from Velasco et al. 2016
_mutation_rate = {
    "1": _overall_rate,
    "2": _overall_rate,
    "3": _overall_rate,
    "4": _overall_rate,
    "5": _overall_rate,
    "6": _overall_rate,
    "7": _overall_rate,
    "8": _overall_rate,
    "MT": _overall_rate,
    "Pt":_overall_rate,
    "Test":_overall_rate
}



_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
    ploidy=_ploidy,
    citations=[
_VelascoEtAl, _CastaneraEtAl
    ],
)
stdpopsim.utils.append_common_synonyms(_genome)

_species = stdpopsim.Species(
    id="PruDul",
    ensembl_id="prunus_dulcis",
    name="Prunus dulcis",
    common_name="Prunus dulcis",
    genome=_genome,
    generation_time=10,
    ploidy=_species_ploidy,
    # the estimated population size from the plot of Velasco et al. 2016
    population_size=200,
    citations=[_VelascoEtAl, _CastaneraEtAl,_DAmico_WillmanEtal, _WangEtal], #_MasGomezEtAl,
)

stdpopsim.register_species(_species)
