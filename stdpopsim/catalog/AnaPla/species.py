import stdpopsim

from . import genome_data

_HuangEtAl2006 = stdpopsim.Citation(
    author="Huang et al.",
    year=2006,
    doi="https://dx.doi.org/10.1534%2Fgenetics.105.053256",
    reasons={stdpopsim.CiteReason.REC_RATE},
)

_LavretskyEtAl2019 = stdpopsim.Citation(
    author="Lavretsky et al.",
    year=2019,
    doi="https://doi.org/10.1002/ece3.4981",
    reasons={stdpopsim.CiteReason.REC_RATE},
)

_LavretskyEtAl2020 = stdpopsim.Citation(
    author="Lavretsky et al.",
    year=2020,
    doi="https://doi.org/10.1111/mec.15343",
    reasons={stdpopsim.CiteReason.GEN_TIME, stdpopsim.CiteReason.MUT_RATE},
)

_GuoEtAl2020 = stdpopsim.Citation(
    author="Guo et al.",
    year=2021,
    doi="https://doi.org/10.24272/j.issn.2095-8137.2020.133",
    reasons={stdpopsim.CiteReason.POP_SIZE},
)

# link to ENSEMBL project, before publication
# Not used now
_MallardAssembly_unpub = stdpopsim.Citation(
    author="Ensemble project: PRJNA554956",
    year=2020,
    doi="https://www.ebi.ac.uk/ena/browser/view/PRJNA554956",
    reasons={stdpopsim.CiteReason.ASSEMBLY},
)

_ZhuEtAl2021 = stdpopsim.Citation(
    author="Zhu et al.",
    year=2021,
    doi="https://doi.org/10.1038/s41467-021-26272-1",
    reasons={stdpopsim.CiteReason.ASSEMBLY},
)

# def _print_recomb_rates():
#     """
#     Code to produce the following rates, which extends the approach taken by
#     Lavretsky et al 2019 in estimating average recombination rates for
#     chromosomes 1-7. Following this approach, we took the total map lengths
#     of chromosomes 1-18, computed using microsatellites by Huang et al 2019
#     and reported in their Table 1. We divide this rate (in morgans) by the
#     physical length of chromosome in the genome build to obtain an average
#     rate for that chromosome. A genome-wide average rate is computed by
#     taking the ratio of the sum of map lengths and the sum of physical lengths
#     of chromosomes 1-18. This average rate is assumed for all other chromosomes
#     """
#     from genome_data import data as mallard
#
#     chroms = mallard["chromosomes"]
#     # from Huang et al 2006, Table 1.
#     # See: https://academic.oup.com/view-large/280279985
#     huang_chroms = {
#          "1": {"M": 3.170},
#          "2": {"M": 2.260},
#          "3": {"M": 1.123},
#          "4": {"M": 0.934},
#          "5": {"M": 0.786},
#          "6": {"M": 1.201},
#          "7": {"M": 0.981},
#          "8": {"M": 0.427},
#          "9": {"M": 0.359},
#         "10": {"M": 0.574},
#         "11": {"M": 0.333},
#         "12": {"M": 0.703},
#         "13": {"M": 0.307},
#         "14": {"M": 0.175},
#         "15": {"M": 0.069},
#         "16": {"M": 0.062},
#         "17": {"M": 0.048},
#         "18": {"M": 0.021},
#     }
#     for c in huang_chroms:
#         h = huang_chroms[c]
#         h["bp"] = mallard["chromosomes"][c]["length"]
#         h["rate"] = h["M"] / h["bp"]
#     total_M = sum([x["M"] for x in huang_chroms.values()])
#     total_bp = sum([x["bp"] for x in huang_chroms.values()])
#     mean_rate = total_M / total_bp
#     print(f"_total_M = {total_M}")
#     print(f"_total_bp = {total_bp}")
#     print(f"_default_recombination_rate = {mean_rate}")
#     print("_recombination_rate = {")
#     for c in chroms:
#         if c in huang_chroms:
#             print(f"\"{c}\": {huang_chroms[c]['rate']:.2e},")
#         else:
#             print(f'"{c}": _default_recombination_rate,')
#     print("}")

_total_M = 13.533
_total_bp = 943487373
_default_recombination_rate = 1.43e-08
_recombination_rate = {
    "1": 1.52e-08,
    "2": 1.39e-08,
    "3": 9.38e-09,
    "4": 1.20e-08,
    "5": 1.21e-08,
    "6": 3.04e-08,
    "7": 2.59e-08,
    "8": 1.28e-08,
    "9": 1.34e-08,
    "10": 2.50e-08,
    "11": 1.50e-08,
    "12": 3.15e-08,
    "13": 1.41e-08,
    "14": 8.61e-09,
    "15": 3.79e-09,
    "16": 3.86e-09,
    "17": 3.13e-09,
    "18": 1.58e-09,
    "19": _default_recombination_rate,
    "20": _default_recombination_rate,
    "21": _default_recombination_rate,
    "22": _default_recombination_rate,
    "23": _default_recombination_rate,
    "24": _default_recombination_rate,
    "25": _default_recombination_rate,
    "26": _default_recombination_rate,
    "27": _default_recombination_rate,
    "28": _default_recombination_rate,
    "29": _default_recombination_rate,
    "30": _default_recombination_rate,
    "Z": _default_recombination_rate,
    "31": _default_recombination_rate,
    "32": _default_recombination_rate,
    "33": _default_recombination_rate,
    "34": _default_recombination_rate,
    "35": _default_recombination_rate,
    "36": _default_recombination_rate,
    "37": _default_recombination_rate,
    "38": _default_recombination_rate,
    "39": _default_recombination_rate,
    "40": _default_recombination_rate,
}


# value used in Lavretsky et al 2020, as obtained for nuclear genes in other
# ducks by Peters, Zhuravlev, Fefelov, Humphries, & Omland, 2008
# The per-year rate of 1.2e-9 was converted to a per-generation rate
# by using assumptions on average generation time of 4 years (see below)
_overall_mutation_rate = 4.83e-9  # per generation

_mutation_rate = {
    "1": _overall_mutation_rate,
    "2": _overall_mutation_rate,
    "3": _overall_mutation_rate,
    "4": _overall_mutation_rate,
    "5": _overall_mutation_rate,
    "6": _overall_mutation_rate,
    "7": _overall_mutation_rate,
    "8": _overall_mutation_rate,
    "9": _overall_mutation_rate,
    "10": _overall_mutation_rate,
    "11": _overall_mutation_rate,
    "12": _overall_mutation_rate,
    "13": _overall_mutation_rate,
    "14": _overall_mutation_rate,
    "15": _overall_mutation_rate,
    "16": _overall_mutation_rate,
    "17": _overall_mutation_rate,
    "18": _overall_mutation_rate,
    "19": _overall_mutation_rate,
    "20": _overall_mutation_rate,
    "21": _overall_mutation_rate,
    "22": _overall_mutation_rate,
    "23": _overall_mutation_rate,
    "24": _overall_mutation_rate,
    "25": _overall_mutation_rate,
    "26": _overall_mutation_rate,
    "27": _overall_mutation_rate,
    "28": _overall_mutation_rate,
    "29": _overall_mutation_rate,
    "30": _overall_mutation_rate,
    "Z": _overall_mutation_rate,
    "31": _overall_mutation_rate,
    "32": _overall_mutation_rate,
    "33": _overall_mutation_rate,
    "34": _overall_mutation_rate,
    "35": _overall_mutation_rate,
    "36": _overall_mutation_rate,
    "37": _overall_mutation_rate,
    "38": _overall_mutation_rate,
    "39": _overall_mutation_rate,
    "40": _overall_mutation_rate,
}

_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
    citations=[
        _LavretskyEtAl2019,
        _HuangEtAl2006,
        _LavretskyEtAl2020,
        _ZhuEtAl2021,
    ],
)
stdpopsim.utils.append_common_synonyms(_genome)

_species = stdpopsim.Species(
    id="AnaPla",
    ensembl_id="anas_platyrhynchos",
    name="Anas platyrhynchos",
    common_name="Mallard",
    # description="The 'mallard' species complex consists of 14 hybridizing and "
    # "recently diverged species living around the world, ranging from the holarctic "
    # "mallard with >15M individuals today in North America alone to "
    # "endangered endemics in Hawaii and New Zealand. The assembly, "
    # "recombination rates, and default Ne were estimtaed with wild Chinese "
    # "mallards.",
    genome=_genome,
    # generation time estimate from Lavertsky et al. (2020):
    # Generation time (G) was calculated as G = \alpha + (s/(1 − s)),
    # where \alpha is the age of maturity and s is the expected adult
    # survival rate (Sather et al., 2005). The age of maturity for mallard-
    # like ducks generally is one year (i.e., \alpha = 1), and the average
    # adult survival rate is 0.54 (range: 0.34–0.74) and 0.54 (range: 0.4–0.70)
    # for mallards and black ducks, respectively (Nichols, Obrecht, & Hines, 1987).
    # Using an overall survival rate average of 0.54 for the two species, we
    # estimated the generation time to be 4.0 years.
    generation_time=4,
    # choosing Ne based on theta = 4 Ne u from Guo et al 2021
    # theta = 0.003 (Figure 1), u as above (the paper uses a rate from chicken)
    population_size=156000,
    citations=[
        _LavretskyEtAl2020,
        _GuoEtAl2020,
    ],
)

stdpopsim.register_species(_species)
