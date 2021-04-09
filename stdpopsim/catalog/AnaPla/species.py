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
    doi="https://doi.org/10.1111/mec.15343",
    reasons={stdpopsim.CiteReason.GEN_TIME, stdpopsim.CiteReason.MUT_RATE},
)

_GuoEtAl2020 = stdpopsim.Citation(
    author="Guo et al.",
    year=2020,
    doi="https://doi.org/10.24272/j.issn.2095-8137.2020.133",
    reasons={stdpopsim.CiteReason.POP_SIZE},
)


# def _print_recomb_rates():
#     """
#     Code to produce the following rates, which follows Lavretsky et al 2019
#     in taking the total map lengths from chromosomes 1-7 from a microsat study
#     (Huang et al 2006) to give average rates. This assumes the total lengths for
#     chromsomes 1-7 is as in Huang et al, and applies the average per-bp rate to
#     the rest of the chromosomes.
#     """
#     from genome_data import data as mallard
#
#     chroms = mallard["chromosomes"]
#     # from Huang et al 2006
#     huang_chroms = {
#         "1": {"M": 3.17},
#         "2": {"M": 2.26},
#         "3": {"M": 1.12},
#         "4": {"M": 0.93},
#         "5": {"M": 0.79},
#         "6": {"M": 1.20},
#         "7": {"M": 0.98},
#     }
#     for c in huang_chroms:
#         h = huang_chroms[c]
#         h["bp"] = mallard["chromosomes"][c]["length"]
#         h["rate"] = h["M"] / h["bp"]
#     total_M = sum([x["M"] for x in huang_chroms.values()])
#     total_bp = sum([x["bp"] for x in huang_chroms.values()])
#     mean_rate = total_M / total_bp
#     print(f"_default_recombination_rate = {mean_rate}")
#     print("_recombination_rate = {")
#     for c in chroms:
#         if c in huang_chroms:
#             print(f"\"{c}\": {huang_chroms[c]['rate']:.2e},")
#         else:
#             print(f'"{c}": _default_recombination_rate,')
#     print("}")


_default_recombination_rate = 1.47e-08


_recombination_rate = {
    "1": 1.52e-08,
    "2": 1.39e-08,
    "3": 9.35e-09,
    "4": 1.20e-08,
    "5": 1.22e-08,
    "6": 3.03e-08,
    "7": 2.59e-08,
    "8": _default_recombination_rate,
    "9": _default_recombination_rate,
    "10": _default_recombination_rate,
    "11": _default_recombination_rate,
    "12": _default_recombination_rate,
    "13": _default_recombination_rate,
    "14": _default_recombination_rate,
    "15": _default_recombination_rate,
    "16": _default_recombination_rate,
    "17": _default_recombination_rate,
    "18": _default_recombination_rate,
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


# value used in Lavretsky et al 2019, as obtained for nuclear genes in other
# ducks by Peters, Zhuravlev, Fefelov, Humphries, & Omland, 2008
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
    ],
)

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
    generation_time=4,
    # choosing Ne based on theta = 4 Ne u from Guo et al 2020
    # theta = 0.003 (Figure 1), u as above (the paper uses a rate from chicken)
    population_size=156000,
    citations=[
        _LavretskyEtAl2019,
        _GuoEtAl2020,
    ],
)

stdpopsim.register_species(_species)
