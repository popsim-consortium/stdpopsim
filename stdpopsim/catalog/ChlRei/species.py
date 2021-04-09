import stdpopsim

from . import genome_data

# IMPORTANT: these are recombination rate values with the rate of sex
# incorporated - e.g. they are lower than the physical recombination rates on account
# of being the 'realized' recombination rates in nature, since C. reinhardtii only
# undergoes sex/recombination once every ~840 generations (Hasan and Ness 2020).
# when simulating natural populations with these values, the infrequency
# of sex should be implicitly taken into account. otherwise, to get estimates
# of physical rates, divide all values by 0.001194. see Discussion of Hasan
# and Ness 2020 for more details
_recombination_rate = {
    "1": 1.21e-10,
    "2": 1.49e-10,
    "3": 1.52e-10,
    "4": 1.47e-10,
    "5": 1.70e-10,
    "6": 1.17e-10,
    "7": 8.66e-11,
    "8": 1.39e-10,
    "9": 1.12e-10,
    "10": 1.97e-10,
    "11": 1.63e-10,
    "12": 9.15e-11,
    "13": 1.43e-10,
    "14": 1.9e-10,
    "15": 3.93e-10,
    "16": 1.71e-10,
    "17": 1.83e-10,
}

_overall_rate = 9.63 * 1e-10  # SNM rate
# there is also an indel mutation rate in Ness 2015 if that is eventually supported
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
}

_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
    citations=[
        stdpopsim.Citation(
            author="Merchant et al",
            year=2007,
            doi="https://doi.org/10.1126/science.1143609",
            reasons={stdpopsim.CiteReason.ASSEMBLY},  # v5 - v6 assembly still en route!
        ),
        stdpopsim.Citation(
            author="Hasan and Ness",
            year=2020,
            doi="https://doi.org/10.1093/gbe/evaa057",
            reasons={stdpopsim.CiteReason.REC_RATE},
        ),
        stdpopsim.Citation(
            author="Ness et al",
            year=2015,
            doi="https://doi.org/10.1101/gr.191494.115",
            reasons={stdpopsim.CiteReason.MUT_RATE},
        ),
    ],
)

_species = stdpopsim.Species(
    id="ChlRei",
    ensembl_id="chlamydomonas_reinhardtii",
    name="Chlamydomonas reinhardtii",
    common_name="Chlamydomonas reinhardtii",
    genome=_genome,
    generation_time=1 / 876,
    population_size=1.4 * 1e-7,
    citations=[
        stdpopsim.Citation(
            author="Ness et al.",
            year=2016,
            doi="https://doi.org/10.1093/molbev/msv272",
            reasons={stdpopsim.CiteReason.POP_SIZE},  # Quebec population
        ),
        stdpopsim.Citation(
            author="Vitova et al.",
            year=2011,
            doi="https://doi.org/10.1007/s00425-011-1427-7",
            reasons={stdpopsim.CiteReason.GEN_TIME},
        ),
    ],
)

stdpopsim.register_species(_species)
