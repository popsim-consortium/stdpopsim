import stdpopsim

from . import genome_data

# IMPORTANT: These are recombination rate values with the rate of sex
# incorporated - e.g. they are lower than the physical recombination rates on account
# of being the 'realized' recombination rates in nature, since C. reinhardtii only
# undergoes sex/recombination once every ~840 generations (Hasan and Ness 2020).
# When simulating natural populations with these values, the infrequency
# of sex should be implicitly taken into account. Otherwise, to get estimates
# of physical rates, divide all values by 0.001194. See Discussion of Hasan
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

# Generic and chromosome-specific ploidy
_species_ploidy = 1
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
}

_mutation_rate = {
    "1": 9.74e-10,
    "2": 8.62e-10,
    "3": 9.5e-10,
    "4": 9.66e-10,
    "5": 1.17e-9,
    "6": 9.12e-10,
    "7": 9.14e-10,
    "8": 8.98e-10,
    "9": 9.17e-10,
    "10": 9.27e-10,
    "11": 1.03e-9,
    "12": 9.55e-10,
    "13": 7.56e-10,
    "14": 8.96e-10,
    "15": 6.91e-10,
    "16": 9.59e-10,
    "17": 1.05e-9,
}

_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
    ploidy=_ploidy,
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
            doi="https://doi.org/10.6084/m9.figshare.14608239.v1",
            reasons={stdpopsim.CiteReason.REC_RATE},
        ),
        stdpopsim.Citation(
            author="Ness et al",
            year=2015,
            doi="https://doi.org/10.6084/m9.figshare.14700156.v1",
            reasons={stdpopsim.CiteReason.MUT_RATE},
        ),
    ],
)
stdpopsim.utils.append_common_synonyms(_genome)

_species = stdpopsim.Species(
    id="ChlRei",
    ensembl_id="chlamydomonas_reinhardtii",
    name="Chlamydomonas reinhardtii",
    common_name="Chlamydomonas reinhardtii",
    separate_sexes=False,
    genome=_genome,
    generation_time=1 / 876,
    population_size=1.4 * 1e-7,
    ploidy=_species_ploidy,
    citations=[
        stdpopsim.Citation(
            author="Ness et al.",
            year=2016,
            doi="https://doi.org/10.1093/molbev/msv272",
            reasons={stdpopsim.CiteReason.POP_SIZE},  # Quebec population
        ),
        stdpopsim.Citation(
            author="Vítová et al",
            year=2011,
            doi="https://doi.org/10.1007/s00425-011-1427-7",
            reasons={stdpopsim.CiteReason.GEN_TIME},
        ),
    ],
)

stdpopsim.register_species(_species)
