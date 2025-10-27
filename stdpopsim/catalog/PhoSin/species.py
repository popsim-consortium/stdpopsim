import stdpopsim

from . import genome_data

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
    "20": _species_ploidy,
    "21": _species_ploidy,
    "X": _species_ploidy,
}

##########################
# There is no direct estimate of recombination rate
# for Vaquita or any closely related species.
# Other studies that simulate data for Vaquita assume
# a default mutation rate of 1e-8, which is what we
# use here.
# For example, Morin et al. (2021) assume a recombination
# rate of 1e-8 when running PSMC, and Robinson et al. (2022)
# use a recombination rate of 1e-8 when simulating data using
# SLiM.
##########################
_overall_recombination_rate = 1e-8
_recombination_rate = {
    "1": _overall_recombination_rate,
    "2": _overall_recombination_rate,
    "3": _overall_recombination_rate,
    "4": _overall_recombination_rate,
    "5": _overall_recombination_rate,
    "6": _overall_recombination_rate,
    "7": _overall_recombination_rate,
    "8": _overall_recombination_rate,
    "9": _overall_recombination_rate,
    "10": _overall_recombination_rate,
    "11": _overall_recombination_rate,
    "12": _overall_recombination_rate,
    "13": _overall_recombination_rate,
    "14": _overall_recombination_rate,
    "15": _overall_recombination_rate,
    "16": _overall_recombination_rate,
    "17": _overall_recombination_rate,
    "18": _overall_recombination_rate,
    "19": _overall_recombination_rate,
    "20": _overall_recombination_rate,
    "21": _overall_recombination_rate,
    "X": _overall_recombination_rate,
}

##########################
# Genome-wide average mutation rate is estimated by Robinson et al. (2022)
# by using divergence from two closely related species: harbor porpoise
# (Phocoena phocoena) and Indo-Pacific finless porpoise (Neophocaena phocaenoides).
# Using these two species and a range of assumptions results in a range of
# plausible estimates depicted in Fig S12 in the paper.
# Eventually, the representative value used in most analyses in the paper was
# 5.83e-9, which is on the high side of the range, to take into account previous
# higher estimates of the mutation rate, cited by Morin et al. (2021).
##########################
_overall_mutation_rate = 5.83e-9
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
    "X": _overall_mutation_rate,
}

_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
    ploidy=_ploidy,
    citations=[
        stdpopsim.Citation(
            author="Morin et al.",
            year=2021,
            doi="https://doi.org/10.1111/1755-0998.13284",
            reasons={
                stdpopsim.CiteReason.ASSEMBLY,
                stdpopsim.CiteReason.REC_RATE,
            },
        ),
        stdpopsim.Citation(
            author="Robinson et al.",
            year=2022,
            doi="https://doi.org/10.1126/science.abm1742",
            reasons={stdpopsim.CiteReason.MUT_RATE},
        ),
    ],
)

_species = stdpopsim.Species(
    id="PhoSin",
    ensembl_id="phocoena_sinus",
    name="Phocoena sinus",
    common_name="Vaquita",
    separate_sexes=True,
    genome=_genome,
    ploidy=_species_ploidy,
    ##########################
    # Robinson et al. (2022) and Morin et al. (2021) assumed an average
    # generation time of 11.9 years in their analyses/calibrations.
    # They cite Taylor et al. (2007) as a primary source.
    # Table 1 in that paper specifies estimates of generation
    # times for 58 cetacean species. It seems like the value selected was
    # the one inferred for the sister species of harbor porpoise (Phocoena
    # phocoena). The value estimated for Vaquita is slightly smaller (11.4),
    # but maybe it reflects a recent reducion in generation time (???)
    ##########################
    generation_time=11.9,
    ##########################
    # Demographic models inferred for Vaquita by Morin et al. (2021) and
    # Robinson et al. (2022) typically estimate small Ne of 2,000 - 5,000
    # in the past 500,000 years, likely followed by some sharp bottleneck.
    # We selected as a representative Ne a value of 3500, which is also
    # consistent with the average Ne inferred in the 2-epoch model of
    # Robinson et al. (2022)
    ##########################
    population_size=3500,
    citations=[
        stdpopsim.Citation(
            author="Robinson et al.",
            year=2022,
            doi="https://doi.org/10.1126/science.abm1742",
            reasons={
                stdpopsim.CiteReason.POP_SIZE,
                stdpopsim.CiteReason.GEN_TIME,
            },
        ),
        stdpopsim.Citation(
            author="Taylor et al.",
            year=2007,
            doi="https://aquadocs.org/handle/1834/41281",
            reasons={stdpopsim.CiteReason.GEN_TIME},
        ),
    ],
)

stdpopsim.register_species(_species)
