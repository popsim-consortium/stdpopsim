import stdpopsim

from . import genome_data


# Chromosome-specific recombination rate
# In order to derive an estimate of recombination rate for each chromosome
# we use information provided in Table 1 of Littrell et al. (2018), specifically
# this is an update of Jensen-Seaman map in 2004, whose values are highly similar.
# https://pmc.ncbi.nlm.nih.gov/articles/PMC6027877/

_recombination_rate = {
    "1": 0.50e-09,
    "2": 0.48e-09,
    "3": 0.59e-09,
    "4": 0.57e-09,
    "5": 0.61e-09,
    "6": 0.60e-09,
    "7": 0.66e-09,
    "8": 0.68e-09,
    "9": 0.64e-09,
    "10": 0.81e-09,
    "11": 0.68e-09,
    "12": 1.01e-09,
    "13": 0.57e-09,
    "14": 0.63e-09,
    "15": 0.65e-09,
    "16": 0.68e-09,
    "17": 0.75e-09,
    "18": 0.67e-09,
    "19": 0.88e-09,
    "20": 0.91e-09,
    "X": 0.37e-09,
    "Y": 0,
    "MT": 0,
}

# Chromosome-specific mutation rate
# According to the paper below the mutation rate in MT
# is about 10-20 times higher, therefore we set it to 15 times

_overall_rate = 2.96e-9
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
    "20": _overall_rate,
    "X": _overall_rate,
    "Y": _overall_rate,
    "MT": _overall_rate
    * 15,  # https://journals.physiology.org/doi/full/10.1152/ajpcell.00234.2006
}


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
        stdpopsim.Citation(
            author="Littrell et al.",
            year=2018,
            doi="https://doi.org/10.1534/g3.118.200187",
            reasons={stdpopsim.CiteReason.ASSEMBLY},
        ),
        stdpopsim.Citation(
            author="Deinum et al.",
            year=2015,
            doi="https://doi.org/10.1093/molbev/msv126",
            reasons={stdpopsim.CiteReason.REC_RATE},
        ),
    ],
)

_species = stdpopsim.Species(
    id="RatNor",
    ensembl_id="rattus_norvegicus",
    name="Rattus norvegicus",
    common_name="Rat",
    genome=_genome,
    ploidy=_species_ploidy,
    generation_time=0.5,
    population_size=1.24e5,
    citations=[
        stdpopsim.Citation(
            author="Deinum et al.",
            year=2015,
            doi="https://doi.org/10.1093/molbev/msv126",
            reasons={stdpopsim.CiteReason.POP_SIZE},
        ),
    ],
)

stdpopsim.register_species(_species)
