import stdpopsim

from . import genome_data


# Chromosome-specific recombination rate
# In order to derive an estimate of recombination rate for each chromosome
# we use information provided in Table 1 of Littrell et al. (2018), specifically
# this is an update of Jensen-Seaman map in 2004, whose values are highly similar.
# https://pmc.ncbi.nlm.nih.gov/articles/PMC6027877/

_recombination_rate = {
    "1": 5.0e-09,
    "2": 4.8e-09,
    "3": 5.9e-09,
    "4": 5.7e-09,
    "5": 6.1e-09,
    "6": 6.0e-09,
    "7": 6.6e-09,
    "8": 6.8e-09,
    "9": 6.4e-09,
    "10": 8.1e-09,
    "11": 6.8e-09,
    "12": 10.1e-09,
    "13": 5.7e-09,
    "14": 6.3e-09,
    "15": 6.5e-09,
    "16": 6.8e-09,
    "17": 7.5e-09,
    "18": 6.7e-09,
    "19": 8.8e-09,
    "20": 9.1e-09,
    "X": 3.7e-09,
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
    "MT": _overall_rate * 15,
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
            author="Howe et al.",
            year=2021,
            doi="https://doi.org/10.12688/wellcomeopenres.16854.1",
            reasons={stdpopsim.CiteReason.ASSEMBLY},
        ),
        stdpopsim.Citation(
            author="Littrell et al.",
            year=2018,
            doi="https://doi.org/10.1534/g3.118.200187",
            reasons={stdpopsim.CiteReason.REC_RATE},
        ),
        stdpopsim.Citation(
            author="Deinum et al.",
            year=2015,
            doi="https://doi.org/10.1093/molbev/msv126",
            reasons={stdpopsim.CiteReason.MUT_RATE},
        ),
        # Mitochondrial mutation rate
        stdpopsim.Citation(
            author="Schlick et al.",
            year=2006,
            doi="https://doi.org/10.1152/ajpcell.00234.2006",
            reasons={stdpopsim.CiteReason.MUT_RATE},
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
        stdpopsim.Citation(
            author="Deinum et al.",
            year=2015,
            doi="https://doi.org/10.1093/molbev/msv126",
            reasons={stdpopsim.CiteReason.GEN_TIME},
        ),
    ],
)

stdpopsim.register_species(_species)
