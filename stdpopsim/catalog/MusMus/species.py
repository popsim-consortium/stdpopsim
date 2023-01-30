import stdpopsim

from . import genome_data

# Chromosome-specific recombination rate
_recombination_rate = {
    "1": 2.75e-09,
    "2": 3.32e-09,
    "3": 3.51e-09,
    "4": 3.56e-09,
    "5": 3.86e-09,
    "6": 3.44e-09,
    "7": 3.15e-09,
    "8": 4.34e-09,
    "9": 4.87e-09,
    "10": 4.28e-09,
    "11": 5.71e-09,
    "12": 4.32e-09,
    "13": 5.08e-09,
    "14": 4.19e-09,
    "15": 5.83e-09,
    "16": 7.62e-09,
    "17": 6.30e-09,
    "18": 7.13e-09,
    "19": 1.45e-08,
    "X": 2.72e-09,
    "Y": 0,
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
    "X": _species_ploidy,
    "Y": 1,
    "MT": 1,
}

# Generic and chromosome-specific mutation rate
_overall_rate = (5.4e-09,)
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
    "X": _overall_rate,
    "Y": _overall_rate,
    "MT": 3.7e-08,
}

_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    citations=[
        stdpopsim.Citation(
            author="Cox et al.",
            year=2009,
            doi="https://doi.org/10.1534/genetics.109.105486",
            reasons={stdpopsim.CiteReason.REC_RATE},
        ),
        stdpopsim.Citation(
            author="Uchimura et al.",
            year=2015,
            doi="http://www.genome.org/cgi/doi/10.1101/gr.186148.114",
            reasons={stdpopsim.CiteReason.MUT_RATE},
        ),
        stdpopsim.Citation(
            author="Goios et al.",
            year=2007,
            doi="http://www.genome.org/cgi/doi/10.1101/gr.5941007",
            reasons={stdpopsim.CiteReason.MUT_RATE},
        ),
    ],
)

_species = stdpopsim.Species(
    id="MusMus",
    ensembl_id="mus_musculus",
    name="Mus musculus",
    common_name="Mouse",
    genome=_genome,
    generation_time=0.75,
    population_size=500000,
    citations=[
        stdpopsim.Citation(
            author="Fujiwara et al.",
            year=2022,
            doi="https://doi.org/10.1093/gbe/evac068",
            reasons={stdpopsim.CiteReason.GEN_TIME, stdpopsim.CiteReason.POP_SIZE},
        ),
        stdpopsim.Citation(
            author="Phifer-Rixey, Harr et al.",
            year=2020,
            doi="https://doi.org/10.1186/s12862-020-01666-9",
            reasons={stdpopsim.CiteReason.GEN_TIME},
        ),
        stdpopsim.Citation(
            author="Phifer-Rixey, Bonhomme, et al.",
            year=2012,
            doi="https://doi.org/10.1093/molbev/mss105",
            reasons={stdpopsim.CiteReason.POP_SIZE},
        ),
    ],
)

stdpopsim.register_species(_species)
