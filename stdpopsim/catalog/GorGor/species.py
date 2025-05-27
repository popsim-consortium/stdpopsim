import stdpopsim
from . import genome_data

_recombination_rate = dict()
_mutation_rate = dict()
_ploidy = dict()

# genome-wide recombination rate
# associated with all recombining chromosomes
_Ne = 25200
_mean_recombination_rate = 1.193e-08
_recombination_rate_data = {
    "1": _mean_recombination_rate,
    "2A": _mean_recombination_rate,
    "2B": _mean_recombination_rate,
    "3": _mean_recombination_rate,
    "4": _mean_recombination_rate,
    "5": _mean_recombination_rate,
    "6": _mean_recombination_rate,
    "7": _mean_recombination_rate,
    "8": _mean_recombination_rate,
    "9": _mean_recombination_rate,
    "10": _mean_recombination_rate,
    "11": _mean_recombination_rate,
    "12": _mean_recombination_rate,
    "13": _mean_recombination_rate,
    "14": _mean_recombination_rate,
    "15": _mean_recombination_rate,
    "16": _mean_recombination_rate,
    "17": _mean_recombination_rate,
    "18": _mean_recombination_rate,
    "19": _mean_recombination_rate,
    "20": _mean_recombination_rate,
    "21": _mean_recombination_rate,
    "22": _mean_recombination_rate,
    "X": _mean_recombination_rate,
    "MT": 0,
}

# genome-wide average mutation rate
# associated with all chromosomes
_mean_mutation_rate = 1.235e-8
_mutation_rate = {
    "1": _mean_mutation_rate,
    "2A": _mean_mutation_rate,
    "2B": _mean_mutation_rate,
    "3": _mean_mutation_rate,
    "4": _mean_mutation_rate,
    "5": _mean_mutation_rate,
    "6": _mean_mutation_rate,
    "7": _mean_mutation_rate,
    "8": _mean_mutation_rate,
    "9": _mean_mutation_rate,
    "10": _mean_mutation_rate,
    "11": _mean_mutation_rate,
    "12": _mean_mutation_rate,
    "13": _mean_mutation_rate,
    "14": _mean_mutation_rate,
    "15": _mean_mutation_rate,
    "16": _mean_mutation_rate,
    "17": _mean_mutation_rate,
    "18": _mean_mutation_rate,
    "19": _mean_mutation_rate,
    "20": _mean_mutation_rate,
    "21": _mean_mutation_rate,
    "22": _mean_mutation_rate,
    "X": _mean_mutation_rate,
    "MT": _mean_mutation_rate,
}

# species ploidy and chromosome-specific ploidy
_species_ploidy = 2
_ploidy = {
    "1": _species_ploidy,
    "2A": _species_ploidy,
    "2B": _species_ploidy,
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
    "22": _species_ploidy,
    "X": _species_ploidy,
    "MT": 1,
}

_chromosomes = []
for name, data in genome_data.data["chromosomes"].items():
    _chromosomes.append(
        stdpopsim.Chromosome(
            id=name,
            length=data["length"],
            synonyms=data["synonyms"],
            mutation_rate=_mutation_rate[name],
            recombination_rate=_recombination_rate_data[name],
            ploidy=_ploidy[name],
        )
    )

_genome = stdpopsim.Genome(
    chromosomes=_chromosomes,
    assembly_name=genome_data.data["assembly_name"],
    assembly_accession=genome_data.data["assembly_accession"],
    citations=[
        stdpopsim.Citation(
            author="Besenbacher et al.",
            year=2019,
            doi="https://doi.org/10.1038/s41559-018-0778-x",
            reasons={stdpopsim.CiteReason.MUT_RATE},
        ),
        stdpopsim.Citation(
            author="Stevison et al.",
            year=2015,
            doi="https://doi.org/10.1093/molbev/msv331",
            reasons={stdpopsim.CiteReason.REC_RATE},
        ),
        stdpopsim.Citation(
            doi="https://doi.org/10.1126/science.aae0344",
            year=2016,
            author="Gordon et al.",
            reasons={stdpopsim.CiteReason.ASSEMBLY},
        ),
    ],
)
stdpopsim.utils.append_common_synonyms(_genome)

_species = stdpopsim.Species(
    id="GorGor",
    ensembl_id="gorilla_gorilla",
    name="Gorilla gorilla",
    common_name="Gorilla",
    genome=_genome,
    generation_time=19.0,
    population_size=_Ne,
    ploidy=_species_ploidy,
    citations=[
        stdpopsim.Citation(
            doi="https://doi.org/10.1038/s41559-018-0778-x",
            year=2019,
            author="Besenbacher et al.",
            reasons={stdpopsim.CiteReason.GEN_TIME},
        ),
        stdpopsim.Citation(
            doi="https://doi.org/10.1534/genetics.166.3.1375",
            year=2004,
            author="Yu et al.",
            reasons={stdpopsim.CiteReason.POP_SIZE},
        ),
    ],
)

stdpopsim.register_species(_species)
