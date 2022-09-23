import stdpopsim

from . import genome_data

_kuderna2017 = stdpopsim.Citation(
    doi="https://doi.org/10.1093/gigascience/gix098",
    year=2017,
    author="Kuderna et al.",
    reasons={stdpopsim.CiteReason.ASSEMBLY},
)

_besenbacher2019 = stdpopsim.Citation(
    doi="https://doi.org/10.1038/s41559-018-0778-x",
    year=2019,
    author="Besenbacher et al.",
    reasons={stdpopsim.CiteReason.MUT_RATE},
)

_stevison2015 = stdpopsim.Citation(
    doi="https://doi.org/10.1093/molbev/msv331",
    year=2015,
    author="Stevison et al.",
    reasons={stdpopsim.CiteReason.REC_RATE},
)

_langergraber2012 = stdpopsim.Citation(
    doi="https://doi.org/10.1073/pnas.1211740109",
    year=2012,
    author="Langergraber et al.",
    reasons={stdpopsim.CiteReason.GEN_TIME},
)

# Generic and chromosome-specific ploidy
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
    "Y": 1,
}

_chromosomes = []
for name, data in genome_data.data["chromosomes"].items():
    _chromosomes.append(
        stdpopsim.Chromosome(
            id=name,
            length=data["length"],
            synonyms=data["synonyms"],
            mutation_rate=1.6e-8,
            recombination_rate=1.2e-8,
            ploidy=_ploidy[name],
        )
    )

_genome = stdpopsim.Genome(
    chromosomes=_chromosomes,
    assembly_name=genome_data.data["assembly_name"],
    assembly_accession=genome_data.data["assembly_accession"],
    citations=[
        _kuderna2017.because(stdpopsim.CiteReason.ASSEMBLY),
        _besenbacher2019.because(stdpopsim.CiteReason.MUT_RATE),
        _stevison2015.because(stdpopsim.CiteReason.REC_RATE),
    ],
)
stdpopsim.utils.append_common_synonyms(_genome)

_species = stdpopsim.Species(
    id="PanTro",
    ensembl_id="pan_troglodytes",
    name="Pan troglodytes",
    common_name="Chimpanzee",
    genome=_genome,
    generation_time=24.6,
    population_size=16781,
    ploidy=_species_ploidy,
    citations=[
        _langergraber2012.because(stdpopsim.CiteReason.GEN_TIME),
        _stevison2015.because(stdpopsim.CiteReason.POP_SIZE),
    ],
)

stdpopsim.register_species(_species)
