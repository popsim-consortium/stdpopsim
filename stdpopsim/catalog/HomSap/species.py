import stdpopsim

from . import genome_data

# Mean recombination rate was computed from the "HapMapII_GRCh38" map
# by dividing total map length by total sequence length over non-missing
# windows. (See: `msprime.RateMap.mean_rate`)
_recombination_rate_data = {
    "1": 1.1523470111585671e-08,
    "2": 1.1042947599187705e-08,
    "3": 1.1258496243179589e-08,
    "4": 1.148200397734684e-08,
    "5": 1.1244324085272748e-08,
    "6": 1.1265873463291707e-08,
    "7": 1.177126324439884e-08,
    "8": 1.1604882118427952e-08,
    "9": 1.219866363355389e-08,
    "10": 1.3333710355319759e-08,
    "11": 1.1721276661036567e-08,
    "12": 1.3098075878935492e-08,
    "13": 1.3061005794618377e-08,
    "14": 1.3629801391728677e-08,
    "15": 1.738760443308643e-08,
    "16": 1.4831470988338458e-08,
    "17": 1.55382549450431e-08,
    "18": 1.4645528336185072e-08,
    "19": 1.8384797447801312e-08,
    "20": 1.6788552013988714e-08,
    "21": 1.724434082825555e-08,
    "22": 2.1057233894035443e-08,
    "X": 1.1848323598825364e-08,
    "Y": 0.0,
    "MT": 0.0,
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
    "20": _species_ploidy,
    "21": _species_ploidy,
    "22": _species_ploidy,
    "X": _species_ploidy,
    "Y": 1,
    "MT": 1,
}

_genome2001 = stdpopsim.Citation(
    doi="http://dx.doi.org/10.1038/35057062",
    year=2001,
    author="International Human Genome Sequencing Consortium",
    reasons={stdpopsim.CiteReason.ASSEMBLY},
)

_hapmap2007 = stdpopsim.Citation(
    doi="https://doi.org/10.1038/nature06258",
    year=2007,
    author="The International HapMap Consortium",
)

_takahata1993 = stdpopsim.Citation(
    doi="https://doi.org/10.1093/oxfordjournals.molbev.a039995",
    year=1993,
    author="Takahata",
    reasons={stdpopsim.CiteReason.POP_SIZE},
)

_tian2019 = stdpopsim.Citation(
    doi="https://doi.org/10.1016/j.ajhg.2019.09.012",
    year=2019,
    author="Tian, Browning, and Browning",
    reasons={stdpopsim.CiteReason.MUT_RATE},
)

_tremblay2000 = stdpopsim.Citation(
    doi="https://doi.org/10.1086/302770",
    year=2000,
    author="Tremblay and Vézina",
    reasons={stdpopsim.CiteReason.GEN_TIME},
)

_chromosomes = []
for name, data in genome_data.data["chromosomes"].items():
    _chromosomes.append(
        stdpopsim.Chromosome(
            id=name,
            length=data["length"],
            synonyms=data["synonyms"],
            mutation_rate=1.29e-8,
            recombination_rate=_recombination_rate_data[name],
            ploidy=_ploidy[name],
        )
    )

_genome = stdpopsim.Genome(
    chromosomes=_chromosomes,
    assembly_name=genome_data.data["assembly_name"],
    assembly_accession=genome_data.data["assembly_accession"],
    citations=[
        _genome2001,
        _tian2019.because(stdpopsim.CiteReason.MUT_RATE),
        _hapmap2007.because(stdpopsim.CiteReason.REC_RATE),
    ],
)
stdpopsim.utils.append_common_synonyms(_genome)

_species = stdpopsim.Species(
    id="HomSap",
    ensembl_id="homo_sapiens",
    name="Homo sapiens",
    common_name="Human",
    separate_sexes=True,
    genome=_genome,
    generation_time=30,
    population_size=10**4,
    ploidy=_species_ploidy,
    citations=[
        _tremblay2000.because(stdpopsim.CiteReason.GEN_TIME),
        _takahata1993.because(stdpopsim.CiteReason.POP_SIZE),
    ],
)

stdpopsim.register_species(_species)
