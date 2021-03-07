import stdpopsim

from . import genome_data

# mean_recombination_rate was computed across all windows of the HapMapII_GRCh37 map
_recombination_rate_data = {
    "1": 1.1485597641285933e-08,
    "2": 1.1054289277533446e-08,
    "3": 1.1279585624662551e-08,
    "4": 1.1231162636001008e-08,
    "5": 1.1280936570022824e-08,
    "6": 1.1222852661225285e-08,
    "7": 1.1764614397655721e-08,
    "8": 1.1478465778920576e-08,
    "9": 1.1780701596308656e-08,
    "10": 1.3365134257075317e-08,
    "11": 1.1719334320833283e-08,
    "12": 1.305017186986983e-08,
    "13": 1.0914860554958317e-08,
    "14": 1.119730771394731e-08,
    "15": 1.3835785893339787e-08,
    "16": 1.4834607113882717e-08,
    "17": 1.582489036239487e-08,
    "18": 1.5075956950023575e-08,
    "19": 1.8220141872466202e-08,
    "20": 1.7178269031631664e-08,
    "21": 1.3045214034879191e-08,
    "22": 1.4445022767788226e-08,
    "X": 1.164662223273842e-08,
    "Y": 0.0,
    "MT": 0.0,
}

_genome2001 = stdpopsim.Citation(
    doi="http://dx.doi.org/10.1038/35057062",
    year="2001",
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
    year="1993",
    author="Takahata",
    reasons={stdpopsim.CiteReason.POP_SIZE},
)

_tian2019 = stdpopsim.Citation(
    doi="https://doi.org/10.1016/j.ajhg.2019.09.012",
    year="2019",
    author="Tian, Browning, and Browning",
    reasons={stdpopsim.CiteReason.MUT_RATE},
)

_tremblay2000 = stdpopsim.Citation(
    doi="https://doi.org/10.1086/302770",
    year="2000",
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
        )
    )

_genome = stdpopsim.Genome(
    chromosomes=_chromosomes,
    mutation_rate_citations=[_tian2019.because(stdpopsim.CiteReason.MUT_RATE)],
    recombination_rate_citations=[_hapmap2007.because(stdpopsim.CiteReason.REC_RATE)],
    assembly_name=genome_data.data["assembly_name"],
    assembly_accession=genome_data.data["assembly_accession"],
    assembly_citations=[_genome2001],
)

_species = stdpopsim.Species(
    id="HomSap",
    name="Homo sapiens",
    common_name="Human",
    genome=_genome,
    generation_time=30,
    generation_time_citations=[_tremblay2000.because(stdpopsim.CiteReason.GEN_TIME)],
    population_size=10 ** 4,
    population_size_citations=[_takahata1993.because(stdpopsim.CiteReason.POP_SIZE)],
)

stdpopsim.register_species(_species)
