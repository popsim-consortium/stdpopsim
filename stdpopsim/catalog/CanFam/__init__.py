import stdpopsim
from . import genome_data

# Recombination rates are the means calculated from the Campbell2016 map. ChrX
# isn't included in that genetic map, so the weighted genome-wide average is
# provided here instead.
_recombination_rate_data = {
    "1": 7.636001498077e-09,
    "2": 8.798516284201015e-09,
    "3": 8.000868549297668e-09,
    "4": 8.052302321538563e-09,
    "5": 9.344333839828404e-09,
    "6": 8.192191165856811e-09,
    "7": 7.293465322857858e-09,
    "8": 8.291312822214086e-09,
    "9": 9.287722798450121e-09,
    "10": 9.107151930130527e-09,
    "11": 7.639449804041734e-09,
    "12": 7.761061128067527e-09,
    "13": 8.413015225299971e-09,
    "14": 9.028123091776737e-09,
    "15": 7.856754982030383e-09,
    "16": 8.614063697877811e-09,
    "17": 9.718834385033715e-09,
    "18": 1.029925446520415e-08,
    "19": 1.0425051705950442e-08,
    "20": 9.990971109010329e-09,
    "21": 1.0339033506095636e-08,
    "22": 8.61504863805126e-09,
    "23": 9.1266350068389e-09,
    "24": 1.1146043069685935e-08,
    "25": 1.1543746646856006e-08,
    "26": 1.2084571786110922e-08,
    "27": 1.1260328390864541e-08,
    "28": 1.2463587044656355e-08,
    "29": 1.1013629677730708e-08,
    "30": 1.1687642441683382e-08,
    "31": 1.1397713284329192e-08,
    "32": 1.1555927931648279e-08,
    "33": 1.3339402745926785e-08,
    "34": 1.0483812411227089e-08,
    "35": 1.4299102611645524e-08,
    "36": 1.187517782077471e-08,
    "37": 1.3834580623461596e-08,
    "38": 1.4363726512881696e-08,
    "X": 9.506483722244087e-09,
    "MT": 0,
}

_LindbladTohEtAl = stdpopsim.Citation(
        # Genome sequence, comparative analysis and haplotype structure of the
        # domestic dog.
        author="Lindblad-Toh et al.",
        year=2005,
        doi="https://doi.org/10.1038/nature04338")

_SkoglundEtAl = stdpopsim.Citation(
        # Ancient wolf genome reveals an early divergence of domestic dog
        # ancestors and admixture into high-latitude breeds.
        author="Skoglund et al.",
        year=2015,
        doi="https://doi.org/10.1016/j.cub.2015.04.019")

_FranzEtAl = stdpopsim.Citation(
        # Genomic and archaeological evidence suggest a dual origin of
        # domestic dogs.
        author="Franz et al.",
        year=2016,
        doi="https://doi.org/10.1126/science.aaf3161")

_CampbellEtAl = stdpopsim.Citation(
        # A Pedigree-Based Map of Recombination in the Domestic Dog Genome.
        author="Campbell et al.",
        year=2016,
        doi="https://doi.org/10.1534/g3.116.034678")

_chromosomes = []
for name, data in genome_data.data["chromosomes"].items():
    _chromosomes.append(stdpopsim.Chromosome(
        id=name, length=data["length"],
        synonyms=data["synonyms"],
        mutation_rate=4e-9,  # based on non-CpG sites only
        recombination_rate=_recombination_rate_data[name]))

_genome = stdpopsim.Genome(
    chromosomes=_chromosomes,
    assembly_name=genome_data.data["assembly_name"],
    assembly_accession=genome_data.data["assembly_accession"],
    mutation_rate_citations=[
        _SkoglundEtAl.because(stdpopsim.CiteReason.MUT_RATE),
        _FranzEtAl.because(stdpopsim.CiteReason.MUT_RATE),
        ],
    recombination_rate_citations=[
        _CampbellEtAl.because(stdpopsim.CiteReason.REC_RATE)
        ],
    assembly_citations=[
        _LindbladTohEtAl.because(stdpopsim.CiteReason.ASSEMBLY)
        ],
    )

_species = stdpopsim.Species(
    id="CanFam",
    name="Canis familiaris",
    common_name="Dog",
    genome=_genome,
    generation_time=3,
    generation_time_citations=[
        # Everyone uses 3 years because everyone else uses it.
        # It's likely higher, at least in wolves:
        # https://pubs.er.usgs.gov/publication/70187564
        ],
    population_size=13000,  # ancestral dog size
    population_size_citations=[
        _LindbladTohEtAl.because(stdpopsim.CiteReason.POP_SIZE)
        ],
    )

stdpopsim.register_species(_species)

_gm = stdpopsim.GeneticMap(
    species=_species,
    id="Campbell2016_CanFam3_1",
    description="Pedigree-based crossover map from 237 individuals",
    long_description="""
        Sex-averaged crossover frequency map based on 163,400 autosomal SNPs
        genotyped in a pedigree of 237 Labrador Retriever x Greyhound crosses.
        Genotypes were phased without respect to the pedigree, using SHAPEIT2,
        recombinations were called using duoHMM, and genetic distances were
        obtained using Haldane's map function.
        """,
    url="https://stdpopsim.s3-us-west-2.amazonaws.com/genetic_maps/"
        "CanFam/dog_genetic_maps.tar.gz",
    file_pattern="chr{id}_average_canFam3.1.txt",
    citations=[
            _CampbellEtAl.because(stdpopsim.CiteReason.GEN_MAP)
        ],
    )
_species.add_genetic_map(_gm)
