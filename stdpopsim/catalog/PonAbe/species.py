import stdpopsim
from . import genome_data


# Average recombination rates per chromosome taken from the Pongo abelii
# recombination map inferred in Nater et al (doi:10.1016/j.cub.2017.09.047),
# then lifted over from PonAbe2 to PonAbe3. PonAbe3 is nearly 400Mb shorter
# so this makes a noticable difference to the calculated averages.
# For chromosome X, we use the weighted genome-wide average.
_recombination_rate_data = {
    "1": 5.7209738520466636e-09,
    "2A": 6.746942995474612e-09,
    "2B": 6.088330307532069e-09,
    "3": 5.776152240522997e-09,
    "4": 6.061756714136327e-09,
    "5": 5.8783946576216106e-09,
    "6": 5.896045046261895e-09,
    "7": 6.82979118342936e-09,
    "8": 6.343129574251662e-09,
    "9": 7.047556965001838e-09,
    "10": 6.948119058906713e-09,
    "11": 5.791845481175716e-09,
    "12": 6.043159266363315e-09,
    "13": 6.4081697383466255e-09,
    "14": 6.110753350074597e-09,
    "15": 6.631500663980084e-09,
    "16": 6.7606578552415915e-09,
    "17": 8.01982524161532e-09,
    "18": 6.400648881415157e-09,
    "19": 9.096468143371218e-09,
    "20": 6.314242320359848e-09,
    "21": 7.588672333415156e-09,
    "22": 9.836284925758108e-09,
    "X": 6.402884398829772e-09,
    "MT": 0,
}

# High-resolution comparative analysis of great ape genomes
_kronenberg2018 = stdpopsim.Citation(
    author="Kronenberg et al.",
    year=2018,
    doi="https://doi.org/10.1126/science.aar6343",
    reasons={stdpopsim.CiteReason.ASSEMBLY},
)

# Comparative and demographic analysis of orang-utan genomes
_locke2011 = stdpopsim.Citation(
    author="Locke et al.",
    year=2011,
    doi="http://doi.org/10.1038/nature09687",
    reasons={stdpopsim.CiteReason.POP_SIZE},
)

# Morphometric, Behavioral, and Genomic Evidence for a New Orangutan Species
_nater2017 = stdpopsim.Citation(
    author="Nater et al.",
    year=2017,
    doi="https://doi.org/10.1016/j.cub.2017.09.047",
    reasons={stdpopsim.CiteReason.MUT_RATE, stdpopsim.CiteReason.REC_RATE},
)

# Orangutans: Geographic Variation in Behavioral Ecology and Conservation
# CHAPTER 5 Orangutan life history variation
_wich2008 = stdpopsim.Citation(
    author="Wich et al.",
    year=2008,
    doi="https://doi.org/10.1093/acprof:oso/9780199213276.003.0005",
    reasons={stdpopsim.CiteReason.GEN_TIME},
)


_chromosomes = []
for name, data in genome_data.data["chromosomes"].items():
    _chromosomes.append(
        stdpopsim.Chromosome(
            id=name,
            length=data["length"],
            synonyms=data["synonyms"],
            # Nater et al. 2017 used mu=1.5e-8 per generation, based on the
            # assumption that it's similar to humans and chimps.
            mutation_rate=1.5e-8,
            recombination_rate=_recombination_rate_data[name],
        )
    )

_genome = stdpopsim.Genome(
    chromosomes=_chromosomes,
    assembly_name=genome_data.data["assembly_name"],
    assembly_accession=genome_data.data["assembly_accession"],
    citations=[_kronenberg2018],
)
stdpopsim.utils.append_common_synonyms(_genome)

_species = stdpopsim.Species(
    id="PonAbe",
    ensembl_id="pongo_abelii",
    name="Pongo abelii",
    common_name="Sumatran orangutan",
    genome=_genome,
    # generation time used by Nater et al., citing Wich et al.
    generation_time=25,
    # Locke et al. inferred ancestral Ne
    population_size=1.79e4,
    citations=[_locke2011, _wich2008],
)

stdpopsim.register_species(_species)
