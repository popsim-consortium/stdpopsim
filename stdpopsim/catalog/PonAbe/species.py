import stdpopsim
from . import genome_data


# Average recombination rates per chromosome taken from the Pongo abelii
# recombination map inferred in Nater et al (doi:10.1016/j.cub.2017.09.047)
# For chromosome X, we use the rate assumed in Locke et al.
_recombination_rate_data = {
    "1": 5.19e-9,
    "2a": 5.43e-9,
    "2b": 5.38e-9,
    "3": 5.36e-9,
    "4": 5.43e-9,
    "5": 5.20e-9,
    "6": 5.22e-9,
    "7": 5.73e-9,
    "8": 5.67e-9,
    "9": 5.34e-9,
    "10": 5.91e-9,
    "11": 5.29e-9,
    "12": 5.44e-9,
    "13": 4.91e-9,
    "14": 4.70e-9,
    "15": 4.82e-9,
    "16": 6.12e-9,
    "17": 7.26e-9,
    "18": 4.57e-9,
    "19": 7.56e-9,
    "20": 5.83e-9,
    "21": 4.98e-9,
    "22": 6.03e-9,
    "X": 9.50e-9,
    "MT": 0,
}

_locke2011 = stdpopsim.Citation(
    author="Locke et al.", year=2011, doi="http://doi.org/10.1038/nature09687"
)

_nater2017 = stdpopsim.Citation(
    author="Nater et al.", year=2017, doi="https://doi.org/10.1016/j.cub.2017.09.047"
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
    mutation_rate_citations=[_nater2017.because(stdpopsim.CiteReason.MUT_RATE)],
)

_species = stdpopsim.Species(
    id="PonAbe",
    ensembl_id="pongo_abelii",
    name="Pongo abelii",
    common_name="Sumatran orangutan",
    genome=_genome,
    # generation time used by Locke et al. without further citation
    generation_time=20,
    generation_time_citations=[_locke2011.because(stdpopsim.CiteReason.GEN_TIME)],
    # Locke et al. inferred ancestral Ne
    population_size=1.79e4,
    population_size_citations=[_locke2011.because(stdpopsim.CiteReason.POP_SIZE)],
)

stdpopsim.register_species(_species)
