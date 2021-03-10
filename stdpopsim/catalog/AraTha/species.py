import stdpopsim
from . import genome_data


_mean_recombination_rate = 200 / 124000 / 2 / 1e6

_recombination_rate_data = {str(j): _mean_recombination_rate for j in range(1, 6)}
_recombination_rate_data["Mt"] = 0
_recombination_rate_data["Pt"] = 0  # JK Is this correct??


# mutation rate from Ossowski 2010 Science
# recombination value from Huber et al 2014 MBE
# rho=200/Mb, assume Ne=124,000, rho=2*Ne*r
_chromosomes = []
for name, data in genome_data.data["chromosomes"].items():
    _chromosomes.append(
        stdpopsim.Chromosome(
            id=name,
            length=data["length"],
            synonyms=data["synonyms"],
            mutation_rate=7e-9,
            recombination_rate=_recombination_rate_data[name],
        )
    )

_genome = stdpopsim.Genome(
    chromosomes=_chromosomes,
    assembly_name=genome_data.data["assembly_name"],
    assembly_accession=genome_data.data["assembly_accession"],
    citations=[
        stdpopsim.Citation(
            author="Ossowski et al.",
            year=2010,
            doi="https://doi.org/10.1126/science.1180677",
            reasons={stdpopsim.CiteReason.MUT_RATE},
        ),
        stdpopsim.Citation(
            author="Huber et al.",
            year=2014,
            doi="https://doi.org/10.1093/molbev/msu247",
            reasons={stdpopsim.CiteReason.REC_RATE},
        ),
        stdpopsim.Citation(
            doi="https://doi.org/10.1093/nar/gkm965",
            year=2007,
            author="Swarbreck et al.",
            reasons={stdpopsim.CiteReason.ASSEMBLY},
        ),
    ],
)

_species = stdpopsim.Species(
    id="AraTha",
    ensembl_id="arabidopsis_thaliana",
    name="Arabidopsis thaliana",
    common_name="A. thaliana",
    genome=_genome,
    generation_time=1.0,
    population_size=10 ** 4,
    citations=[
        stdpopsim.Citation(
            doi="https://doi.org/10.1890/0012-9658(2002)083[1006:GTINSO]2.0.CO;2",
            year=2002,
            author="Donohue",
            reasons={stdpopsim.CiteReason.GEN_TIME},
        ),
        stdpopsim.Citation(
            doi="https://doi.org/10.1016/j.cell.2016.05.063",
            year=2016,
            author="1001GenomesConsortium",
            reasons={stdpopsim.CiteReason.POP_SIZE},
        ),
    ],
)

stdpopsim.register_species(_species)
