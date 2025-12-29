import stdpopsim
from . import genome_data


# genome-wide recombination rate from Huber et al 2014 MBE
# associated with all recombining chromosomes
_rho = 200 / 1e6  # 200/Mb
_Ne = 124000
_mean_recombination_rate = _rho / (2 * _Ne)
_recombination_rate = {str(j): _mean_recombination_rate for j in range(1, 6)}
_recombination_rate["Mt"] = 0
_recombination_rate["Pt"] = 0  # JK Is this correct??

# genome-wide average mutation rate from Ossowski 2010 Science
# associated with all chromosomes
_mean_mutation_rate = 7e-9
_mutation_rate = {str(j): _mean_mutation_rate for j in range(1, 6)}
_mutation_rate["Mt"] = _mean_mutation_rate
_mutation_rate["Pt"] = _mean_mutation_rate

# species ploidy and chromosome-specific ploidy
_species_ploidy = 2
_ploidy = {str(j): _species_ploidy for j in range(1, 6)}
_ploidy["Mt"] = 1
_ploidy["Pt"] = 1

_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
    ploidy=_ploidy,
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

# add common chromosome synonyms to _genome object
stdpopsim.utils.append_common_synonyms(_genome)

_species = stdpopsim.Species(
    id="AraTha",
    ensembl_id="arabidopsis_thaliana",
    name="Arabidopsis thaliana",
    common_name="A. thaliana",
    separate_sexes=False,
    genome=_genome,
    generation_time=1.0,
    population_size=10**4,
    ploidy=_species_ploidy,
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
