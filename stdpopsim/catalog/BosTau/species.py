import collections
import stdpopsim

from . import genome_data

# De novo assembly of the cattle reference genome with single-molecule sequencing.
_RosenEtAl = stdpopsim.Citation(
    author="Rosen et al.",
    year=2020,
    doi="https://doi.org/10.1093/gigascience/giaa021",
    reasons={stdpopsim.CiteReason.ASSEMBLY},
)

# Frequency of mosaicism points towards mutation-prone early cleavage
# cell divisions in cattle.
_HarlandEtAl = stdpopsim.Citation(
    author="Harland et al.",
    year=2017,
    doi="https://doi.org/10.1101/079863",
    reasons={stdpopsim.CiteReason.MUT_RATE},
)

# Cattle Sex-Specific Recombination and Genetic Control from a
# Large Pedigree Analysis.
_MaEtAl = stdpopsim.Citation(
    author="Ma et al.",
    year=2015,
    doi="https://doi.org/10.1371/journal.pgen.1005387",
    reasons={stdpopsim.CiteReason.REC_RATE},
)

# Inferring Demography from Runs of Homozygosity in Whole-Genome Sequence,
# with Correction for Sequence Errors.
_MacLeodEtAl = stdpopsim.Citation(
    author="MacLeod et al.",
    year=2013,
    doi="https://doi.org/10.1093/molbev/mst125",
    reasons={
        stdpopsim.CiteReason.GEN_TIME,
        stdpopsim.CiteReason.POP_SIZE,
    },
)

# Recombination rate has been derived from dairy cattle crossovers
# per meiosis, by taking the average between females and males and then
# dividing by the whole genome length (equal to the sum of chromosome
# lengths used above).
# From Ma et al. (2015), 25.5 crossovers per meiosis in males and
# 23.2 crossovers per meiosis in females, gives an average of 24.35
# crossovers per meiosis. The sum of chromosome lengths is 2628394923 bp.
# 24.35 / 2628394923 = 9.26e-9 per bp per generation.
_genome_wide_recombination_rate = 9.26e-9

_recombination_rate_data = collections.defaultdict(
    lambda: _genome_wide_recombination_rate
)
# Set some exceptions for non-recombining chrs.
_recombination_rate_data["MT"] = 0

_chromosomes = []
for name, data in genome_data.data["chromosomes"].items():
    _chromosomes.append(
        stdpopsim.Chromosome(
            id=name,
            length=data["length"],
            synonyms=data["synonyms"],
            # Harland et al. (2017), sex-averaged estimate per bp per generation.
            mutation_rate=1.2e-8,
            recombination_rate=_recombination_rate_data[name],
        )
    )

_genome = stdpopsim.Genome(
    chromosomes=_chromosomes,
    citations=[
        _RosenEtAl,
        _HarlandEtAl,
        _MaEtAl,
    ],
)
stdpopsim.utils.append_common_synonyms(_genome)

_species = stdpopsim.Species(
    id="BosTau",
    ensembl_id="bos_taurus",
    name="Bos taurus",
    common_name="Cattle",
    genome=_genome,
    generation_time=5,
    population_size=90,  # most recent Ne in _MacLeodEtAl
    citations=[_MacLeodEtAl],
)

stdpopsim.register_species(_species)
