import collections

import stdpopsim
from . import genome_data

# De novo assembly of the cattle reference genome with single-molecule sequencing.
_RosenEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1093/gigascience/giaa021",
    year="2020",
    author="Rosen et al.",
)

# Frequency of mosaicism points towards mutation-prone early cleavage
# cell divisions in cattle.
_HarlandEtAl = stdpopsim.Citation(
    author="Harland et al.",
    year="2017",
    # BioRxiv preprint
    doi="https://doi.org/10.1101/079863",
)

# Cattle Sex-Specific Recombination and Genetic Control from a
# Large Pedigree Analysis.
_MaEtAl = stdpopsim.Citation(
    author="Ma et al.",
    year="2015",
    doi="https://doi.org/10.1371/journal.pgen.1005387",
)

# Inferring Demography from Runs of Homozygosity in Whole-Genome Sequence,
# with Correction for Sequence Errors.
_MacLeodEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1093/molbev/mst125",
    year="2013",
    author="MacLeod et al.",
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
    mutation_rate_citations=[
        _HarlandEtAl.because(stdpopsim.CiteReason.MUT_RATE),
    ],
    recombination_rate_citations=[_MaEtAl.because(stdpopsim.CiteReason.REC_RATE)],
    assembly_citations=[_RosenEtAl.because(stdpopsim.CiteReason.ASSEMBLY)],
)

_species = stdpopsim.Species(
    id="BosTau",
    ensembl_id="bos_taurus",
    name="Bos Taurus",
    common_name="Cattle",
    genome=_genome,
    generation_time=5,
    generation_time_citations=[_MacLeodEtAl.because(stdpopsim.CiteReason.GEN_TIME)],
    population_size=62000,
    population_size_citations=[_MacLeodEtAl.because(stdpopsim.CiteReason.POP_SIZE)],
)

stdpopsim.register_species(_species)
