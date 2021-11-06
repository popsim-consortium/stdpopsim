import stdpopsim
from . import genome_data


_LiAndStephan = stdpopsim.Citation(
    author="Li et al.",
    year=2006,
    doi="https://doi.org/10.1371/journal.pgen.0020166",
    reasons={stdpopsim.CiteReason.GEN_TIME, stdpopsim.CiteReason.POP_SIZE},
)

_SchriderEtAl = stdpopsim.Citation(
    author="Schrider et al.",
    year=2013,
    doi="https://doi.org/10.1534/genetics.113.151670",
)

_DosSantosEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1093/nar/gku1099",
    year=2015,
    author="dos Santos et al.",
    reasons={stdpopsim.CiteReason.ASSEMBLY},
)

_HoskinsEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1101/gr.185579.114",
    year=2015,
    author="Hoskins et al.",
    reasons={stdpopsim.CiteReason.ASSEMBLY},
)

_ComeronEtAl = stdpopsim.Citation(
    author="Comeron et al",
    doi="https://doi.org/10.1371/journal.pgen.1002905",
    year=2012,
)

# Mean chromosomal rates, calculated by taking the
# average of all rates in the
# "ComeronCrossoverV2_dm6" genetic map, weighted by
# their respective interval lengths.
# Chromosome 4 isn't in this map, so the average of
# 2L, 2R, 3L and 3R weighted by their respective
# chromosome lengths was used instead.
_recombination_rate_data = {
    "2L": 2.40462600791e-08,
    "2R": 2.23458641776e-08,
    "3L": 1.79660308862e-08,
    "3R": 1.71642045777e-08,
    "4": 0,
    "X": 2.89650687913e-08,
    "Y": 0,
    "mitochondrion_genome": 0,
}


_chromosomes = []
for name, data in genome_data.data["chromosomes"].items():
    _chromosomes.append(
        stdpopsim.Chromosome(
            id=name,
            length=data["length"],
            synonyms=data["synonyms"],
            mutation_rate=5.49e-9,  # _SchriderEtAl de novo mutation rate
            recombination_rate=_recombination_rate_data[name],
        )
    )

_genome = stdpopsim.Genome(
    chromosomes=_chromosomes,
    assembly_name=genome_data.data["assembly_name"],
    assembly_accession=genome_data.data["assembly_accession"],
    citations=[
        _SchriderEtAl.because(stdpopsim.CiteReason.MUT_RATE),
        _DosSantosEtAl,
        _HoskinsEtAl,
        _ComeronEtAl.because(stdpopsim.CiteReason.REC_RATE),
    ],
)
stdpopsim.utils.append_common_synonyms(_genome)

_species = stdpopsim.Species(
    id="DroMel",
    ensembl_id="drosophila_melanogaster",
    name="Drosophila melanogaster",
    common_name="D. melanogaster",
    genome=_genome,
    generation_time=0.1,
    # Population size is the older of two population sizes estimated by
    # Li and Stephan in a two-epoch model of African populations.
    # N_A0 is given as 8.603e6, and N_A1 (used here) is 5 times smaller.
    population_size=1720600,
    citations=[_LiAndStephan],
)

stdpopsim.register_species(_species)
