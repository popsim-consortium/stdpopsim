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

# Mean chromosomal rates, calculated from the Comeron 2012 map.
# Chromosome 4 isn't in this map, so the weighted mean of 2L, 2R, 3L and 3R
# was used instead.
_recombination_rate_data = {
    "2L": 2.4125016027908946e-08,
    "2R": 2.2366522822806982e-08,
    "3L": 1.7985794693631893e-08,
    "3R": 1.7165556232922828e-08,
    "4": 2.0085234464525437e-08,
    "X": 2.9151053903465754e-08,
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
