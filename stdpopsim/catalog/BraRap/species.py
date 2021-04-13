import stdpopsim

from . import genome_data

# I still need an estimate...
_mean_recombination_rate = 0

_recombination_rate = {
    "A01": _mean_recombination_rate,
    "A02": _mean_recombination_rate,
    "A03": _mean_recombination_rate,
    "A04": _mean_recombination_rate,
    "A05": _mean_recombination_rate,
    "A06": _mean_recombination_rate,
    "A07": _mean_recombination_rate,
    "A08": _mean_recombination_rate,
    "A09": _mean_recombination_rate,
    "A10": _mean_recombination_rate,
}

_mean_mutation_rate = 9e-9

_mutation_rate = {
    "A01": _mean_mutation_rate,
    "A02": _mean_mutation_rate,
    "A03": _mean_mutation_rate,
    "A04": _mean_mutation_rate,
    "A05": _mean_mutation_rate,
    "A06": _mean_mutation_rate,
    "A07": _mean_mutation_rate,
    "A08": _mean_mutation_rate,
    "A09": _mean_mutation_rate,
    "A10": _mean_mutation_rate,
}

_recombination_rate_citation = stdpopsim.Citation(
    author="", year=-1, doi="", reasons={stdpopsim.CiteReason.REC_RATE}
)

_mutation_rate_citation = stdpopsim.Citation(
    doi="https://doi.org/10.1105/tpc.114.126391",
    author="Kagale et al.",
    year=2014,
    reasons={stdpopsim.CiteReason.MUT_RATE},
)

_assembly_citation = stdpopsim.Citation(
    author="Wang et al.",
    year=2011,
    doi="https://doi.org/10.1038/ng.919",
    reasons={stdpopsim.CiteReason.ASSEMBLY},
)

_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
    citations=[
        _assembly_citation,
        _recombination_rate_citation,
        _mutation_rate_citation,
    ],
)

_species = stdpopsim.Species(
    id="BraRap",
    ensembl_id="brassica_rapa",
    name="Brassica rapa",
    common_name="Brassica rapa",
    genome=_genome,
    # Brassica rapa is an anual crop
    generation_time=1,
    # Population size inferred using dadi from the SFS in Qi et al
    population_size=40e3,
    citations=[
        stdpopsim.Citation(
            author="Qi et al.",
            year=2017,
            doi="https://doi.org/10.1111/mec.14131",
            reasons={stdpopsim.CiteReason.POP_SIZE},
        ),
    ],
)

stdpopsim.register_species(_species)
