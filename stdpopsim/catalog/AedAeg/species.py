import stdpopsim

from . import genome_data

# These are in Table 1 of Juneja et al:
_recombination_rate = {"1": 0.306, "2": 0.249, "3": 0.291, "MT": 0}

_JunejaEtAl = stdpopsim.Citation(
    doi="https://doi:10.1371/journal.pntd.0002652",
    year=2014,
    author="Juneja et al.",
    reasons={stdpopsim.CiteReason.REC_RATE},
)

# This rate is for Drosophila (same as used in the A. gambiae tutorial).
_overall_rate = 5.49e-9
_mutation_rate = {
    "1": _overall_rate,
    "2": _overall_rate,
    "3": _overall_rate,
    "MT": _overall_rate,
}

_Ag1000G = stdpopsim.Citation(
    doi="https://doi.org/10.1038/nature24995",
    year=2017,
    author="Ag1000G Consortium",
    reasons={stdpopsim.CiteReason.MUT_RATE},
)

_NeneEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1126/science.1138878",
    year=2007,
    author="Nene et al.",
    reasons={stdpopsim.CiteReason.ASSEMBLY},
)


_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
    citations=[
        _NeneEtAl,
        _Ag1000G,
        _JunejaEtAl,
    ],
)


_CrawfordEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1186/s12915-017-0351-0",
    year=2017,
    author="Crawford et al.",
    reasons={stdpopsim.CiteReason.GEN_TIME, stdpopsim.CiteReason.POP_SIZE},
)

_species = stdpopsim.Species(
    id="AedAeg",
    ensembl_id="aedes_aegypti_lvpagwg",
    name="Aedes aegypti",
    common_name="Aedes aegypti (LVP_AGWG)",
    genome=_genome,
    generation_time=1 / 10,
    population_size=1e6,
    citations=[_CrawfordEtAl],
)

stdpopsim.register_species(_species)
