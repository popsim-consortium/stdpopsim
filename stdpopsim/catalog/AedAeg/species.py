import stdpopsim

from . import genome_data

# These are in Table 1 of Juneja et al:
_recombination_rate = {"1": 0.306e-8, "2": 0.249e-8, "3": 0.291e-8, "MT": 0}

# Generic and chromosome-specific ploidy
_species_ploidy = 2
_ploidy = {
    "1": _species_ploidy,
    "2": _species_ploidy,
    "3": _species_ploidy,
    "MT": 1,
}

_JunejaEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1371/journal.pntd.0002652",
    year=2014,
    author="Juneja et al.",
    reasons={stdpopsim.CiteReason.REC_RATE},
)


_CrawfordEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1186/s12915-017-0351-0",
    year=2017,
    author="Crawford et al.",
    reasons={
        stdpopsim.CiteReason.GEN_TIME,
        stdpopsim.CiteReason.POP_SIZE,
        stdpopsim.CiteReason.MUT_RATE,
    },
)

_KeightleyEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1101/gr.091231.109",
    year=2009,
    author="Keightley et al.",
    reasons={
        stdpopsim.CiteReason.MUT_RATE,
    },
)


# This rate is for DroMel, as used in Crawford
# TODO: mt should probably be lower
_overall_rate = 3.5e-9  # per generation, from Keightley et al 2009
_mutation_rate = {
    "1": _overall_rate,
    "2": _overall_rate,
    "3": _overall_rate,
    "MT": _overall_rate,
}

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
    ploidy=_ploidy,
    citations=[
        _NeneEtAl,
        _JunejaEtAl,
        _CrawfordEtAl,
        _KeightleyEtAl,
    ],
)
stdpopsim.utils.append_common_synonyms(_genome)

_species = stdpopsim.Species(
    id="AedAeg",
    ensembl_id="aedes_aegypti_lvpagwg",
    name="Aedes aegypti",
    common_name="Yellow fever mosquito",
    separate_sexes=True,
    genome=_genome,
    generation_time=1 / 15,
    ploidy=_species_ploidy,
    # the estimated population size today the modern Senegal forest population
    population_size=1e6,
    citations=[_CrawfordEtAl],
)

stdpopsim.register_species(_species)
