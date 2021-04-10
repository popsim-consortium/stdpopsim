import stdpopsim

from . import genome_data


_Ag1000G = stdpopsim.Citation(
    doi="https://doi.org/10.1038/nature24995",
    year=2017,
    author="Ag1000G Consortium",
    reasons={
        stdpopsim.CiteReason.MUT_RATE,
        stdpopsim.CiteReason.GEN_TIME,
        stdpopsim.CiteReason.POP_SIZE,
    },
)

_PombiEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.4269/ajtmh.2006.75.901",
    year=2006,
    author="Pombi et al.",
    reasons={stdpopsim.CiteReason.REC_RATE},
)

_ZhengEtAl = stdpopsim.Citation(
    # there's no DOI?!?!
    doi="https://www.genetics.org/content/143/2/941",
    year=1996,
    author="Zheng et al.",
    reasons={stdpopsim.CiteReason.REC_RATE},
)

_SharakhovaEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1186/gb-2007-8-1-r5",
    year=2006,
    author="Sharakhova et al.",
    reasons={stdpopsim.CiteReason.ASSEMBLY},
)

_KeightleyEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1101/gr.091231.109",
    year=2009,
    author="Keightley et al.",
    reasons={stdpopsim.CiteReason.MUT_RATE},
)

# values are in Pombi et al, but except for the X,
# they refer to Zheng et al for the data
# and have no data for 2L so we set equal to 2R
# the value for X is the average from Table 1 in Zheng
_recombination_rate = {
    "2L": 1.3e-8,
    "2R": 1.3e-8,
    "3L": 1.3e-8,
    "3R": 1.6e-8,
    "X": 2.04e-8,
    "Mt": 0,
}


# the rate for DroMel as inferred by Keightley et al 2009,
# which was used for the stairwayplot demographic inference
_overall_rate = 3.5e-9
_mutation_rate = {
    "2L": _overall_rate,
    "2R": _overall_rate,
    "3L": _overall_rate,
    "3R": _overall_rate,
    "X": _overall_rate,
    "Mt": _overall_rate,
}

_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
    citations=[
        _SharakhovaEtAl,
        _Ag1000G,
        _PombiEtAl,
        _KeightleyEtAl,
        _ZhengEtAl,
    ],
)
stdpopsim.utils.append_common_synonyms(_genome)

_species = stdpopsim.Species(
    id="AnoGam",
    ensembl_id="anopheles_gambiae",
    name="Anopheles gambiae",
    common_name="Anopheles gambiae",
    genome=_genome,
    generation_time=1 / 11,
    # based on theta = 4 Ne u in Gabon population, rounded
    population_size=1e6,
    citations=[_Ag1000G],
)

stdpopsim.register_species(_species)
