import stdpopsim

from . import genome_data

# Currently specifying a genome wide recombination rate
# reported by Renaut et al., 2013. This study, and other
# more recent studies, have produced good genetic maps 
# for sunflower (but it's not easy to find the data).
_overall_rate = 4e-9
_recombination_rate = {
    "1": _overall_rate,
    "2": _overall_rate,
    "3": _overall_rate,
    "4": _overall_rate,
    "5": _overall_rate,
    "6": _overall_rate,
    "7": _overall_rate,
    "8": _overall_rate,
    "9": _overall_rate,
    "10": _overall_rate,
    "11": _overall_rate,
    "12": _overall_rate,
    "13": _overall_rate,
    "14": _overall_rate,
    "15": _overall_rate,
    "16": _overall_rate,
    "17": _overall_rate,
}

# Average mutation rate reported by Andrew et al., 2013.
_overall_rate = 6.1e-9
_mutation_rate = {
    "1": _overall_rate,
    "2": _overall_rate,
    "3": _overall_rate,
    "4": _overall_rate,
    "5": _overall_rate,
    "6": _overall_rate,
    "7": _overall_rate,
    "8": _overall_rate,
    "9": _overall_rate,
    "10": _overall_rate,
    "11": _overall_rate,
    "12": _overall_rate,
    "13": _overall_rate,
    "14": _overall_rate,
    "15": _overall_rate,
    "16": _overall_rate,
    "17": _overall_rate,
}

#
_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
    citations=[
        stdpopsim.Citation(
            doi="https://doi.org/10.1038/nature22380",
            year=2017,
            author="Badouin et al.",
            reasons={stdpopsim.CiteReason.ASSEMBLY},
        ),
        stdpopsim.Citation(
            doi="https://doi.org/10.1038/ncomms2833",
            year=2013,
            author="Renaut et al.",
            reasons={stdpopsim.CiteReason.REC_RATE},
        ),
        stdpopsim.Citation(
            doi="https://doi.org/10.1111/mec.12038",
            year=2013,
            author="Andrew et al.",
            reasons={stdpopsim.CiteReason.MUT_RATE},
        ),
    ],
)

_species = stdpopsim.Species(
    id="HelAnn",
    ensembl_id="helianthus_annuus",
    name="Helianthus annuus",
    common_name="Helianthus annuus",
    genome=_genome,
    generation_time=1.0,
    population_size=673968,
    citations=[
        stdpopsim.Citation(
            doi="https://doi.org/10.1093/molbev/msq270",
            year=2010,
            author="Strasburg et al.",
            reasons={stdpopsim.CiteReason.POP_SIZE},
        ),
        stdpopsim.Citation(
            doi="https://doi.org/10.1111/j.1558-5646.2008.00415.x",
            year=2008,
            author="Strasburg and Rieseberg",
            reasons={stdpopsim.CiteReason.GEN_TIME},
        ),
    ],
)

stdpopsim.register_species(_species)

