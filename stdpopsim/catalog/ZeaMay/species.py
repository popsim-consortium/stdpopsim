import stdpopsim

from . import genome_data

# Mean per-chromosome recombination rates (per bp per generation), computed
# as the total genetic length of each chromosome divided by its physical
# length in Zm-B73-REFERENCE-NAM-5.0. Total genetic lengths (cM) are those
# of the Ogut et al. (2015) NAM consensus map, which is also the
# per-chromosome length that the FineMap genetic map is calibrated to
# (Ross-Ibarra 2026):
#   chr1 210.4, chr2 161.2, chr3 163.4, chr4 151.8, chr5 157.0,
#   chr6 111.4, chr7 138.4, chr8 137.4, chr9 131.2, chr10 113.0
_recombination_rate = {
    "1": 6.8211e-09,
    "2": 6.6154e-09,
    "3": 6.8650e-09,
    "4": 6.0640e-09,
    "5": 6.9361e-09,
    "6": 6.1426e-09,
    "7": 7.4485e-09,
    "8": 7.5324e-09,
    "9": 8.0488e-09,
    "10": 7.4130e-09,
}

# Genome-wide mutation rate per bp per generation, estimated from a
# mutation accumulation experiment (Stitzer et al., unpublished; presented
# in the talk cited below). The same rate was assumed when inferring the
# Magalang et al. demographic model.
_overall_mutation_rate = 7.3e-8
_mutation_rate = {c: _overall_mutation_rate for c in genome_data.data["chromosomes"]}

_species_ploidy = 2
_ploidy = {c: _species_ploidy for c in genome_data.data["chromosomes"]}

# Unpublished results presented in a recorded talk; to be replaced with
# preprint/publication DOIs when available.
_talk_url = "https://www.youtube.com/watch?v=9K2iWaLZQLY"

_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
    ploidy=_ploidy,
    citations=[
        stdpopsim.Citation(
            author="Hufford et al.",
            year=2021,
            doi="https://doi.org/10.1126/science.abg5289",
            reasons={stdpopsim.CiteReason.ASSEMBLY},
        ),
        stdpopsim.Citation(
            author="Ogut et al.",
            year=2015,
            doi="https://doi.org/10.1038/hdy.2014.123",
            reasons={stdpopsim.CiteReason.REC_RATE},
        ),
        stdpopsim.Citation(
            author="Stitzer et al.",
            year=2026,
            doi=_talk_url,
            reasons={stdpopsim.CiteReason.MUT_RATE},
        ),
    ],
)

# Population size is the time-weighted harmonic mean of the population
# sizes in the Magalang et al. demographic model (see
# demographic_models.py) between 10 and 1,000,000 generations ago.
# Maize is an annual plant, so generation time is 1 year.
_species = stdpopsim.Species(
    id="ZeaMay",
    ensembl_id="zea_mays",
    name="Zea mays",
    common_name="Maize",
    separate_sexes=False,
    genome=_genome,
    ploidy=_species_ploidy,
    generation_time=1,
    population_size=97815,
    citations=[
        stdpopsim.Citation(
            author="Magalang et al.",
            year=2026,
            doi=_talk_url,
            reasons={
                stdpopsim.CiteReason.POP_SIZE,
                stdpopsim.CiteReason.GEN_TIME,
            },
        ),
    ],
)

stdpopsim.register_species(_species)
