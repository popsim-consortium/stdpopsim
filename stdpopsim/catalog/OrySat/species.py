import collections
import stdpopsim

from . import genome_data

_KawaharaEtAl = stdpopsim.Citation(
    author="Kawahara et al.",
    year=2013,
    doi="https://doi.org/10.1186/1939-8433-6-4",
    reasons={stdpopsim.CiteReason.ASSEMBLY},
)

# multiple first authors
_DuEtAl = stdpopsim.Citation(
    author="Du et al.",
    year=2017,
    doi="https://doi.org/10.1038/ncomms15324",
    reasons={stdpopsim.CiteReason.ASSEMBLY},
)

_YangEtAl = stdpopsim.Citation(
    author="Yang et al.",
    year=2015,
    doi="https://doi.org/10.1038/nature14649",
    reasons={stdpopsim.CiteReason.MUT_RATE},
)

_SiEtAl = stdpopsim.Citation(
    author="Si et al.",
    year=2015,
    doi="https://doi.org/10.1111/nph.13319",
    reasons={stdpopsim.CiteReason.REC_RATE},
)

_HuangEtAl = stdpopsim.Citation(
    author="Huang et al.",
    year=2012,
    doi="https://doi.org/10.1038/nature11532",
    reasons={stdpopsim.CiteReason.POP_SIZE},
)

# Ne is reported for Oryza sativa japonica
# It is calculated from pi reported in
# Huang et al. 2012
# Ne = 0.0006 / (4 * 3.2e-9)
# No sensible reason to ever model the
# demography of the japonica and indica
# subspecies together, so using the mean
# of japonica and indica would make less
# sense

_CaicedoEtAl = stdpopsim.Citation(
    author="Caicedo et al.",
    year=2007,
    doi="https://doi.org/10.1371/journal.pgen.0030163",
    reasons={stdpopsim.CiteReason.GEN_TIME, stdpopsim.CiteReason.REC_RATE},
)

# Recombination rate was calculated from sequences of rice
# F2 plants as 4.53 cM/Mb, which translates to 4.53e-8 per bp per generation
# The effective recombination rate, adjusting for the high observed
# selfing rate in rice, can be calculated using the formula
# r(1 − σ/[2 − σ]), where σ is the selfing rate,
# as applied in Caicedo et al.
# Using a selfing rate of 0.99 as in Caicedo et al.:
_genome_wide_recombination_rate = 8.97e-10

_recombination_rate_data = collections.defaultdict(
    lambda: _genome_wide_recombination_rate
)
# Set some exceptions for non-recombining chrs.
_recombination_rate_data["Mt"] = 0
_recombination_rate_data["Pt"] = 0

# Generic and chromosome-specific ploidy
_species_ploidy = 2
_ploidy = {
    "1": _species_ploidy,
    "2": _species_ploidy,
    "3": _species_ploidy,
    "4": _species_ploidy,
    "5": _species_ploidy,
    "6": _species_ploidy,
    "7": _species_ploidy,
    "8": _species_ploidy,
    "9": _species_ploidy,
    "10": _species_ploidy,
    "11": _species_ploidy,
    "12": _species_ploidy,
    "Mt": 1,
    "Pt": 1,
}

_chromosomes = []
for name, data in genome_data.data["chromosomes"].items():
    _chromosomes.append(
        stdpopsim.Chromosome(
            id=name,
            length=data["length"],
            synonyms=data["synonyms"],
            mutation_rate=3.2e-9,
            recombination_rate=_recombination_rate_data[name],
            ploidy=_ploidy[name],
        )
    )

_genome = stdpopsim.Genome(
    chromosomes=_chromosomes,
    citations=[
        _KawaharaEtAl,
        _DuEtAl,
        _YangEtAl,
        _SiEtAl,
    ],
)
stdpopsim.utils.append_common_synonyms(_genome)

_species = stdpopsim.Species(
    id="OrySat",
    ensembl_id="oryza_sativa",
    name="Oryza sativa",
    common_name="Asian rice",
    separate_sexes=False,
    genome=_genome,
    generation_time=1,
    population_size=46875,
    ploidy=_species_ploidy,
    citations=[_CaicedoEtAl],
)

stdpopsim.register_species(_species)
