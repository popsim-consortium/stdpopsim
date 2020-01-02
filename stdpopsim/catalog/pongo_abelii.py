"""
Orangutans.
"""
import math

import msprime

import stdpopsim

import logging

logger = logging.getLogger(__name__)

###########################################################
#
# Genome definition
#
###########################################################

# List of chromosomes.

# length information can be found here
# <https://www.ncbi.nlm.nih.gov/genome/325>
# Orangutan July 201 (WUGSC 2.0.2/ponAbe2) Browser Sequences
# Average recombination rates per chromosome taken from the Pongo abelii
# recombination map inferred in Nater et al (doi:10.1016/j.cub.2017.09.047)
# For chromosome X, we use the rate assumed in Locke et al.
_chromosome_data = """\
chr1     229942017   5.19e-9
chr2a    113028656   5.43e-9
chr2b    135000294   5.38e-9
chr3     202140232   5.36e-9
chr4     198332218   5.43e-9
chr5     183952662   5.20e-9
chr6     174210431   5.22e-9
chr7     157549271   5.73e-9
chr8     153482349   5.67e-9
chr9     135191526   5.34e-9
chr10    133410057   5.91e-9
chr11    132107971   5.29e-9
chr12    136387465   5.44e-9
chr13    117095149   4.91e-9
chr14    108868599   4.70e-9
chr15    99152023    4.82e-9
chr16    77800216    6.12e-9
chr17    73212453    7.26e-9
chr18    94050890    4.57e-9
chr19    60714840    7.56e-9
chr20    62736349    5.83e-9
chr21    48394510    4.98e-9
chr22    46535552    6.03e-9
chrX     156195299   9.50e-9
"""

_locke2011 = stdpopsim.Citation(
    author="Locke et al.",
    year=2011,
    doi="http://doi.org/10.1038/nature09687"
)

_nater2017 = stdpopsim.Citation(
    author="Nater et al.",
    year=2017,
    doi="https://doi.org/10.1016/j.cub.2017.09.047"
)

_chromosomes = []
for line in _chromosome_data.splitlines():
    name, length, mean_rr = line.split()[:3]
    _chromosomes.append(stdpopsim.Chromosome(
        id=name, length=int(length),
        mutation_rate=1.5e-8,
        recombination_rate=float(mean_rr)))

_genome = stdpopsim.Genome(
    chromosomes=_chromosomes,
    mutation_rate_citations=[
        _nater2017.because(stdpopsim.CiteReason.MUT_RATE)])

_species = stdpopsim.Species(
    id="PonAbe",
    name="Pongo abelii",
    common_name="Sumatran orangutan",
    genome=_genome,
    generation_time=20,
    generation_time_citations=[
        _locke2011.because(stdpopsim.CiteReason.GEN_TIME)],
    population_size=1.79e4,
    population_size_citations=[
        _locke2011.because(stdpopsim.CiteReason.POP_SIZE)])

stdpopsim.register_species(_species)


###########################################################
#
# Genetic maps
#
###########################################################

# There are two genetic maps available for Orangutan species: one for Pongo
# abelii (Sumatran orangutan) and one for Pongo pygmaeus (Bornean orangutan).
# Both recombination maps were inferred using LDhat in Nater et al. (2017),
# doi: 10.1016/j.cub.2017.09.047. Both recombination maps are mapped to PonAbe2.
# Recombination maps from Nater et al. were converted from rho/kbp to cM using
# Watterson's estimator of theta to estimate Ne = 41,000 (Sumatra) and
# Ne = 27,000 (Borneo). See supporting information in Nater et al. for details.


_gm_pa = stdpopsim.GeneticMap(
    species=_species,
    id="NaterPA_PonAbe2",
    description="From Nater et al. (2017) for Pongo abelii",
    long_description="""
        This genetic map is from the Nater et al. (2017) study, inferred using
        LDhat from n=15 whole-genome sequenced Sumatran orangutan individuals.
        See https://doi.org/10.1016/j.cub.2017.09.047 for more details.
        """,
    url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/genetic_maps/PonAbe/NaterPA_PonAbe2.tar.gz"),  # NOQA
    file_pattern="Nater_et_al_PA_{id}_PonAbe2.txt",
    citations=[
        stdpopsim.Citation(
            doi="https://doi.org/10.1016/j.cub.2017.09.047",
            year=2017,
            author="Nater et al.",
            reasons={stdpopsim.CiteReason.GEN_MAP}),
    ]
)
_species.add_genetic_map(_gm_pa)

_gm_pp = stdpopsim.GeneticMap(
    species=_species,
    id="NaterPP_PonAbe2",
    description="From Nater et al. (2017) for Pongo pygmaeus",
    long_description="""
        This genetic map is from the Nater et al. (2017) study, inferred using
        LDhat from n=20 whole-genome sequenced Bornean orangutan individuals.
        See https://doi.org/10.1016/j.cub.2017.09.047 for more details.
        """,
    url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/genetic_maps/PonPyg/NaterPP_PonAbe2.tar.gz"),  # NOQA
    file_pattern="Nater_et_al_PP_{id}_PonAbe2.txt",
    citations=[
        stdpopsim.Citation(
            doi="https://doi.org/10.1016/j.cub.2017.09.047",
            year=2017,
            author="Nater et al.",
            reasons={stdpopsim.CiteReason.GEN_MAP})]
)
_species.add_genetic_map(_gm_pp)


###########################################################
#
# Demographic models
#
###########################################################


def _orangutan():
    id = "TwoSpecies_2L11"
    description = "Two population orangutan model"
    long_description = """
        The two orang-utan species, Sumatran (Pongo abelii) and Bornean (Pongo
        pygmaeus) inferred from the joint-site frequency spectrum with ten
        individuals from each population. This model is an isolation-with-
        migration model, with exponential growth or decay in each population
        after the split. The Sumatran population grows in size, while the
        Bornean population slightly declines.
    """

    citations = [_locke2011.because(stdpopsim.CiteReason.DEM_MODEL)]

    populations = [
        stdpopsim.Population(
            "Bornean", "Pongo pygmaeus (Bornean) population"),
        stdpopsim.Population(
            "Sumatran", "Pongo abelii (Sumatran) population")
    ]

    # Parameters from paper:
    # ancestral size, before split
    Na = 17934

    # time of split
    T_split_years = 403149
    # get split time in units of generations
    generation_time = _species.generation_time
    T_split = T_split_years / generation_time

    # proportion of ancestral pop to branch B
    s = 0.592

    # sizes at split
    Na_B = 17934*s
    Na_S = 17934*(1-s)

    # present sizes
    N_B = 8805
    N_S = 37661

    # get growth rates
    r_B = -1*math.log(Na_B/N_B)/T_split
    r_S = -1*math.log(Na_S/N_S)/T_split

    # migration rates
    m_S_B = 0.395 / 2 / Na
    m_B_S = 0.239 / 2 / Na

    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=B and 1=S
    # initially.

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        citations=citations,
        populations=populations,
        generation_time=generation_time,
        population_configurations=[
            msprime.PopulationConfiguration(initial_size=N_B, growth_rate=r_B,
                                            metadata=populations[0].asdict()),  # NOQA
            msprime.PopulationConfiguration(initial_size=N_S, growth_rate=r_S,
                                            metadata=populations[1].asdict())  # NOQA
        ],
        migration_matrix=[
            [0, m_B_S],  # NOQA
            [m_S_B,       0],  # NOQA
        ],
        demographic_events=[
            # Merge populations and turn off migration, change to size Na
            msprime.MassMigration(
                time=T_split, source=1, destination=0, proportion=1.0),
            msprime.MigrationRateChange(time=T_split, rate=0),
            msprime.PopulationParametersChange(
                time=T_split, initial_size=Na, growth_rate=0, population_id=0),
        ],
        )


_species.add_demographic_model(_orangutan())
