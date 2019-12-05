"""
Orangutans.
"""
import math

import msprime

import stdpopsim.models as models

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
# <http://https://www.ncbi.nlm.nih.gov/genome/325>
# Orangutan Jan. 2018 (Susie_PABv2/ponAbe3) Browser Sequences
# Note: no chrY
# Average recombination rate estimated in Locke et al (doi:10.1038/nature09687)
# which reports 0.95 +/- 0.72 cM/Mb
_chromosome_data = """\
chr1    227913704  0.95e-8
chr2A	109511694  0.95e-8
chr2B	129937803  0.95e-8
chr3	193656255  0.95e-8
chr4	189387572  0.95e-8
chr5	179185813  0.95e-8
chr6	169501136  0.95e-8
chr7	145408105  0.95e-8
chr8	144036388  0.95e-8
chr9	112206110  0.95e-8
chr10	132178492  0.95e-8
chr11	128122151  0.95e-8
chr12	132184051  0.95e-8
chr13	98475126   0.95e-8
chr14	88963417   0.95e-8
chr15	82547911   0.95e-8
chr16	68237989   0.95e-8
chr17	75914007   0.95e-8
chr18	75923960   0.95e-8
chr19	57575784   0.95e-8
chr20	60841859   0.95e-8
chr21	34683425   0.95e-8
chr22	35308119   0.95e-8
chrX	151242693  0.95e-8
"""

_chromosomes = []
for line in _chromosome_data.splitlines():
    name, length, mean_rr = line.split()[:3]
    _chromosomes.append(stdpopsim.Chromosome(
        id=name, length=int(length),
        mutation_rate=2.0e-8,
        recombination_rate=float(mean_rr)))

_genome = stdpopsim.Genome(chromosomes=_chromosomes)

_species = stdpopsim.Species(
    id="pongo",
    name="Pongo",
    genome=_genome,
    generation_time=20,
    population_size=1.79e4)

stdpopsim.register_species(_species)


###########################################################
#
# Genetic maps
#
###########################################################


# To do


###########################################################
#
# Demographic models
#
###########################################################


class LockeEtAlPongoIM(models.Model):
    id = "orangutan"
    name = "Two population orangutan model"
    description = """i
        The two orang-utan species, Sumatran (Pongo abelii) and Bornean (Pongo
        pygmaeus) inferred from the joint-site frequency spectrum with ten
        individuals from each population. This model is an isolation-with-
        migration model, with exponential growth or decay in each population
        after the split. The Sumatran population grows in size, while the
        Bornean population slightly declines.
    """

    citations = [stdpopsim.Citation(
        author="Locke et al.",
        year=2011,
        doi="http://doi.org/10.1038/nature09687")
    ]

    populations = [
        stdpopsim.Population(
            "Bornean", "Pongo pygmaeus (Bornean) population"),
        stdpopsim.Population(
            "Sumatran", "Pongo abelii (Sumatran) population")
    ]

    def __init__(self):
        super().__init__()

        # Parameters from paper:
        # ancestral size, before split
        Na = 17934

        # time of split
        T_split_years = 403149
        # get split time in units of generations
        generation_time = 20
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

        self.population_configurations = [
            msprime.PopulationConfiguration(initial_size=N_B, growth_rate=r_B,
                                            metadata=self.populations[0].asdict()),  # NOQA
            msprime.PopulationConfiguration(initial_size=N_S, growth_rate=r_S,
                                            metadata=self.populations[1].asdict())  # NOQA
        ]
        self.migration_matrix = [
            [      0, m_B_S],  # NOQA
            [m_S_B,       0],  # NOQA
        ]
        self.demographic_events = [
            # Merge populations and turn off migration, change to size Na
            msprime.MassMigration(
                time=T_split, source=1, destination=0, proportion=1.0),
            msprime.MigrationRateChange(time=T_split, rate=0),
            msprime.PopulationParametersChange(
                time=T_split, initial_size=Na, growth_rate=0, population_id=0),
        ]
