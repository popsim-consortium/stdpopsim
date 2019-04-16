"""
Genome, genetic map and demographic model definitions for humans.
"""
import math

import msprime

import stdpopsim.models as models
import stdpopsim.genomes as genomes
import stdpopsim.genetic_maps as genetic_maps


###########################################################
#
# Genetic maps
#
###########################################################


class HapmapII_GRCh37(genetic_maps.GeneticMap):
    """
    Usage: `hapmap = homo_sapiens.HapmapII_GRCh37()` (note the
    parentheses).

    The Phase II HapMap Genetic map (lifted over to GRCh37) used in
    1000 Genomes. Please see the `README
    <ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20110106_recombination_hotspots/README_hapmapII_GRCh37_map>`_
    for more details.
    """
    url = (
        "http://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/"
        "20110106_recombination_hotspots/"
        "HapmapII_GRCh37_RecombinationHotspots.tar.gz")
    file_pattern = "genetic_map_GRCh37_{name}.txt"


genetic_maps.register_genetic_map(HapmapII_GRCh37())


###########################################################
#
# Genome definition
#
###########################################################

# List of chromosomes.

# FIXME: add mean mutation rate data to this table.
# Name  Length  mean_recombination_rate mean_mutation_rate

# length information can be found here
# <http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz>

# mean_recombination_rate was computed across all windows of the GRCh37 genetic map
# <ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20110106_recombination_hotspots>
_chromosome_data = """\
chr1 	 249250621 	 1.1485597641285933e-08
chr2 	 243199373 	 1.1054289277533446e-08
chr3 	 198022430 	 1.1279585624662551e-08
chr4 	 191154276 	 1.1231162636001008e-08
chr5 	 180915260 	 1.1280936570022824e-08
chr6 	 171115067 	 1.1222852661225285e-08
chr7 	 159138663 	 1.1764614397655721e-08
chr8 	 146364022 	 1.1478465778920576e-08
chr9 	 141213431 	 1.1780701596308656e-08
chr10 	 135534747 	 1.3365134257075317e-08
chr11 	 135006516 	 1.1719334320833283e-08
chr12 	 133851895 	 1.305017186986983e-08
chr13 	 115169878 	 1.0914860554958317e-08
chr14 	 107349540 	 1.119730771394731e-08
chr15 	 102531392 	 1.3835785893339787e-08
chr16 	 90354753 	 1.4834607113882717e-08
chr17 	 81195210 	 1.582489036239487e-08
chr18 	 78077248 	 1.5075956950023575e-08
chr19 	 59128983 	 1.8220141872466202e-08
chr20 	 63025520 	 1.7178269031631664e-08
chr21 	 48129895 	 1.3045214034879191e-08
chr22 	 51304566 	 1.4445022767788226e-08
chrX 	 155270560 	 1.164662223273842e-08
chrY 	 59373566 	 0.0
"""

_chromosomes = []
for line in _chromosome_data.splitlines():
    name, length, mean_rr = line.split()[:3]
    _chromosomes.append(genomes.Chromosome(
        name=name, length=int(length),
        default_mutation_rate=1e-8,  # WRONG!,
        default_recombination_rate=float(mean_rr)))


#: :class:`stdpopsim.Genome` definition for humans.
genome = genomes.Genome(
    species="homo_sapiens",
    chromosomes=_chromosomes,
    default_genetic_map=HapmapII_GRCh37.name)


###########################################################
#
# Demographic models
#
###########################################################


class GutenkunstThreePopOutOfAfrica(models.Model):
    """
    The three population Out-of-Africa model from Gutenkunst et al.

    .. todo:: document this model, including the original publications
        and clear information about what the different population indexes
        mean.

    """

    def __init__(self):
        super().__init__()
        # First we set out the maximum likelihood values of the various parameters
        # given in Table 1.
        N_A = 7300
        N_B = 2100
        N_AF = 12300
        N_EU0 = 1000
        N_AS0 = 510
        # Times are provided in years, so we convert into generations.
        generation_time = 25
        T_AF = 220e3 / generation_time
        T_B = 140e3 / generation_time
        T_EU_AS = 21.2e3 / generation_time
        # We need to work out the starting (diploid) population sizes based on
        # the growth rates provided for these two populations
        r_EU = 0.004
        r_AS = 0.0055
        N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
        N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
        # Migration rates during the various epochs.
        m_AF_B = 25e-5
        m_AF_EU = 3e-5
        m_AF_AS = 1.9e-5
        m_EU_AS = 9.6e-5
        # Population IDs correspond to their indexes in the population
        # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
        # initially.
        self.population_configurations = [
            msprime.PopulationConfiguration(initial_size=N_AF),
            msprime.PopulationConfiguration(initial_size=N_EU, growth_rate=r_EU),
            msprime.PopulationConfiguration(initial_size=N_AS, growth_rate=r_AS)
        ]
        self.migration_matrix = [
            [      0, m_AF_EU, m_AF_AS],  # noqa
            [m_AF_EU,       0, m_EU_AS],  # noqa
            [m_AF_AS, m_EU_AS,       0],  # noqa
        ]
        self.demographic_events = [
            # CEU and CHB merge into B with rate changes at T_EU_AS
            msprime.MassMigration(
                time=T_EU_AS, source=2, destination=1, proportion=1.0),
            msprime.MigrationRateChange(time=T_EU_AS, rate=0),
            msprime.MigrationRateChange(
                time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
            msprime.MigrationRateChange(
                time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
            msprime.PopulationParametersChange(
                time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),
            # Population B merges into YRI at T_B
            msprime.MassMigration(
                time=T_B, source=1, destination=0, proportion=1.0),
            # Size changes to N_A at T_AF
            msprime.PopulationParametersChange(
                time=T_AF, initial_size=N_A, population_id=0)
        ]


class TennessenTwoPopOutOfAfrica(models.Model):
    """
    The model is derived from the Tennesen et al.
    `analysis <https://doi.org/10.1126/science.1219240>`_  of the jSFS from
    European Americans and African Americans.

    Model parameters are taken from Fig. S5 in
    `Fu et al. (2013) <https://doi.org/10.1038/nature11690>`_.

    .. todo:: document this model, including the original publications
        and clear information about what the different population indexes
        mean.

    """

    def __init__(self):
        super().__init__()

        generation_time = 25
        T_AF = 148e3 / generation_time
        T_OOA = 51e3 / generation_time
        T_EU0 = 23e3 / generation_time
        T_EG = 5115 / generation_time

        # Growth rates
        r_EU0 = 0.00307
        r_EU = 0.0195
        r_AF = 0.0166

        # population sizes
        N_A = 7310
        N_AF1 = 14474
        N_B = 1861
        N_EU0 = 1032
        N_EU1 = N_EU0 / math.exp(-r_EU0 * (T_EU0-T_EG))

        # migration rates
        m_AF_B = 15e-5
        m_AF_EU = 2.5e-5

        # present Ne
        N_EU = N_EU1 / math.exp(-r_EU * T_EG)
        N_AF = N_AF1 / math.exp(-r_AF * T_EG)

        self.population_configurations = [
            msprime.PopulationConfiguration(initial_size=N_AF, growth_rate=r_AF),
            msprime.PopulationConfiguration(initial_size=N_EU, growth_rate=r_EU)
        ]

        self.migration_matrix = [
            [0, m_AF_EU],
            [m_AF_EU, 0],
        ]

        self.demographic_events = [
            msprime.MigrationRateChange(
                time=T_EG, rate=m_AF_EU, matrix_index=(0, 1)),
            msprime.MigrationRateChange(
                time=T_EG, rate=m_AF_EU, matrix_index=(1, 0)),
            msprime.PopulationParametersChange(
                time=T_EG, growth_rate=r_EU0, initial_size=N_EU1, population_id=1),
            msprime.PopulationParametersChange(
                time=T_EG, growth_rate=0, initial_size=N_AF1, population_id=0),
            msprime.MigrationRateChange(
                time=T_EU0, rate=m_AF_B, matrix_index=(0, 1)),
            msprime.MigrationRateChange(
                time=T_EU0, rate=m_AF_B, matrix_index=(1, 0)),
            msprime.PopulationParametersChange(
                time=T_EU0, initial_size=N_B, growth_rate=0, population_id=1),
            msprime.MassMigration(
                time=T_OOA, source=1, destination=0, proportion=1.0),
            msprime.PopulationParametersChange(
                time=T_AF, initial_size=N_A, population_id=0)
        ]


class BrowningAmerica(models.Model):
    """
    Demographic model for American admixture, taken from
    `Browning et al. <http://dx.doi.org/10.1371/journal.pgen.1007385>`_.
    This model extends the Gravel et al. (2011) model to simulate an admixed
    population with admixture occurring 12 generations ago. The admixed population
    had an initial size of 30,000 and grew at a rate of 5% per generation,
    with 1/6 of the population of African ancestry, 1/3 European, and 1/2 Asian.
    This code was ported over from
    `Supplementary File 1 <https://doi.org/10.1371/journal.pgen.1007385.s005>`_
    """
    def __init__(self):
        super().__init__()
        N0 = 7310  # initial population size
        Thum = 5920  # time (gens) of advent of modern humans
        Naf = 14474  # size of african population
        Tooa = 2040  # number of generations back to Out of Africa
        Nb = 1861  # size of out of Africa population
        mafb = 1.5e-4  # migration rate Africa and Out-of-Africa
        Teu = 920  # number generations back to Asia-Europe split
        Neu = 1032  # bottleneck population sizes
        Nas = 554
        mafeu = 2.5e-5  # mig. rates
        mafas = 7.8e-6
        meuas = 3.11e-5
        reu = 0.0038  # growth rate per generation in Europe
        ras = 0.0048  # growth rate per generation in Asia
        Tadmix = 12  # time of admixture
        Nadmix = 30000  # initial size of admixed population
        radmix = .05  # growth rate of admixed population
        # pop0 is Africa, pop1 is Europe, pop2 is Asia,  pop3 is admixed
        self.population_configurations = [
            msprime.PopulationConfiguration(
                initial_size=Naf, growth_rate=0.0),
            msprime.PopulationConfiguration(
                initial_size=Neu*math.exp(reu*Teu), growth_rate=reu),
            msprime.PopulationConfiguration(
                initial_size=Nas*math.exp(ras*Teu), growth_rate=ras),
            msprime.PopulationConfiguration(
                initial_size=Nadmix*math.exp(radmix*Tadmix), growth_rate=radmix)
        ]

        self.migration_matrix = [
            [0, mafeu, mafas, 0],
            [mafeu, 0, meuas, 0],
            [mafas, meuas, 0, 0],
            [0, 0, 0, 0]
        ]
        # Admixture event, 1/6 Africa, 2/6 Europe, 3/6 Asia
        admixture_event = [
            msprime.MassMigration(
                time=Tadmix, source=3, destination=0, proportion=1.0/6.0),
            msprime.MassMigration(
                time=Tadmix+0.0001, source=3, destination=1, proportion=2.0/5.0),
            msprime.MassMigration(
                time=Tadmix+0.0002, source=3, destination=2, proportion=1.0)
        ]
        # Asia and Europe split
        eu_event = [
            msprime.MigrationRateChange(
                time=Teu, rate=0.0),
            msprime.MassMigration(
                time=Teu+0.0001, source=2, destination=1, proportion=1.0),
            msprime.PopulationParametersChange(
                time=Teu+0.0002, initial_size=Nb, growth_rate=0.0, population_id=1),
            msprime.MigrationRateChange(
                time=Teu+0.0003, rate=mafb, matrix_index=(0, 1)),
            msprime.MigrationRateChange(
                time=Teu+0.0003, rate=mafb, matrix_index=(1, 0))
        ]
        # Out of Africa event
        ooa_event = [
            msprime.MigrationRateChange(
                time=Tooa, rate=0.0),
            msprime.MassMigration(
                time=Tooa+0.0001, source=1, destination=0, proportion=1.0)
        ]
        # initial population size
        init_event = [
            msprime.PopulationParametersChange(
                time=Thum,
                initial_size=N0,
                population_id=0)
        ]
        self.demographic_events = admixture_event + eu_event + ooa_event + init_event
