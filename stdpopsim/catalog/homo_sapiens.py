"""
Catalog definitions for
"""
import math
import logging

import msprime

import stdpopsim

logger = logging.getLogger(__name__)

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
    _chromosomes.append(stdpopsim.Chromosome(
        name=name, length=int(length),
        mutation_rate=1e-8,  # WRONG!,
        recombination_rate=float(mean_rr)))

_genome = stdpopsim.Genome(chromosomes=_chromosomes)

_species = stdpopsim.Species(
    id="homsap",
    name="Homo sapiens",
    genome=_genome,
    # TODO reference for these
    generation_time=25,
    population_size=10**4)

stdpopsim.register_species(_species)


###########################################################
#
# Genetic maps
#
###########################################################


_gm = stdpopsim.GeneticMap(
    species=_species,
    name="HapmapII_GRCh37",
    url=(
        "http://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/"
        "20110106_recombination_hotspots/"
        "HapmapII_GRCh37_RecombinationHotspots.tar.gz"),
    file_pattern="genetic_map_GRCh37_{name}.txt",
    description=(
        "The Phase II HapMap Genetic map (lifted over to GRCh37) used in "
        "1000 Genomes. Please see the `README <ftp://ftp-trace.ncbi.nih.gov/1000genomes"
        "/ftp/technical/working/20110106_recombination_hotspots"
        "/README_hapmapII_GRCh37_map>`_ for more details."),
    citations=[stdpopsim.Citation(
        doi="https://doi.org/10.1038/nature06258",
        year=2007,
        author="1000 Genomes Project consortium")]
    )
_species.add_genetic_map(_gm)

_gm = stdpopsim.GeneticMap(
    species=_species,
    name="Decode_2010_sex_averaged",
    url=(
        "http://sesame.uoregon.edu/~adkern/stdpopsim/decode/"
        "decode_2010_sex-averaged_map.tar.gz"),
    file_pattern="genetic_map_decode_2010_sex-averaged_{name}.txt",
    description=(
        "Decode fine scale genetic map from Kong, A et al. Fine scale "
        "recombination rate differences between sexes, populations and "
        "individuals. Nature (28 October 2010). "
        "Please see https://www.decode.com/addendum/ for more details."),
    citations=[stdpopsim.Citation(
        year=2010,
        author="Kong et al",
        doi="https://doi.org/10.1038/nature09525")]
    )
_species.add_genetic_map(_gm)


###########################################################
#
# Demographic models
#
###########################################################

# population definitions that are reused.
_yri_population = stdpopsim.Population(
    name="YRI",
    description="1000 Genomes YRI (Yorubans)")
_ceu_population = stdpopsim.Population(
    name="CEU",
    description=(
        "1000 Genomes CEU (Utah Residents (CEPH) with Northern and "
        "Western European Ancestry"))
_chb_population = stdpopsim.Population(
    name="CHB",
    description="1000 Genomes CHB (Han Chinese in Beijing, China)")


_tennessen_et_al = stdpopsim.Citation(
    author="Tennessen et al.",
    year=2012,
    doi="https://doi.org/10.1126/science.1219240")


# TODO: remove this superclass
class HomoSapiensModel(stdpopsim.Model):
    species = _species
    """
    TODO: documentation
    """
    def __init__(self):
        super().__init__()
        self.generation_time = _species.generation_time


# TODO we want to move away from defining these as classes which
# can be instantiated and rather to creating *instances* of the Model
# class which has this behaviour. However, it's not clear what form
# the refactored versions would have. Marking these classes with a
# __ to emphasise that they're not supposed to be used externally
# like this, but should be found via the species catalog.


class _GutenkunstThreePopOutOfAfrica(HomoSapiensModel):
    id = "ooa_3"
    name = "Three population out-of-Africa"
    description = """
        The three population Out-of-Africa model from `Gutenkunst et al. <https://
        doi.org/10.1371/journal.pgen.1000695>`_ It describes the ancestral human
        population in Africa, the out of Africa event, and the subsequent European-Asian
        population split. Model parameters are the maximum likelihood values of the
        various parameters given in Table 1 of Gutenkunst et al.
    """
    populations = [
        _yri_population,
        _ceu_population,
        _chb_population
    ]

    citations = [stdpopsim.Citation(
        author="Gutenkunst et al.",
        year=2009,
        doi="https://doi.org/10.1371/journal.pgen.1000695")
    ]

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

        # self.generation_time = default_generation_time
        T_AF = 220e3 / self.generation_time
        T_B = 140e3 / self.generation_time
        T_EU_AS = 21.2e3 / self.generation_time
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
            msprime.PopulationConfiguration(
                initial_size=N_AF, metadata=self.populations[0].asdict()),
            msprime.PopulationConfiguration(
                initial_size=N_EU, growth_rate=r_EU,
                metadata=self.populations[1].asdict()),
            msprime.PopulationConfiguration(
                initial_size=N_AS, growth_rate=r_AS,
                metadata=self.populations[2].asdict()),
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


_species.add_model(_GutenkunstThreePopOutOfAfrica())


class _TennessenTwoPopOutOfAfrica(HomoSapiensModel):
    id = "ooa_2"
    name = "Two population out-of-Africa"
    description = """
        The model is derived from the Tennesen et al. `analysis <https://doi.org/10.1126/
        science.1219240>`_  of the jSFS from European Americans and African Americans.
        It describes the ancestral human population in Africa, the out of Africa event,
        and two distinct periods of subsequent European population growth over the past
        23kya. Model parameters are taken from Fig. S5 in `Fu et al. (2013) <https://
        doi.org/10.1038 nature11690>`_.
    """

    populations = [
        stdpopsim.Population(name="AFR", description="African Americans"),
        stdpopsim.Population(name="EUR", description="European Americans")
    ]
    citations = [_tennessen_et_al]

    def __init__(self):
        super().__init__()

        T_AF = 148e3 / self.generation_time
        T_OOA = 51e3 / self.generation_time
        T_EU0 = 23e3 / self.generation_time
        T_EG = 5115 / self.generation_time

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
            msprime.PopulationConfiguration(
                initial_size=N_AF, growth_rate=r_AF,
                metadata=self.populations[0].asdict()),
            msprime.PopulationConfiguration(
                initial_size=N_EU, growth_rate=r_EU,
                metadata=self.populations[1].asdict())
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


_species.add_model(_TennessenTwoPopOutOfAfrica())


class _TennessenOnePopAfrica(HomoSapiensModel):
    id = "african"
    name = "African population"
    description = """
        The model is a simplification of the two population Tennesen et al.
        `model <https://doi.org/10.1126/science.1219240>`_ with the European-American
        population removed so that we are modeling the African population in isolation.
    """
    populations = [
        stdpopsim.Population(name="AFR", description="African"),
    ]
    citations = [_tennessen_et_al]

    def __init__(self):
        super().__init__()

        T_AF = 148e3 / self.generation_time
        T_EG = 5115 / self.generation_time

        # Growth rate
        r_AF = 0.0166

        # population sizes
        N_A = 7310
        N_AF1 = 14474

        # present Ne
        N_AF = N_AF1 / math.exp(-r_AF * T_EG)

        self.population_configurations = [
            msprime.PopulationConfiguration(
                initial_size=N_AF, growth_rate=r_AF,
                metadata=self.populations[0].asdict()),
        ]

        self.migration_matrix = [
            [0]
        ]

        self.demographic_events = [
            msprime.PopulationParametersChange(
                time=T_EG, growth_rate=0, initial_size=N_AF1, population_id=0),
            msprime.PopulationParametersChange(
                time=T_AF, initial_size=N_A, population_id=0)
        ]


_species.add_model(_TennessenOnePopAfrica())


class _BrowningAmerica(HomoSapiensModel):
    id = "america"
    name = "American admixture"
    description = """
        Demographic model for American admixture, taken from
        `Browning et al. <http://dx.doi.org/10.1371/journal.pgen.1007385>`_.
        This model extends the `Gravel et al. (2011) <https://doi.org/10.1073/
        pnas.1019276108>`_ model of African/European/Asian demographic history to
        simulate an admixed population with admixture occurring 12 generations ago. The
        admixed population had an initial size of 30,000 and grew at a rate of 5% per
        generation, with 1/6 of the population of African ancestry, 1/3 European, and 1
        2 Asian. This code was ported over from `Supplementary File 1
        <https://doi.org/10.1371/journal.pgen.1007385.s005>`_
    """

    populations = [
        stdpopsim.Population(name="AFR", description="Contemporary African population"),
        stdpopsim.Population(name="EUR", description="Contemporary European population"),
        stdpopsim.Population(name="ASIA", description="Contemporary Asian population"),
        stdpopsim.Population(
            name="ADMIX", description="Ancient admixed population",
            allow_samples=False),
    ]

    citations = [
        stdpopsim.Citation(
            author="Browning et al.",
            year=2011,
            doi="http://dx.doi.org/10.1371/journal.pgen.1007385")
    ]

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
                initial_size=Naf, growth_rate=0.0,
                metadata=self.populations[0].asdict()),
            msprime.PopulationConfiguration(
                initial_size=Neu*math.exp(reu*Teu), growth_rate=reu,
                metadata=self.populations[1].asdict()),
            msprime.PopulationConfiguration(
                initial_size=Nas*math.exp(ras*Teu), growth_rate=ras,
                metadata=self.populations[2].asdict()),
            msprime.PopulationConfiguration(
                initial_size=Nadmix*math.exp(radmix*Tadmix), growth_rate=radmix,
                metadata=self.populations[3].asdict())
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


_species.add_model(_BrowningAmerica())


class _RagsdaleArchaic(HomoSapiensModel):
    id = "ooa_archaic"
    name = "Three population out-of-Africa with archaic admixture"
    description = """
        The three population out-of-African model popularized by Gutenkunst et al. (2009)
        and augmented by archaic contributions to both Eurasian and African populations.
        Two archaic populations split early in human history, before the African
        expansion, and contribute to Eurasian populations (putative Neanderthal branch)
        and to the African branch (a deep diverging branch within Africa). Admixture
        is modeled as symmetric migration between the archaic and modern human branches,
        with contribution ending at a given time in the past.
    """
    populations = [
        _yri_population,
        _ceu_population,
        _chb_population,
        stdpopsim.Population(
            "Neanderthal", "Putative Neanderthals", allow_samples=False),
        stdpopsim.Population(
            "ArchaicAFR", "Putative Archaic Africans", allow_samples=False),
    ]
    citations = [
        stdpopsim.Citation(
            author="Ragsdale and Gravel",
            year=2019,
            doi="https://doi.org/10.1371/journal.pgen.1008204")
    ]

    def __init__(self):
        super().__init__()

        # First we set out the maximum likelihood values of the various parameters
        # given in Table 1 (under archaic admixture).
        N_0 = 3600
        N_YRI = 13900
        N_B = 880
        N_CEU0 = 2300
        N_CHB0 = 650

        # Times are provided in years, so we convert into generations.
        # In the published model, the authors used a generation time of 29 years to
        # convert from genetic to physical units
        self.generation_time = 29
        T_AF = 300e3 / self.generation_time
        T_B = 60.7e3 / self.generation_time
        T_EU_AS = 36.0e3 / self.generation_time
        T_arch_afr_split = 499e3 / self.generation_time
        T_arch_afr_mig = 125e3 / self.generation_time
        T_nean_split = 559e3 / self.generation_time
        T_arch_adm_end = 18.7e3 / self.generation_time

        # We need to work out the starting (diploid) population sizes based on
        # the growth rates provided for these two populations
        r_CEU = 0.00125
        r_CHB = 0.00372
        N_CEU = N_CEU0 / math.exp(-r_CEU * T_EU_AS)
        N_CHB = N_CHB0 / math.exp(-r_CHB * T_EU_AS)

        # Migration rates during the various epochs.
        m_AF_B = 52.2e-5
        m_YRI_CEU = 2.48e-5
        m_YRI_CHB = 0e-5
        m_CEU_CHB = 11.3e-5
        m_AF_arch_af = 1.98e-5
        m_OOA_nean = 0.825e-5

        # Population IDs correspond to their indexes in the population
        # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
        # initially.
        # We also have two archaic populations, putative Neanderthals and
        # archaicAfrican, which are population indices 3=Nean and 4=arch_afr.
        # Their sizes are equal to the ancestral reference population size N_0.
        self.population_configurations = [
            msprime.PopulationConfiguration(
                initial_size=N_YRI, metadata=self.populations[0].asdict()),
            msprime.PopulationConfiguration(
                initial_size=N_CEU, growth_rate=r_CEU,
                metadata=self.populations[1].asdict()),
            msprime.PopulationConfiguration(
                initial_size=N_CHB, growth_rate=r_CHB,
                metadata=self.populations[2].asdict()),
            msprime.PopulationConfiguration(
                initial_size=N_0,
                metadata=self.populations[3].asdict()),
            msprime.PopulationConfiguration(
                initial_size=N_0,
                metadata=self.populations[4].asdict())
        ]
        self.migration_matrix = [                   # noqa
            [      0, m_YRI_CEU, m_YRI_CHB, 0, 0],  # noqa
            [m_YRI_CEU,       0, m_CEU_CHB, 0, 0],  # noqa
            [m_YRI_CHB, m_CEU_CHB,       0, 0, 0],  # noqa
            [      0,         0,         0, 0, 0],  # noqa
            [      0,         0,         0, 0, 0]   # noqa
        ]                                           # noqa
        self.demographic_events = [
            # first event is migration turned on between modern and archaic humans
            msprime.MigrationRateChange(
                time=T_arch_adm_end, rate=m_AF_arch_af, matrix_index=(0, 4)),
            msprime.MigrationRateChange(
                time=T_arch_adm_end, rate=m_AF_arch_af, matrix_index=(4, 0)),
            msprime.MigrationRateChange(
                time=T_arch_adm_end, rate=m_OOA_nean, matrix_index=(1, 3)),
            msprime.MigrationRateChange(
                time=T_arch_adm_end, rate=m_OOA_nean, matrix_index=(3, 1)),
            msprime.MigrationRateChange(
                time=T_arch_adm_end, rate=m_OOA_nean, matrix_index=(2, 3)),
            msprime.MigrationRateChange(
                time=T_arch_adm_end, rate=m_OOA_nean, matrix_index=(3, 2)),

            # CEU and CHB merge into B with rate changes at T_EU_AS
            msprime.MassMigration(
                time=T_EU_AS, source=2, destination=1, proportion=1.0),
            msprime.MigrationRateChange(time=T_EU_AS, rate=0),
            msprime.MigrationRateChange(
                time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
            msprime.MigrationRateChange(
                time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
            msprime.MigrationRateChange(
                time=T_EU_AS, rate=m_AF_arch_af, matrix_index=(0, 4)),
            msprime.MigrationRateChange(
                time=T_EU_AS, rate=m_AF_arch_af, matrix_index=(4, 0)),
            msprime.MigrationRateChange(
                time=T_EU_AS, rate=m_OOA_nean, matrix_index=(1, 3)),
            msprime.MigrationRateChange(
                time=T_EU_AS, rate=m_OOA_nean, matrix_index=(3, 1)),
            msprime.PopulationParametersChange(
                time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),

            # Population B merges into YRI at T_B
            msprime.MassMigration(
                time=T_B, source=1, destination=0, proportion=1.0),
            msprime.MigrationRateChange(time=T_B, rate=0),
            msprime.MigrationRateChange(
                time=T_B, rate=m_AF_arch_af, matrix_index=(0, 4)),
            msprime.MigrationRateChange(
                time=T_B, rate=m_AF_arch_af, matrix_index=(4, 0)),

            # Beginning of migration between African and archaic African populations
            msprime.MigrationRateChange(time=T_arch_afr_mig, rate=0),

            # Size changes to N_0 at T_AF
            msprime.PopulationParametersChange(
                time=T_AF, initial_size=N_0, population_id=0),

            # Archaic African merges with moderns
            msprime.MassMigration(
                time=T_arch_afr_split, source=4, destination=0, proportion=1.0),

            # Neanderthal merges with moderns
            msprime.MassMigration(
                time=T_nean_split, source=3, destination=0, proportion=1.0)
        ]


_species.add_model(_RagsdaleArchaic())


class _SchiffelsZigzag(HomoSapiensModel):
    id = "zigzag"
    name = "Periodic growth and decline."
    description = """
        A validation model used by Schiffels and Durbin (2014) and Terhorst and
        Terhorst, Kamm, and Song (2017) with periods of exponential growth and
        decline in a single population.
        """
    populations = [
        stdpopsim.Population("generic", "Generic expanding and contracting population"),
    ]
    citations = [
        stdpopsim.Citation(
            author="Schiffels and Durbin",
            year=2014,
            doi="https://doi.org/10.1038/ng.3015")
    ]

    def __init__(self):
        super().__init__()

        self.generation_time = 29
        N0 = 14312

        g_1 = 0.023025
        t_1 = 33.333
        n_1 = N0

        g_2 = -0.005756
        t_2 = 133.33
        n_2 = N0/10

        g_3 = 0.0014391
        t_3 = 533.33
        n_3 = N0

        g_4 = -0.00035977
        t_4 = 2133.33
        n_4 = N0/10

        g_5 = 8.99448e-5
        t_5 = 8533.33
        n_5 = N0

        n_ancient = N0/10
        t_ancient = 34133.31

        self.population_configurations = [
            msprime.PopulationConfiguration(
                initial_size=N0)
        ]

        self.migration_matrix = [
            [0]
        ]

        self.demographic_events = [
                msprime.PopulationParametersChange(
                    initial_size=n_1, time=t_1, growth_rate=g_1),
                msprime.PopulationParametersChange(
                    initial_size=n_2, time=t_2, growth_rate=g_2),
                msprime.PopulationParametersChange(
                    initial_size=n_3, time=t_3, growth_rate=g_3),
                msprime.PopulationParametersChange(
                    initial_size=n_4, time=t_4, growth_rate=g_4),
                msprime.PopulationParametersChange(
                    initial_size=n_5, time=t_5, growth_rate=g_5),
                msprime.PopulationParametersChange(
                    time=t_ancient, initial_size=n_ancient, growth_rate=0)
        ]


_species.add_model(_SchiffelsZigzag())
