"""
Genome, genetic map and demographic model definitions for Arabidopsis thaliana.
"""
import msprime
import numpy as np

import stdpopsim

###########################################################
#
# Genome definition
#
###########################################################

# Data for length information based on:
# https://www.arabidopsis.org/portals/genAnnotation/
#   gene_structural_annotation/agicomplete.jsp
# Lengths from TAIR 10 although its unclear what reference the genetic map used
#   -- follow up on this with Salome 2012 authors

_chromosome_data = """\
chr1 30427671
chr2 19698289
chr3 23459830
chr4 18585056
chr5 26975502
"""
# mutation rate from Ossowski 2010 Science
# recombination value from Huber et al 2014 MBE
# rho=200/Mb, assume Ne=124,000, rho=2*Ne*r
_chromosomes = []
for line in _chromosome_data.splitlines():
    name, length = line.split()[:2]
    _chromosomes.append(stdpopsim.Chromosome(
        id=name, length=int(length),
        mutation_rate=7e-9,
        recombination_rate=8.1e-9))

_SwarbreckEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1093/nar/gkm965",
    year="2007",
    author="Swarbreck et al.",
    reasons={stdpopsim.CiteReason.ASSEMBLY}
)

_genome = stdpopsim.Genome(
        chromosomes=_chromosomes,
        assembly_citations=[
            _SwarbreckEtAl])

_species = stdpopsim.Species(
    id="AraTha",
    name="Arabidopsis thaliana",
    common_name="A. thaliana",
    genome=_genome,
    generation_time=1.0,
    generation_time_citations=[stdpopsim.Citation(
        doi="https://doi.org/10.1890/0012-9658(2002)083[1006:GTINSO]2.0.CO;2",
        year="2002",
        author="Donohue",
        reasons={stdpopsim.CiteReason.GEN_TIME})],
    population_size=10**4,
    population_size_citations=[stdpopsim.Citation(
        doi="https://doi.org/10.1016/j.cell.2016.05.063",
        year="2016",
        author="1001GenomesConsortium",
        reasons={stdpopsim.CiteReason.POP_SIZE})]
    )

stdpopsim.register_species(_species)

###########################################################
#
# Genetic maps
#
###########################################################

_gm = stdpopsim.GeneticMap(
    species=_species,
    id="SalomeAveraged_TAIR7",
    description="Crossover frequency map averaged over 17 populations",
    long_description="""
        This map is based on the study of crossover frequencies in over 7000
        plants in 17 F2 populations derived from crosses between 18 A. thaliana
        accessions. Salomé et al provide genetic maps for each of these
        populations. To get a single map for each chromosome, the Haldane map function
        distances were converted to recombination rates (cM/Mb) for each cross
        and then averaged across the 17 populations using loess.
        """,
    url=(
        "http://www.eeb.ucla.edu/Faculty/Lohmueller/data/"
        "uploads/salome2012_maps.tar.gz"),
    file_pattern="arab_{id}_map_loess.txt",
    citations=[stdpopsim.Citation(
        doi="https://doi.org/10.1038/hdy.2011.95",
        author="Salomé et al.",
        year=2012,
        reasons={stdpopsim.CiteReason.GEN_MAP})]
    )
_species.add_genetic_map(_gm)


###########################################################
#
# Demographic models
#
###########################################################


# FIXME this documentation needs to be filled out.
def _sma_1pop():
    # the size during the interval times[k] to times[k+1] = sizes[k]
    times = np.array([
        699, 2796, 6068, 9894, 14370, 19606, 25730, 32894, 41275,
        51077, 62544, 75958, 91648, 110001, 131471, 156584, 185960, 220324,
        260520, 307540, 362541, 426879, 502139, 590173, 693151, 813610,
        954517, 1119341, 1312147, 1537686, 1801500, 2110100])
    sizes = np.array([
        42252426, 42252426, 60323, 72174, 40591, 21158, 21442,
        39942, 78908, 111132, 110745, 96283, 87661, 83932, 83829, 91813,
        111644, 143456, 181571, 217331, 241400, 246984, 238593, 228222,
        217752, 198019, 165210, 121796, 121796, 73989, 73989, 73989])

    # MSMC is accurate from 40Kya-1.6Mya for A.thaliana (Durvasula et al 2017)
    # set the first 7 sizes
    # equal to the size at 8 (~40Kya)
    sizes[:8] = sizes[8]
    # set the last 2 entries equal
    # to the size at 30 (~1.6Mya)
    sizes[30:32] = sizes[30]

    demographic_events = []
    for sz, t in zip(sizes, times):
        demographic_events.append(
            msprime.PopulationParametersChange(
                time=t, initial_size=sz, population_id=0))

    populations = [
        stdpopsim.Population(
            id="SouthMiddleAtlas",
            description="Arabidopsis Thaliana South Middle Atlas population")
    ]

    return stdpopsim.DemographicModel(
        id="SouthMiddleAtlas_1D17",
        description="South Middle Atlas piecewise constant size",
        # TODO more detail here. We should at least explain that we're skipping
        # some of the estimates from MSMC because we don't think they're accurate.
        long_description="""
            Model estimated from two homozygous individuals from
            the South Middle Atlas using MSMC (TODO: more detail).
        """,
        populations=populations,
        citations=[stdpopsim.Citation(
            author="Durvasula et al.",
            year=2017,
            doi="https://doi.org/10.1073/pnas.1616736114",
            reasons={stdpopsim.CiteReason.DEM_MODEL})
        ],
        generation_time=1,
        demographic_events=demographic_events,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=sizes[0], metadata=populations[0].asdict())
            ]
        )


_species.add_demographic_model(_sma_1pop())


def _afr_2epoch():
    N_A = 746148
    N_0 = 100218
    t_1 = 568344
    populations = [
        stdpopsim.Population(
            id="Africa", description="Arabidopsis thaliana African population")
    ]
    return stdpopsim.DemographicModel(
        id="African2Epoch_1H18",
        description="African two epoch model",
        long_description="""
            Model estimated from site frequency spectrum of synonymous
            SNPs from African samples using Williamson et al. (2005)
            methodology.
        """,
        populations=populations,
        citations=[stdpopsim.Citation(
            author="Huber et al.",
            year=2018,
            doi="https://doi.org/10.1038/s41467-018-05281-7",
            reasons={stdpopsim.CiteReason.DEM_MODEL})
        ],
        generation_time=1,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_0, metadata=populations[0].asdict())
            ],
        demographic_events=[
            msprime.PopulationParametersChange(
                time=t_1, initial_size=N_A, population_id=0)
            ]
        )


_species.add_demographic_model(_afr_2epoch())


def _afr_3epoch():
    # values from the supplementary Table 1.
    # the size changed as N_A -> N_2 -> N_3.
    # t_2 is time of 2nd epoch and t_3 of the third epoch
    N_A = 161744
    N_2 = 24076
    N_3 = 203077
    t_2 = 7420
    t_3 = 14534
    populations = [
        stdpopsim.Population(
            id="Africa", description="Arabidopsis thaliana African population")
    ]
    return stdpopsim.DemographicModel(
        id="African3Epoch_1H18",
        description="African three epoch model",
        long_description="""
            Model estimated from site frequency spectrum of synonymous
            SNPs from African samples using Williamson et al. (2005)
            methodology.
        """,
        populations=populations,
        citations=[stdpopsim.Citation(
            author="Huber et al.",
            year=2018,
            doi="https://doi.org/10.1038/s41467-018-05281-7",
            reasons={stdpopsim.CiteReason.DEM_MODEL})
        ],
        generation_time=1,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_3, metadata=populations[0].asdict())],
        demographic_events=[
            msprime.PopulationParametersChange(
                time=t_3, initial_size=N_2, population_id=0),
            msprime.PopulationParametersChange(
                time=t_2 + t_3, initial_size=N_A, population_id=0)]
        )


_species.add_demographic_model(_afr_3epoch())
