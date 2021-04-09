import stdpopsim
from . import genome_data

# MyoLuc2.0
_LindbladTohEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1038/nature10530",
    year=2011,
    author="Lindblad-Toh, K., Garber, M., Zuk, O. et al.",
    reasons={stdpopsim.CiteReason.ASSEMBLY},
)

# Dumont and Payseur - recombination rate
_DumontAndPayseur = stdpopsim.Citation(
    doi="https://doi.org/10.1111/j.1558-5646.2007.00278.x",
    year=2007,
    author="Dumont, B. L. and Payseur, B. A.",
    reasons={stdpopsim.CiteReason.REC_RATE},
)

# Humphrey and Cope - generation time
_HumphreyAndCope = stdpopsim.Citation(
    doi="https://doi.org/10.5962/bhl.title.39539",
    year=1976,
    author="Humphrey, S. R. and Cope, J. B.",
    reasons={stdpopsim.CiteReason.GEN_TIME},
)

# Frick et al. - generation time
_FrickEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1111/j.1365-2656.2009.01615.x",
    year=2009,
    author="Frick, S.F., Reynolds, D. S., and Kunz, T. H.",
    reasons={stdpopsim.CiteReason.GEN_TIME},
)

# Ray et al. - mutation rate
_RayEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1101/gr.071886.107",
    year=2008,
    author="Ray, D.A., Feschotte, C., Pagan, H. J. T., et al.",
    reasons={stdpopsim.CiteReason.MUT_RATE},
)

# Davy et al. - effective population size
# We are settling on the Davy et al (2017) population size of 1682.9 (794.0–261 726.6)
# based on their sequences from Eastern Canada,
# as it is the only reasonable Ne estimate close to Massachusetts,
# the origin of the MyoLuc2.0 genome individual
_DavyEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1139/facets-2017-0022",
    year=2017,
    author="Davy, C. M., Donaldson, M. E., Rico, Y. et al.",
    reasons={stdpopsim.CiteReason.POP_SIZE},
)

# mutation rate:
# Ray et al. 2008 estimate 2.366 × 10^−9 per site per year
# Humphrey and Cope 1976, Frick et al. 2010 state generation time of 2 years
_overall_mutation_rate = 0.000000004732

# no recombination rate estimate available, use "standard" mammalian rate 10^-8
# similar to ~1cM/Mb for horses and cows
# (bats closest relatives) from Dumont and Payseur 2008
_overall_recombination_rate = 0.00000001

_chromosomes = []
for name, data in genome_data.data["chromosomes"].items():
    _chromosomes.append(
        stdpopsim.Chromosome(
            id=name,
            length=data["length"],
            synonyms=data["synonyms"],
            mutation_rate=_overall_mutation_rate,
            recombination_rate=_overall_recombination_rate,
        )
    )

_genome = stdpopsim.Genome(
    chromosomes=_chromosomes,
    citations=[_LindbladTohEtAl, _DumontAndPayseur, _RayEtAl],
)

_species = stdpopsim.Species(
    id="MyoLuc",
    ensembl_id="myotis_lucifugus",
    name="Myotis lucifugus",
    common_name="Microbat",
    genome=_genome,
    generation_time=2.0,
    population_size=1693,  # Davy et al. for Eastern Canada,
    # (794.0–261,726.6)... that's a big spread
    citations=[_DavyEtAl, _FrickEtAl, _HumphreyAndCope],
)

stdpopsim.register_species(_species)
