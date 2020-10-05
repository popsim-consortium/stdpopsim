"""
Catalog definitions for Bos Taurus
"""
import collections

import msprime

import stdpopsim
from . import genome_data

###########################################################
#
# Genome definition
#
###########################################################

# De novo assembly of the cattle reference genome with single-molecule sequencing.
_RosenEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1093/gigascience/giaa021",
    year="2020",
    author="Rosen et al.",
)

# Frequency of mosaicism points towards mutation-prone early cleavage
# cell divisions in cattle.
_HarlandEtAl = stdpopsim.Citation(
    author="Harland et al.",
    year="2017",
    # BioRxiv preprint
    doi="https://doi.org/10.1101/079863",
)

# Cattle Sex-Specific Recombination and Genetic Control from a
# Large Pedigree Analysis.
_MaEtAl = stdpopsim.Citation(
    author="Ma et al.",
    year="2015",
    doi="https://doi.org/10.1371/journal.pgen.1005387",
)

# Inferring Demography from Runs of Homozygosity in Whole-Genome Sequence,
# with Correction for Sequence Errors.
_MacLeodEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1093/molbev/mst125",
    year="2013",
    author="MacLeod et al.",
)

# Recombination rate has been derived from dairy cattle crossovers
# per meiosis, by taking the average between females and males and then
# dividing by the whole genome length (equal to the sum of chromosome
# lengths used above).
# From Ma et al. (2015), 25.5 crossovers per meiosis in males and
# 23.2 crossovers per meiosis in females, gives an average of 24.35
# crossovers per meiosis. The sum of chromosome lengths is 2628394923 bp.
# 24.35 / 2628394923 = 9.26e-9 per bp per generation.
_genome_wide_recombination_rate = 9.26e-9

_recombination_rate_data = collections.defaultdict(
    lambda: _genome_wide_recombination_rate
)
# Set some exceptions for non-recombining chrs.
_recombination_rate_data["MT"] = 0

_chromosomes = []
for name, data in genome_data.data["chromosomes"].items():
    _chromosomes.append(
        stdpopsim.Chromosome(
            id=name,
            length=data["length"],
            synonyms=data["synonyms"],
            # Harland et al. (2017), sex-averaged estimate per bp per generation.
            mutation_rate=1.2e-8,
            recombination_rate=_recombination_rate_data[name],
        )
    )

_genome = stdpopsim.Genome(
    chromosomes=_chromosomes,
    mutation_rate_citations=[
        _HarlandEtAl.because(stdpopsim.CiteReason.MUT_RATE),
    ],
    recombination_rate_citations=[_MaEtAl.because(stdpopsim.CiteReason.REC_RATE)],
    assembly_citations=[_RosenEtAl.because(stdpopsim.CiteReason.ASSEMBLY)],
)

_species = stdpopsim.Species(
    id="BosTau",
    name="Bos Taurus",
    common_name="Cattle",
    genome=_genome,
    generation_time=5,
    generation_time_citations=[_MacLeodEtAl.because(stdpopsim.CiteReason.GEN_TIME)],
    population_size=62000,
    population_size_citations=[_MacLeodEtAl.because(stdpopsim.CiteReason.POP_SIZE)],
)

stdpopsim.register_species(_species)


###########################################################
#
# Demographic models
#
###########################################################


def _HolsteinFriesian_1M13():
    id = "HolsteinFriesian_1M13"
    description = "Piecewise size changes in Holstein-Friesian cattle."
    long_description = """
    The piecewise-constant population size model of Holstein-Friesian cattle
    from MacLeod et al. 2013. Population sizes were estimated from inferred
    runs of homozygosity, with parameter values taken from Figure 4A and Table S1.
    """
    populations = [
        stdpopsim.Population(id="Holstein-Friesian", description="Holstein-Friesian"),
    ]
    citations = [_MacLeodEtAl.because(stdpopsim.CiteReason.DEM_MODEL)]

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=_species.generation_time,
        population_configurations=[
            #      3 generations at    90,     1-     3
            msprime.PopulationConfiguration(initial_size=90,
                                            metadata=populations[0].asdict())
        ],
        # Here 'time' should be in generation notation ie. how many
        # generations ago were that Ne (effective population size)
        demographic_events=[
            #      3 generations at   120,     4-     6
            msprime.PopulationParametersChange(time=4, initial_size=120,
                population_id=0),
            #      6 generations at   250,     7-    12
            msprime.PopulationParametersChange(time=7, initial_size=250,
                population_id=0),
            #      6 generations at   350,    13-    18
            msprime.PopulationParametersChange(time=13, initial_size=350,
                population_id=0),
            #      6 generations at  1000,    19-    24
            msprime.PopulationParametersChange(time=19, initial_size=1000,
                population_id=0),
            #    130 generations at  1500,    25-   154
            msprime.PopulationParametersChange(time=25, initial_size=1500,
                population_id=0),
            #    200 generations at  2000,   155-   454
            msprime.PopulationParametersChange(time=155, initial_size=2000,
                population_id=0),
            #    200 generations at  2500,   455-   654
            msprime.PopulationParametersChange(time=455, initial_size=2500,
                population_id=0),
            #   1100 generations at  3500,   655-  1754
            msprime.PopulationParametersChange(time=655, initial_size=3500,
                population_id=0),
            #    600 generations at  7000,  1755-  2354
            msprime.PopulationParametersChange(time=1755, initial_size=7000,
                population_id=0),
            #   1000 generations at 10000,  2355-  3354
            msprime.PopulationParametersChange(time=2355, initial_size=10000,
                population_id=0),
            #  29800 generations at 17000,  3355- 33154
            msprime.PopulationParametersChange(time=3355, initial_size=17000,
                population_id=0),
            # 900000 generations at 62000, 33155-933154
            msprime.PopulationParametersChange(time=33155, initial_size=62000,
                population_id=0),
        ],
    )

_species.add_demographic_model(_HolsteinFriesian_1M13())


