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


def _HolsteinFriesan_1M13():
    id = "HolsteinFriesian_1M13"
    description = "Piecewise size changes in Holstein-Friesian cattle."
    long_description = """
    The piecewise-constant population size model of Holstein-Friesian cattle
    from MacLeod et al. 2013. Population sizes were estimated from inferred
    runs of homozygosity, with parameter values taken from Figure 4A by visual
    inspection of the plots.
    """
    populations = [
        stdpopsim.Population(id="Holstein_Friesian", description="Holstein Friesian"),
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
            msprime.PopulationConfiguration(
                initial_size=90, growth_rate=0.0166, metadata=populations[0].asdict()
            )
        ],
        # Here 'time' should be in generation notation ie. how many
        # generations ago were that Ne (effective population size)
        # and growth rate.
        # Growth rate is "per generation exponential growth rate":
        # -alpha= [ln(initial_pop_size/next_stage_pop_size)/generation_span_in_years]
        # For example: ln(90/120)/3= -0.095894024
        demographic_events=[
            msprime.PopulationParametersChange(
                time=1,
                initial_size=90,
                growth_rate=-0.095894024,
                population_id=0,
            ),  # Ne 90 to 120
            msprime.PopulationParametersChange(
                time=4, growth_rate=-0.24465639, population_id=0
            ),  # Ne 120 to 250
            msprime.PopulationParametersChange(
                time=7, growth_rate=-0.0560787, population_id=0
            ),  # Ne 250 to 350
            msprime.PopulationParametersChange(
                time=13, growth_rate=-0.1749704, population_id=0
            ),  # Ne 350 to 1000
            msprime.PopulationParametersChange(
                time=19, growth_rate=-0.0675775, population_id=0
            ),  # Ne 1000 to 1500
            msprime.PopulationParametersChange(
                time=25, growth_rate=-0.0022129, population_id=0
            ),  # Ne 1500 to 2000
            msprime.PopulationParametersChange(
                time=155, growth_rate=-0.0007438, population_id=0
            ),  # Ne 2000 to 2500
            msprime.PopulationParametersChange(
                time=455, growth_rate=-0.0016824, population_id=0
            ),  # Ne 2500 to 3500
            msprime.PopulationParametersChange(
                time=655, growth_rate=-0.0006301, population_id=0
            ),  # Ne 3500 to 7000
            msprime.PopulationParametersChange(
                time=1755, growth_rate=-0.0005945, population_id=0
            ),  # Ne 7000 to 10000
            msprime.PopulationParametersChange(
                time=2355, growth_rate=-0.0005306, population_id=0
            ),  # Ne 10000 to 17000
            msprime.PopulationParametersChange(
                time=3355, growth_rate=-0.0000434, population_id=0
            ),  # Ne 17000 to 62000
            msprime.PopulationParametersChange(
                time=33155, growth_rate=-0.0000, population_id=0
            ),  # Ne 62000 (model has "coalesced")
            msprime.PopulationParametersChange(
                time=933155, growth_rate=-0.0, population_id=0
            ),
        ],
    )


_species.add_demographic_model(_HolsteinFriesan_1M13())
