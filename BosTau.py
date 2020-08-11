"""
Catalog definitions for Bos Taurus
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


# Create a string of chromosome lengths for easy parsing
# Length data was extracted from Ensembl genome browser Eg. for Chr1 ....https://www.ensembl.org/Bos_taurus/Location/Chromosome?r=1:421347-158534110
_chromosome_data = """\
chr1    158534110
chr2    136231102
chr3    121005158
chr4    120000601
chr5    120089316
chr6    117806340
chr7    110682743
chr8    113319770
chr9    105454467
chr10   103308737
chr11   106982474
chr12   87216183
chr13   83472345
chr14   82403003
chr15   85007780
chr16   81013979
chr17   73167244
chr18   65820629
chr19   63449741
chr20   71974595
chr21   69862954
chr22   60773035
chr23   52498615
chr24   62317253
chr25   42350435
chr26   51992305
chr27   45612108
chr28   45940150
chr29   51098607
chrX    139009144
"""
# Parse list of chromosomes into a list of Chromosome objects which contain the
# chromosome name, length, mutation rate, and recombination rate
_chromosomes = []

for line in _chromosome_data.splitlines():
    name, length = line.split()[:2]
    _chromosomes.append(stdpopsim.Chromosome(
        id=name, length=int(length),
        mutation_rate=1.21e-8,
        recombination_rate=9.26e-9))

# A citation for the chromosome parameters. Additional citations may be needed if
# the mutation or recombination rates come from other sources. In that case create
# additional citations with the appropriate reasons specified (see API documentation
# for stdpopsim.citations)



_RosenEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1093/gigascience/giaa021",
    year="2020",
    author="Rosen et al.",
    reasons={stdpopsim.CiteReason.ASSEMBLY})

_CoppietersEtAl = stdpopsim.Citation(
        # Rate of de novo mutation in dairy cattle and potential impact of reproductive technologies. Proceedings of the World Congress on Genetics Applied to Livestock Production, 11. 983
        author="Coppieters et al.",
        year=2018,
        doi="")

_MaEtAl = stdpopsim.Citation(
        # Cattle Sex-Specific Recombination and Genetic Control from a Large Pedigree Analysis
        author="Ma et al.",
        year=2015,
        doi="https://doi.org/10.1371/journal.pgen.1005387")


# Create a genome object

_genome = stdpopsim.Genome(
    chromosomes=_chromosomes,
    # Mutation rate is a sex-averaged estimate per base pair per generation
     mutation_rate_citations=[
        _CoppietersEtAl.because(stdpopsim.CiteReason.MUT_RATE),
        ],
    # Recombination rate has been derived from dairy cattle crossovers per meiosis, by taking the average between females and males and then dividing by the 
    # whole genome length (equal to the sum of chromosome lengths used above)
    recombination_rate_citations=[
        _MaEtAl.because(stdpopsim.CiteReason.REC_RATE)
        ],
    assembly_citations=[
        _RosenEtAl.because(stdpopsim.CiteReason.ASSEMBLY)
        ],
  )

# Create a Species Object
_gen_time_citation = stdpopsim.Citation(
    doi="https://doi.org/10.1093/molbev/mst125",
    year="2013",
    author="MacLeod et al.",
    reasons={stdpopsim.CiteReason.GEN_TIME})

_pop_size_citation = stdpopsim.Citation(
    doi="https://doi.org/10.1093/molbev/mst125",
    year="2013",
    author="MacLeod et al.",
    reasons={stdpopsim.CiteReason.POP_SIZE})

_species = stdpopsim.Species(
    id="BosTau",
    name="Bos Taurus",
    common_name="Cattle",
    genome=_genome,
    generation_time=5,
    generation_time_citations=[_gen_time_citation],
    population_size=90,
    population_size_citations=[_pop_size_citation]
)

stdpopsim.register_species(_species)


###########################################################
#
# Demographic models
#
###########################################################

def _inferred_1_M_13():
    id = "IonaInferredDemography"
    description = "Iona MacLeod's Inferred Demographic Model for Bos Taurus"
    long_description = """
    The Runs of Homozygosity-based infer of Demography from MacLeod et al. 2013. Only the Holstein breed has been taken into account for this inference.
    """
    populations = [
        stdpopsim.Population(id="FILL ME", description="FILL ME"),
    ]
    citations = [
        stdpopsim.Citation(
            author="MacLeod et al.",
            year="2013",
            doi="https://doi.org/10.1093/molbev/mst125",
            reasons={stdpopsim.CiteReason.DEM_MODEL})
    ]

    generation_time = 5

    # parameter value definitions based on published values

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=90, growth_rate=0.0166,
                metadata=populations[0].asdict())
        ],
        
        
        demographic_events=[
            msprime.PopulationParametersChange(
                # Here 'time' should be in generation notation ie. how many generations ago were that Ne (effective population size) and growth rate.
                # Growth rate is "per generation exponential growth rate": -alpha= [ln(initial_pop_size/next_stage_pop_size)/generation_span_in_years]
                # For example: ln(90/120)/3= -0.095894024
                time=1, initial_size=90, growth_rate=-0.095894024, population_id=0), #Ne 90 to 120
            msprime.PopulationParametersChange(
                time=4, growth_rate=-0.24465639, population_id=0),  #Ne 120 to 250
            msprime.PopulationParametersChange(
                time=7, growth_rate=-0.0560787, population_id=0),   #Ne 250 to 350 
            msprime.PopulationParametersChange(
                time=13, growth_rate=-0.1749704, population_id=0),  #Ne 350 to 1000
            msprime.PopulationParametersChange(
                time=19, growth_rate=-0.0675775, population_id=0),  #Ne 1000 to 1500
            msprime.PopulationParametersChange(
                time=25, growth_rate=-0.0022129, population_id=0),  #Ne 1500 to 2000
            msprime.PopulationParametersChange(
                time=155, growth_rate=-0.0007438, population_id=0), #Ne 2000 to 2500
            msprime.PopulationParametersChange(
                time=455, growth_rate=-0.0016824, population_id=0), #Ne 2500 to 3500
            msprime.PopulationParametersChange(
                time=655, growth_rate=-0.0006301, population_id=0), #Ne 3500 to 7000
            msprime.PopulationParametersChange(
                time=1755, growth_rate=-0.0005945, population_id=0),    #Ne 7000 to 10000
            msprime.PopulationParametersChange(
                time=2355, growth_rate=-0.0005306, population_id=0),    #Ne 10000 to 17000
            msprime.PopulationParametersChange(
                time=3355, growth_rate=-0.0000434, population_id=0),    #Ne 17000 to 62000
            msprime.PopulationParametersChange(
                time=33155, growth_rate=-0.0000, population_id=0),      #Ne 62000 (model has "coalesced")
            msprime.PopulationParametersChange(
                time=933155, growth_rate=-0.0, population_id=0),

        ],
    )


_species.add_demographic_model(_inferred_1_M_13())
