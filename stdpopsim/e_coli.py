"""
Genome and demographic model definitions for Escherichia coli.
"""
import msprime

import stdpopsim.models as models
import stdpopsim.genomes as genomes


###########################################################
#
# Genome definition
#
###########################################################

_chromosomes = []
_chromosomes.append(genomes.Chromosome(
        name="chr",
        length=4641652,
        default_mutation_rate=1e-5+2e-4,
        default_recombination_rate=0.0))
# mean_conversion_rate=8.9e-11 # not implemented yet!
# mean_conversion_length=542 # not implemented yet!

#: :class:`stdpopsim.Genome` definition for E. Coli.
# Chromosome length data is based on strain K-12.
genome = genomes.Genome(
    species="e_coli",
    chromosomes=_chromosomes,
    default_genetic_map=None)


###########################################################
#
# Demographic models
#
###########################################################

# TODO add a generation time here
# default_generation_time = -1

class LapierreConstant(models.Model):
    """
    Model Name:
        LapierreConstant

    Model Description:
        The constant population size model from `Lapierre et al. 2016 <https://doi.org/
        10.1093/molbev/msw048>`_. The population does not undergo growth or size changes
        making it a simple scenario in which to study the effects of recombination and/
        or mutation on a variety of inference methods.

    Model population indexes:
        - E. coli: 0

    Parameter Table:
        .. csv-table::
            :widths: 15 8 20
            :header: "Parameter Type (units)", "Value", "Description"
            :file: ../docs/parameter_tables/e_coli/LapierreConstant_params.csv

    CLI help:
        python -m stdpopsim e-coli LapierreConstant -h

    Citation:
        Lapierre, M., Blin, C., Lambert, A., Achaz, G. & Rocha, E. P. C. The Impact of
        Selection, Gene Conversion, and Biased Sampling on the Assessment of Microbial
        Demography. Mol Biol Evol 33, 1711–1725 (2016).


    """  # noqa: E501

    author = "Lapierre et al."
    year = "2016"
    doi = "https://doi.org/10.1093/molbev/msw048"

    def __init__(self):
        super().__init__()

        N_e = 1.8e8
        # Single population
        self.population_configurations = [
            msprime.PopulationConfiguration(initial_size=N_e),
        ]
        self.migration_matrix = [[0]]
        self.demographic_events = []
