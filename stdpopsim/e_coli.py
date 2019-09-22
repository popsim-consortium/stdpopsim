"""
Genome and demographic model definitions for Escherichia coli.
"""
import msprime

import stdpopsim.models as models
import stdpopsim.genomes as genomes
import stdpopsim.generic_models as generic_models

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
default_generation_time = 0.00003805175  # 1.0 / (525600 min/year / 20 min/gen)


class EColiModel(models.Model):
    """
    TODO: documentation
    """
    def __init__(self):
        super().__init__()
        self.generation_time = default_generation_time
        self.default_population_size = 10000


class GenericConstantSize(EColiModel, generic_models.ConstantSizeMixin):
    def __init__(self):
        EColiModel.__init__(self)
        generic_models.ConstantSizeMixin.__init__(self, self.default_population_size)


class GenericTwoEpoch(EColiModel, generic_models.TwoEpochMixin):
    def __init__(self, n2=None, t=None):
        EColiModel.__init__(self)
        n1 = self.default_population_size
        if n2 is None:
            n2 = n1 / 2.0
        if t is None:
            t = n1 / 100
        generic_models.TwoEpochMixin.__init__(self, n1, n2, t)


class LapierreConstant(EColiModel):
    name = "LapierreConstant"
    short_description = "Constant size model for E-coli"
    description = """
        The constant population size model from `Lapierre et al. 2016 <https://doi.org/
        10.1093/molbev/msw048>`_. The population does not undergo growth or size changes
        making it a simple scenario in which to study the effects of recombination and/
        or mutation on a variety of inference methods.
    """
    citations = [
        """
        Lapierre, M., Blin, C., Lambert, A., Achaz, G. & Rocha, E. P. C. The Impact of
        Selection, Gene Conversion, and Biased Sampling on the Assessment of Microbial
        Demography. Mol Biol Evol 33, 1711–1725 (2016).
        """
    ]
    populations = [
        models.Population(name="e_coli", description="Single E-coli population"),
    ]
    author = "Lapierre et al."
    year = "2016"
    doi = "https://doi.org/10.1093/molbev/msw048"

    def __init__(self):
        super().__init__()

        N_e = 1.8e8
        # Single population
        self.population_configurations = [
            msprime.PopulationConfiguration(
                initial_size=N_e, metadata=self.populations[0].asdict()),
        ]
        self.migration_matrix = [[0]]
        self.demographic_events = []
