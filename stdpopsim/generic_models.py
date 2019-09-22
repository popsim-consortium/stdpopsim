# generic_models
import msprime


#########
# Generic Model Mixin class
# provides generic models that can be
# initialized with methods
#
#
########


class ConstantSizeMixin(object):
    def __init__(self, N):
        """
        Model Name:
            ConstantSizeMixin

        Model Description:
            Generic model of constant size. Uses default_population_size
            depending on organism

        Model population indexes:
            - Population 0: 0

        CLI help:
            python -m stdpopsim homo-sapiens GenericConstantSize -h
        """  # noqa: E501
        self.population_configurations = [
            msprime.PopulationConfiguration(initial_size=N)
        ]
        self.migration_matrix = [[0]]
        self.demographic_events = []


class TwoEpochMixin(object):
    def __init__(self, N1, N2, t):
        """
        Model Name:
            TwoEpochMixin

        Model Description:
            Generic model of a single population with
            piecewise constant population size with a single change.
            Defaults depend on organism, but population sizes and time
            should be set by user

        Model population indexes:
            - Population 0: 0

        CLI help:
            python -m stdpopsim homo-sapiens GenericTwoEpoch -h
        """  # noqa: E501
        self.population_configurations = [
            msprime.PopulationConfiguration(initial_size=N1)
        ]
        self.migration_matrix = [[0]]
        self.demographic_events = [
            msprime.PopulationParametersChange(
                time=t, initial_size=N2, growth_rate=0, population_id=0),
        ]
