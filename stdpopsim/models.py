"""
Common infrastructure for specifying demographic models.
"""
import sys

import msprime
import numpy as np


# Defaults taken from np.allclose
DEFAULT_ATOL = 1e-05
DEFAULT_RTOL = 1e-08


class UnequalModelsError(Exception):
    """
    Exception raised models by verify_equal to indicate that models are
    not sufficiently close.
    """


def population_configurations_equal(
        pop_configs1, pop_configs2, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL):
    """
    Returns True if the specified lists of msprime PopulationConfiguration
    objects are equal to the specified tolerances.

    See the :func:`.verify_population_configurations_equal` function for
    details on the assumptions made about the objects.
    """
    try:
        verify_population_configurations_equal(
            pop_configs1, pop_configs2, rtol=rtol, atol=atol)
        return True
    except UnequalModelsError:
        return False


def verify_population_configurations_equal(
        pop_configs1, pop_configs2, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL):
    """
    Checks if the specified lists of msprime PopulationConfiguration
    objects are equal to the specified tolerances and raises an UnequalModelsError
    otherwise.

    We make some assumptions here to ensure that the models we specify
    are well-defined: (1) The sample size is not set for PopulationConfigurations
    (2) the initial_size is defined. If these assumptions are violated a
    ValueError is raised.
    """
    for pc1, pc2 in zip(pop_configs1, pop_configs2):
        if pc1.sample_size is not None or pc2.sample_size is not None:
            raise ValueError(
                "Models defined in stdpopsim must not use the 'sample_size' "
                "PopulationConfiguration option")
        if pc1.initial_size is None or pc2.initial_size is None:
            raise ValueError(
                "Models defined in stdpopsim must set the initial_size")
    if len(pop_configs1) != len(pop_configs2):
        raise UnequalModelsError("Different numbers of populations")
    initial_size1 = np.array([pc.initial_size for pc in pop_configs1])
    initial_size2 = np.array([pc.initial_size for pc in pop_configs2])
    if not np.allclose(initial_size1, initial_size2, rtol=rtol, atol=atol):
        raise UnequalModelsError("Initial sizes differ")
    growth_rate1 = np.array([pc.growth_rate for pc in pop_configs1])
    growth_rate2 = np.array([pc.growth_rate for pc in pop_configs2])
    if not np.allclose(growth_rate1, growth_rate2, rtol=rtol, atol=atol):
        raise UnequalModelsError("Growth rates differ")


def demographic_events_equal(
        events1, events2, num_populations, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL):
    """
    Returns True if the specified list of msprime DemographicEvent objects are equal
    to the specified tolerances.
    """
    try:
        verify_demographic_events_equal(
            events1, events2, num_populations, rtol=rtol, atol=atol)
        return True
    except UnequalModelsError:
        return False


def verify_demographic_events_equal(
        events1, events2, num_populations, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL):
    """
    Checks if the specified list of msprime DemographicEvent objects are equal
    to the specified tolerances and raises a UnequalModelsError otherwise.
    """
    # Get the low-level dictionary representations of the events.
    dicts1 = [event.get_ll_representation(num_populations) for event in events1]
    dicts2 = [event.get_ll_representation(num_populations) for event in events2]
    if len(dicts1) != len(dicts2):
        raise UnequalModelsError("Different numbers of demographic events")
    for d1, d2 in zip(dicts1, dicts2):
        if set(d1.keys()) != set(d2.keys()):
            raise UnequalModelsError("Different types of demographic events")
        for key in d1.keys():
            value1 = d1[key]
            value2 = d2[key]
            if isinstance(value1, float):
                if not np.isclose(value1, value2, rtol=rtol, atol=atol):
                    raise UnequalModelsError("Event {} mismatch: {} != {}".format(
                        key, value1, value2))
            else:
                if value1 != value2:
                    raise UnequalModelsError("Event {} mismatch: {} != {}".format(
                        key, value1, value2))


class Population(object):
    """
    Class recording metadata representing a population in a simulation.
    """
    # TODO change this to use the usual id, name combination
    def __init__(self, name, description, allow_samples=True):
        self.name = name
        self.description = description
        self.allow_samples = allow_samples

    def asdict(self):
        """
        Returns a dictionary representing the metadata about this population.
        """
        return {"name": self.name, "description": self.description}


class Model(object):
    """
    Class representating a simulation model that can be run to output a tree sequence.
    Concrete subclasses must define population_configurations, demographic_events
    and migration_matrix instance variables which define the model.

    :ivar id: The unique identifier for this model. Model IDs should be
        short and memorable, perhaps as an abbreviation of the model's
        name.
    :vartype id: str
    :ivar name: The informal name for this model as it would be used in
        written text, e.g., "Three population Out-of-Africa"
    :vartype informal_name: str
    """
    # TODO the infrastructure here is left over from a structure that
    # rigidly used class definitions as a way to define population
    # models. This contructor should take all the instance variables
    # as parameteters, and we should use factory functions to define
    # the model instances that are added to the catalog rather than
    # subclasses.

    def __init__(self):
        self.population_configurations = []
        self.demographic_events = []
        # Defaults to a single population
        self.migration_matrix = [[0]]
        self.generation_time = -1

    @property
    def num_populations(self):
        return len(self.populations)

    @property
    def num_sampling_populations(self):
        return sum(int(pop.allow_samples) for pop in self.populations)

    def equals(self, other, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL):
        """
        Returns True if this model is equal to the specified model to the
        specified numerical tolerance (as defined by numpy.allclose).

        We use the 'equals' method here rather than the equality operator
        because we need to be able to specifiy the numerical tolerances.
        """
        try:
            self.verify_equal(other, rtol=rtol, atol=atol)
            return True
        except (UnequalModelsError, AttributeError):
            return False

    def verify_equal(self, other, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL):
        """
        Equivalent to the :func:`.equals` method, but raises a UnequalModelsError if the
        models are not equal rather than returning False.
        """
        mm1 = np.array(self.migration_matrix)
        mm2 = np.array(other.migration_matrix)
        if mm1.shape != mm2.shape:
            raise UnequalModelsError("Migration matrices different shapes")
        if not np.allclose(mm1, mm2, rtol=rtol, atol=atol):
            raise UnequalModelsError("Migration matrices differ")
        verify_population_configurations_equal(
            self.population_configurations, other.population_configurations,
            rtol=rtol, atol=atol)
        verify_demographic_events_equal(
            self.demographic_events, other.demographic_events,
            len(self.population_configurations),
            rtol=rtol, atol=atol)

    def debug(self, out_file=sys.stdout):
        # Use the demography debugger to print out the demographic history
        # that we have just described.
        dd = msprime.DemographyDebugger(
            population_configurations=self.population_configurations,
            migration_matrix=self.migration_matrix,
            demographic_events=self.demographic_events)
        dd.print_history(out_file)

    def get_samples(self, *args):
        """
        Returns a list of msprime.Sample objects as described by the args and
        keyword args. Positional arguments are interpreted as the number of
        samples to take from the given population.

        .. todo:: Add a description how the positional arguments work and
            perhaps link into a section of the tutorial showing it in action.

        """
        samples = []
        for pop_index, n in enumerate(args):
            samples.extend([msprime.Sample(pop_index, time=0)] * n)
        return samples

    def simulate(self, contig, samples, seed=None):
        """
        Simulates this model for the specified contig (defining the recombination
        map and mutation rate) and samples.
        """
        ts = msprime.simulate(
            samples=samples,
            recombination_map=contig.recombination_map,
            mutation_rate=contig.mutation_rate,
            population_configurations=self.population_configurations,
            migration_matrix=self.migration_matrix,
            demographic_events=self.demographic_events,
            random_seed=seed)
        return ts


# Reusable generic populations
_pop0 = Population(name="pop0", description="Generic population")
_pop1 = Population(name="pop1", description="Generic population")
_popAnc = Population(name="popAnc", description="Generic ancestral population",
                     allow_samples=False)


class PiecewiseConstantSize(Model):
    id = "constant"
    name = "Piecewise constant size"
    description = "Piecewise constant size population model over multiple epochs."
    citations = []
    populations = [_pop0]
    author = None
    year = None
    doi = None

    def __init__(self, N0, *args):
        self.population_configurations = [
            msprime.PopulationConfiguration(
                initial_size=N0, metadata=self.populations[0].asdict())
        ]
        self.migration_matrix = [[0]]
        self.demographic_events = []
        for t, N in args:
            self.demographic_events.append(msprime.PopulationParametersChange(
                time=t, initial_size=N, growth_rate=0, population_id=0))


class GenericIM(Model):
    """
    Takes parameters N0, N1, N2, T, M12, and M21, as list of arguments:
    (NA, N1, N2, T, M12, M21).
    Sampling is disallowed in population index 0, as this is the ancestral
    population.
    """
    id = "IM"
    name = "Isolation with migration"
    description = """
        A generic isolation with migration model where a single ancestral
        population of size NA splits into two populations of constant size N1
        and N2 time T generations ago, with migration rates M12 and M21 between
        the split populations.
        """
    citations = []
    populations = [_pop0, _pop1, _popAnc]
    author = None
    year = None
    doi = None

    def __init__(self, NA, N1, N2, T, M12, M21):
        self.population_configurations = [
            msprime.PopulationConfiguration(
                initial_size=N1, metadata=self.populations[0].asdict()),
            msprime.PopulationConfiguration(
                initial_size=N2, metadata=self.populations[1].asdict()),
            msprime.PopulationConfiguration(
                initial_size=NA, metadata=self.populations[2].asdict())
        ]
        self.migration_matrix = [[0, M12, 0], [M21, 0, 0], [0, 0, 0]]
        self.demographic_events = [
            msprime.MassMigration(
                time=T, source=0, destination=2, proportion=1),
            msprime.MassMigration(
                time=T, source=1, destination=2, proportion=1)
        ]
