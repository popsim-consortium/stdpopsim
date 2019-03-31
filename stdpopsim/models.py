"""
Common infrastructure for specifying demographic models.
"""
import sys

import msprime
import numpy as np


# Defaults taken from np.allclose
DEFAULT_ATOL = 1e-05
DEFAULT_RTOL = 1e-08


def population_configurations_equal(
        pop_configs1, pop_configs2, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL):
    """
    Returns True if the specified lists of msprime PopulationConfiguration
    objects are equal to the specified tolerances.

    We make some assumptions here to ensure that the models we specify
    are well-defined: (1) The sample size is not set for PopulationConfigurations
    (2) the initial_size is defined.
    """
    for pc1, pc2 in zip(pop_configs1, pop_configs2):
        if pc1.sample_size is not None or pc2.sample_size is not None:
            raise ValueError(
                "Models defined in stdpopsim must not use the 'sample_size' "
                "PopulationConfiguration option")
        if pc1.initial_size is None or pc2.initial_size is None:
            raise ValueError(
                "Models defined in stdpopsim must set the initial_size")
    sample_size1 = np.array([pc.sample_size for pc in pop_configs1])
    sample_size2 = np.array([pc.sample_size for pc in pop_configs2])
    initial_size1 = np.array([pc.initial_size for pc in pop_configs1])
    initial_size2 = np.array([pc.initial_size for pc in pop_configs2])
    growth_rate1 = np.array([pc.growth_rate for pc in pop_configs1])
    growth_rate2 = np.array([pc.growth_rate for pc in pop_configs2])
    return (
        len(pop_configs1) == len(pop_configs2)
        and np.all(sample_size1 == sample_size2)
        and np.allclose(initial_size1, initial_size2, rtol=rtol, atol=atol)
        and np.allclose(growth_rate1, growth_rate2, rtol=rtol, atol=atol))


def demographic_events_equal(
        events1, events2, num_populations, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL):
    """
    Returns True if the specified list of msprime DemographicEvent objects are equal
    to the specified tolerances.
    """
    # Get the low-level dictionary representations of the events.
    dicts1 = [event.get_ll_representation(num_populations) for event in events1]
    dicts2 = [event.get_ll_representation(num_populations) for event in events2]
    if len(dicts1) != len(dicts2):
        return False
    for d1, d2 in zip(dicts1, dicts2):
        if set(d1.keys()) != set(d2.keys()):
            return False
        for key in d1.keys():
            value1 = d1[key]
            value2 = d2[key]
            if isinstance(value1, float):
                if not np.isclose(value1, value2, rtol=rtol, atol=atol):
                    return False
            else:
                if value1 != value2:
                    return False
    return True


class Model(object):
    """
    Class representing a simulation model that can be run in msprime.

    .. todo:: Document this class.
    """
    def __init__(self):
        self.population_configurations = []
        self.demographic_events = []
        # Defaults to a single population
        self.migration_matrix = [[0]]

    def debug(self, out_file=sys.stdout):
        # Use the demography debugger to print out the demographic history
        # that we have just described.
        dd = msprime.DemographyDebugger(
            population_configurations=self.population_configurations,
            migration_matrix=self.migration_matrix,
            demographic_events=self.demographic_events)
        dd.print_history(out_file)

    def asdict(self):
        return {
            "population_configurations": self.population_configurations,
            "migration_matrix": self.migration_matrix,
            "demographic_events": self.demographic_events}

    def equals(self, other, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL):
        """
        Returns True if this model is equal to the specified model to the
        specified numerical tolerance (as defined by numpy.allclose).

        We use the 'equals' method here rather than the equality operator
        because we need to be able to specifiy the numerical tolerances.
        """
        ret = False
        try:
            mm1 = np.array(self.migration_matrix)
            mm2 = np.array(other.migration_matrix)
            ret = (
                mm1.shape == mm2.shape
                and np.allclose(mm1, mm2, rtol=rtol, atol=atol)
                and population_configurations_equal(
                    self.population_configurations, other.population_configurations,
                    rtol=rtol, atol=atol)
                and demographic_events_equal(
                    self.demographic_events, other.demographic_events,
                    len(self.population_configurations),
                    rtol=rtol, atol=atol))
        except AttributeError:
            # Anything that's not duck-typeable to a Model is considered not equal
            pass
        return ret


def all_models():
    """
    Returns the list of all Model classes that have been defined.
    """
    return [cls() for cls in Model.__subclasses__()]
