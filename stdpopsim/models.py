"""
Common infrastructure for specifying demographic models.
"""
import sys
import inspect

import msprime
import numpy as np
import scipy.linalg

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

# ###############################################################
# UNDER CONSTRUCTION - GROUND TRUTH N_t
# ###############################################################

    def population_size_trajectory(self, end, num_steps=10):
        """
        This function computes the defined effective population size
        for each sub-population, for any demographic model.

        `end` specifies at which generation to stop
        collecting effective population sizes among subpopulations.
        The sampling will always begin at generation 0

        `num_steps` parameter will determine how many points from
        0 to `end` to sample.

        Returns a tuple where the first value
        is the steps (in generations) and the second is a
        num_steps-by-num_pops ndarray that will contain the effective
        population sizes for each subpopulation individually, for each step,
        as defined by the model which inherits from this class.
        """

        num_pops = self.num_pops()
        N_t = np.zeros([num_steps, num_pops])
        steps = np.linspace(0, end, num_steps)
        for j, t in enumerate(steps):
            N = self.pop_size_and_migration_at_t(t)["population_sizes"]
            N_t[j] = N
        return steps, N_t

    def coalescence_rate_trajectory(self, end, num_samples,
                                    num_steps=10, min_pop_size=1):
        """
        This function will calculate the ground truth
        coalescent rate trajectory, r (Ne = 1 / 2r), for a population with
        multiple sub-populations. This takes a
        demographic model and returns a
        function of time that gives haploid coalescence rate.

        `num_samples` should be a list of the same length as the number
        of populations, so that num_samples[j] is the number of
        samples in subpopulation j.

        `num_steps` parameter will determine how many points along the
        time axis are returned.

        `min_pop_size` is the smallest a population size can be represented as during
        the computation of coalescent rates. In models where populations
        grow exponentially, very small effective population sizes skew the
        calculation when we want to compute coalescent rates within a population (1/2Ne).

        Returns a tuple where the first value
        is the steps (in generations) and the second
        is an array whose jth element is the coalescence rate at time
        steps[j], i.e., the probability per unit time that a randomly
        chosen pair of lineages from the original samples coalesces
        around that time, conditioned on the pair not having yet coalesced.
        The third is the probablility that the pair of lineages
        has not yet coalesced.
        """

        num_pops = self.num_pops()
        assert(len(num_samples) == num_pops)
        P = np.zeros([num_pops**2, num_pops**2])
        IA = np.array(range(num_pops**2)).reshape([num_pops, num_pops])
        Identity = np.eye(num_pops)
        for x in range(num_pops):
            for y in range(num_pops):
                K_x = num_samples[x]
                K_y = num_samples[y]
                P[IA[x, y], IA[x, y]] = K_x * (K_y - (x == y))
        P = P / np.sum(P)
        r = np.zeros(num_steps)
        p_t = np.zeros(num_steps)
        steps = np.linspace(0, end, num_steps)
        dt = steps[1] - steps[0]
        mass_migration_objects = []
        mass_migration_steps = []
        for demo in self.demographic_events:
            if type(demo) == msprime.simulations.MassMigration:
                mass_migration_objects.append(demo)
                mass_migration_steps.append(int(demo.time/dt))
        for j in range(num_steps):
            N_M = self.pop_size_and_migration_at_t(steps[j])
            N = N_M["population_sizes"]
            M = N_M["migration_matrix"]
            C = np.zeros([num_pops**2, num_pops**2])
            for idx in range(num_pops):
                C[IA[idx, idx], IA[idx, idx]] = 1 / (2 * max(min_pop_size, N[idx]))
            for idx, row in enumerate(M):
                M[idx][idx] = -1 * sum(row)
            G = (np.kron(M, Identity) + np.kron(Identity, M)) - C
            if j in mass_migration_steps:
                idx = mass_migration_steps.index(j)
                ds = mass_migration_objects[idx].time - (j * dt)
                a = mass_migration_objects[idx].source
                b = mass_migration_objects[idx].dest
                p = mass_migration_objects[idx].proportion
                S = np.eye(num_pops**2, num_pops**2)
                for x in range(num_pops):
                    if x == a:
                        S[IA[a, a], IA[a, b]] = S[IA[a, a], IA[b, a]] = p * (1 - p)
                        S[IA[a, a], IA[b, b]] = p ** 2
                        S[IA[a, a], IA[a, a]] = (1 - p) ** 2
                    else:
                        S[IA[x, a], IA[x, b]] = S[IA[a, x], IA[b, x]] = p
                        S[IA[x, a], IA[x, a]] = S[IA[a, x], IA[a, x]] = 1 - p
                P = np.matmul(P, scipy.linalg.expm(ds * G))
                P = np.matmul(P, S)
                P = np.matmul(P, scipy.linalg.expm((dt - ds) * G))
            else:
                P = np.matmul(P, scipy.linalg.expm(dt * G))
            p_t[j] = np.sum(P)
            r[j] = np.sum(np.matmul(P, C)) / np.sum(P)
        return steps, r, p_t

    def num_pops(self):
        """
        This function returns the number populations
        defined by the demographic model
        """
        ddb = msprime.DemographyDebugger(**self.asdict())
        return len(ddb.epochs[0].populations)

    def pop_size_and_migration_at_t(self, t):
        """
        Given a time, t (in generations), find and return
        a dictionary containing:
        1: "population_sizes" the vector N which should represent effective
        population size for each populations at time t and,
        2: "migration_matrix" The migration matrix for for the populations at time t
        """

        ddb = msprime.DemographyDebugger(**self.asdict())
        epochs = ddb.epochs
        j = 0
        while epochs[j].end_time <= t:
            j += 1
        N = ddb.population_size_history[:, j]
        for i, pop in enumerate(epochs[j].populations):
            s = t - epochs[j].start_time
            g = pop.growth_rate
            N[i] *= np.exp(-1 * g * s)
        return {"population_sizes": N, "migration_matrix": epochs[j].migration_matrix}

# ###############################################################
# END CONSTRUCTION - GROUND TRUTH N_t
# ###############################################################

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


def all_models():
    """
    Returns the list of all Model classes that are defined within the stdpopsim
    module.
    """
    ret = []
    for cls in Model.__subclasses__():
        mod = inspect.getmodule(cls).__name__
        if mod.startswith("stdpopsim"):
            ret.append(cls())
    return ret
