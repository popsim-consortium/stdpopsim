"""
Common infrastructure for specifying demographic models.
"""
import sys

import attr
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


def sampling_times_equal(
        populations1, populations2, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL):
    """
    Returns True if the sampling times for the specified lists of Populations
    objects are equal to the specified tolerances.
    """
    try:
        verify_sampling_times_equal(
            populations1, populations2, rtol=rtol, atol=atol)
        return True
    except UnequalModelsError:
        return False


def verify_sampling_times_equal(
        populations1, populations2, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL):
    """
    Check the the sampling times for lists of population objects are the same
    """
    # Retrieve sampling times
    sampling_times1 = np.array([p.sampling_time for p in populations1])
    sampling_times2 = np.array([p.sampling_time for p in populations2])
    # Check for equal vector length
    if len(sampling_times1) != len(sampling_times2):
        raise UnequalModelsError("Different numbers of populations")
    # Get indicies where None values present and compare
    s1_none = np.where(sampling_times1 == None) # noqa
    s2_none = np.where(sampling_times2 == None) # noqa
    if not np.array_equal(s1_none, s2_none):
        raise UnequalModelsError("None-valued sampling times differ")
    # Subset out only non-None values and cast to float so they can be compared
    s1_subset = sampling_times1[sampling_times1 != np.array(None)].astype(float)
    s2_subset = sampling_times2[sampling_times2 != np.array(None)].astype(float)
    if not np.allclose(s1_subset, s2_subset, rtol=rtol, atol=atol):
        raise UnequalModelsError("Positive real-valued sampling times differ")


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

    :ivar ~.id: The id of the population.
    :vartype ~.id: str
    :ivar ~.description: a short description of the population
    :vartype ~.description: str
    :ivar sampling_time: an integer value which indicates how many
        generations prior to the present individuals should samples should
        be drawn from this population. If `None`, sampling not allowed from this
        population (default = 0).
    :vartype sampling_time: int
    """
    def __init__(self, id, description, sampling_time=0):
        self.id = id
        self.description = description
        self.sampling_time = sampling_time

    @property
    def allow_samples(self):
        return self.sampling_time is not None

    def asdict(self):
        """
        Returns a dictionary representing the metadata about this population.
        """
        return {"id": self.id, "description": self.description,
                "sampling_time": self.sampling_time}


@attr.s(kw_only=True)
class DemographicModel(object):
    """
    Class representing a demographic model.

    This class is indended to be used by model implementors. To instead
    obtain a pre-specified model, see :class:`Species.get_demographic_model`.

    :ivar ~.id: The unique identifier for this model. DemographicModel IDs should be
        short and memorable, perhaps as an abbreviation of the model's name.
    :vartype ~.id: str
    :ivar ~.description: A short description of this model as it would be used in
        written text, e.g., "Three population Out-of-Africa". This should
        describe the model itself and not contain author or year information.
    :vartype ~.description: str
    :ivar long_description: A concise, but detailed, summary of the model.
    :vartype long_description: str
    :ivar generation_time: Mean inter-generation interval, in years.
    :vartype generation_time: int
    :ivar populations: TODO
    :vartype populations: list of :class:`.Population`
    :ivar qc_model: An independent implementation of the model, against which
        the model's accuracy is validated. This should not be set by the user,
        and may be None if no QC implementation exists yet.
    :vartype qc_model: :class:`.DemographicModel` or None

    :ivar citations: TODO
    :vartype citations: list of :class:`.Citation`
    :ivar demographic_events: TODO
    :vartype demographic_events: list of :class:`msprime.DemographicEvent`
    :ivar population_configurations: TODO
    :vartype population_configurations: list of :class:`msprime.PopulationConfiguration`
    :ivar migration_matrix: TODO
    :vartype migration_matrix: list of list of int
    """

    # required attributes
    id = attr.ib(type=str)
    description = attr.ib(type=str)
    long_description = attr.ib(type=str)
    generation_time = attr.ib(type=int)

    # optional attributes
    citations = attr.ib(factory=list)
    demographic_events = attr.ib(factory=list)
    population_configurations = attr.ib(factory=list)
    migration_matrix = attr.ib()
    populations = attr.ib()

    qc_model = attr.ib(default=None)

    @populations.default
    def _populations_default(self):
        npops = len(self.population_configurations)
        pops = []
        for i in range(npops):
            # Try to get population constructor info from PopConfig metadata
            kwargs = self.population_configurations[i].metadata
            # This is here so we don't have to change msprime, may remove later
            if kwargs is None:
                kwargs = dict()
            if "id" not in kwargs:
                kwargs["id"] = f"pop{i}"
            if "description" not in kwargs:
                kwargs["description"] = f"Population {i}"
            pops.append(Population(**kwargs))
        return(pops)

    @migration_matrix.default
    def _migration_matrix_default(self):
        npops = len(self.population_configurations)
        return([[0 for j in range(npops)] for i in range(npops)])

    @property
    def num_populations(self):
        return len(self.populations)

    @property
    def num_sampling_populations(self):
        return sum(int(pop.allow_samples) for pop in self.populations)

    @staticmethod
    def empty(**kwargs):
        """
        Return a model with the mandatory attributes filled out.
        """
        kwargs.update(
                id=kwargs.get("id", ""),
                description=kwargs.get("description", ""),
                long_description=kwargs.get("long_description", ""),
                populations=kwargs.get("populations", []),
                generation_time=kwargs.get("generation_time", -1))
        return DemographicModel(**kwargs)

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
        verify_sampling_times_equal(
            self.populations, other.populations,
            rtol=rtol, atol=atol)

    def register_qc(self, qc_model):
        """
        Register a QC model implementation for this model.
        """
        if not isinstance(qc_model, self.__class__):
            raise ValueError(
                    f"Cannot register non-DemographicModel '{qc_model}' as QC model")
        if self.qc_model is not None:
            raise ValueError(f"QC model already registered for {self.id}.")
        self.qc_model = qc_model

    def debug(self, out_file=sys.stdout):
        # Use the demography debugger to print out the demographic history
        # that we have just described.
        dd = self.get_demography_debugger()
        dd.print_history(out_file)

    def get_samples(self, *args):
        """
        Returns a list of msprime.Sample objects, with the number of samples
        from each population determined by the positional arguments.
        For instance, ``model.get_samples(2, 5, 7)`` would return a list of 14 samples,
        two of which are from the model's first population (i.e., with population ID
        ``model.populations[0].id``), five are from the model's second population,
        and seven are from the model's third population.
        The number of of arguments must be less than or equal to the number of
        "sampling" populations, ``model.num_sampling_populations``;
        if the number of arguments is less than the number of sampling populations,
        then remaining numbers are treated as zero.
        """
        samples = []
        for pop_index, n in enumerate(args):
            if self.populations[pop_index].allow_samples:
                sample = msprime.Sample(
                                        pop_index,
                                        time=self.populations[pop_index].sampling_time)
                samples.extend([sample] * n)
            elif n > 0:
                raise ValueError("Samples requested from non-sampling population"
                                 f" {pop_index}")
        return samples

    def get_demography_debugger(self):
        """
        Returns an :class:`msprime.DemographyDebugger` instance initialized
        with the parameters for this model. Please see the msprime documentation
        for details on how to use a DemographyDebugger.

        :return: A DemographyDebugger instance for this DemographicModel.
        :rtype: msprime.DemographyDebugger
        """
        ddb = msprime.DemographyDebugger(
            population_configurations=self.population_configurations,
            migration_matrix=self.migration_matrix,
            demographic_events=self.demographic_events)
        return ddb


# Reusable generic populations
_pop0 = Population(id="pop0", description="Generic population")
_pop1 = Population(id="pop1", description="Generic population")
_popAnc = Population(id="popAnc", description="Generic ancestral population",
                     sampling_time=None)


class PiecewiseConstantSize(DemographicModel):
    """
    Class representing a generic simulation model that can be run to output a
    tree sequence. This is a piecewise constant size model, which allows for
    instantaneous population size change over multiple epochs in a single population.

    :ivar N0: The initial effective population size
    :vartype N0: float
    :ivar args: Each subsequent argument is a tuple (t, N) which gives the
        time at which the size change takes place and the population size.

    The usage is best illustrated by an example:

    .. code-block:: python

        model1 = stdpopsim.PiecewiseConstantSize(N0, (t1, N1)) # One change
        model2 = stdpopsim.PiecewiseConstantSize(N0, (t1, N1), (t2, N2)) # Two changes
    """

    id = "PiecewiseConstant"
    description = "Piecewise constant size population model over multiple epochs."
    citations = []
    populations = [_pop0]
    author = None
    year = None
    doi = None
    generation_time = 1

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


class IsolationWithMigration(DemographicModel):
    """
    Class representing a generic simulation model that can be run to output a tree
    sequence. A generic isolation with migration model where a single ancestral
    population of size NA splits into two populations of constant size N1
    and N2 time T generations ago, with migration rates M12 and M21 between
    the split populations. Sampling is disallowed in population index 0,
    as this is the ancestral population.

    :ivar NA: The initial ancestral effective population size
    :vartype NA: float
    :ivar N1: The effective population size of population 1
    :vartype N1: float
    :ivar N2: The effective population size of population 2
    :vartype N2: float
    :ivar T: Time of split between populations 1 and 2 (in generations)
    :vartype T: float
    :ivar M12: Migration rate from population 1 to 2
    :vartype M12: float
    :ivar M21: Migration rate from population 2 to 1
    :vartype M21: float


    Example usage:

    .. code-block:: python

        model1 = stdpopsim.IsolationWithMigration(NA, N1, N2, T, M12, M21)

    """
    id = "IsolationWithMigration"
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
    generation_time = 1

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
