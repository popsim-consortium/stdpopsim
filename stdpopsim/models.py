"""
Common infrastructure for specifying demographic models.
"""
import sys
import inspect
import string
import pathlib

import msprime
import numpy as np

import stdpopsim.citations as citations


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


_model_docstring_template = string.Template("""
Name:
    $name

Description:
    $description

Populations:
$populations

CLI help:
    $cli_help

Citations:
$citations

Parameter Table:
    .. csv-table::
        :widths: 15 8 20
        :header: "Parameter Type (units)", "Value", "Description"
        :file: $parameters_csv_file
""")


class Population(object):
    """
    Class recording metadata representing a population in a simulation.
    """
    def __init__(self, name, description, allow_samples=True):
        self.name = name
        self.description = description
        self.allow_samples = allow_samples

    def asdict(self):
        """
        Returns a dictionary representing the metadata about this population.
        """
        return {"name": self.name, "description": self.description}


def cleanup(docs):
    """
    Remove all newlines and trailing whitespace from the specified string.
    """
    return ' '.join([line.strip() for line in docs.strip().splitlines() if line.strip()])


class Model(citations.CitableMixin):
    """
    Class representating a simulation model that can be run to output a tree sequence.
    Concrete subclasses must define population_configurations, demographic_events
    and migration_matrix instance variables which define the model.

    Docstrings are handled in a non-standard way by these subclasses. The docstrings
    are used to generate online documentation and are composed of smaller pieces
    of information that are also used in the CLI documentation. To avoid repeating
    these smaller elements, we define them separately in each subclass, and
    call <SubclassName>._write_docstring() immediately after it is defined to
    build the docstring folling the standard template.
    """

    @classmethod
    def _write_docstring(cls):
        species = "homo_sapiens"
        base_dir = pathlib.Path(__file__).resolve().parents[1]
        parameters_csv_file = (
            base_dir / "docs" / "parameter_tables" / species / f"{cls.name}.csv")
        cli_help = f"python -m stdpopsim {species} {cls.name} -h"
        citations = "".join(
            f"    - {cleanup(citation)}\n" for citation in cls.citations)
        populations = "".join(
            f"    - {index}: {population.name}: {cleanup(population.description)}\n"
            for index, population in enumerate(cls.populations))
        cls.__doc__ = _model_docstring_template.substitute(
            name=cleanup(cls.name),
            description=cleanup(cls.description), populations=populations,
            cli_help=cli_help, citations=citations,
            parameters_csv_file=str(parameters_csv_file))

    def __init__(self):
        self.population_configurations = []
        self.demographic_events = []
        # Defaults to a single population
        self.migration_matrix = [[0]]
        self.generation_time = -1

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

    # TODO deprecate/remove this, as the recommended interface is now 'run'.
    def asdict(self):
        return {
            "population_configurations": self.population_configurations,
            "migration_matrix": self.migration_matrix,
            "demographic_events": self.demographic_events}

    def run(self, chromosome, samples):
        """
        Runs this model for the specified chromosome (defining the recombination
        map and mutation rate) and samples.
        """
        ts = msprime.simulate(
            samples=samples,
            recombination_map=chromosome.recombination_map,
            mutation_rate=chromosome.mutation_rate,
            population_configurations=self.population_configurations,
            migration_matrix=self.migration_matrix,
            demographic_events=self.demographic_events)
        return ts


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
