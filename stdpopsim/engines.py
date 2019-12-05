import logging

import attr
import msprime
import stdpopsim

logger = logging.getLogger(__name__)

_registered_engines = {}


def register_engine(engine):
    """
    Registers the specified simulation engine.
    """
    logger.debug(f"Registering simulation engine '{engine.id}'")
    _registered_engines[engine.id] = engine


def get_engine(id):
    """
    Returns the simulation engine with the specified id.
    """
    if id not in _registered_engines:
        raise ValueError(f"Simulation engine '{id}' not registered")
    return _registered_engines[id]


def all_engines():
    """
    Returns an iterator over all registered simulation engines.
    """
    for sim_engine in _registered_engines.values():
        yield sim_engine


@attr.s(frozen=True)
class Engine(object):
    """
    Abstract class representing a simulation engine.

    To implement a new simulation engine, one should inherit from this
    class. At a minimum, the ``id``, ``name`` and ``citations`` attributes must
    be set, and the :func:`simulate()` and :func:`get_version()` methods must
    be implemented. See msprime example in ``engines.py``.

    :ivar id: The unique identifier for the simulation engine.
    :vartype id: str
    :ivar name: The name for this engine as it would be used in written text
        such as the CLI or error messages.
    :vartype name: str
    :ivar citations: A list of citations for the simulation engine.
    :vartype citations: list of :class:`.Citation`
    """

    def simulate(self, model=None, contig=None, samples=None, **kwargs):
        """
        Simulates the model for the specified contig and samples.

        :param model: The demographic model to simulate.
        :type model: :class:`.Model`
        :param contig: The contig, defining the length and recombination
            rate(s).
        :type contig: :class:`msprime.simulations.Contig`
        :param samples: The samples to be obtained from the simulation.
        :type samples: list of :class:`msprime.simulations.Sample`
        :return: A succinct tree sequence.
        :rtype: :class:`tskit.trees.TreeSequence`
        """
        raise NotImplementedError()

    def get_version(self):
        """
        Returns the version of the engine.

        :rtype: str
        """
        raise NotImplementedError()

    def add_arguments(self, parser):
        """
        Defines engine specific command line parameters.

        If implemented, this function should call :func:`parser.add_argument()`
        to define engine specific command line parameters. The engine's
        :func:`simulate()` method should accept keyword arguments matching the
        names defined here.

        :param parser: A CLI argument parser group.
        :type parser: :class:`argparse._ArgumentGroup`
        """
        pass


class _MsprimeEngine(Engine):
    id = "msprime"
    name = "msprime"
    citations = [
            stdpopsim.Citation(
                doi="https://doi.org/10.1371/journal.pcbi.1004842",
                year="2016",
                author="Kelleher et al."),
            ]

    def simulate(self, model=None, contig=None, samples=None, seed=None,
                 **kwargs):
        ts = msprime.simulate(
                samples=samples,
                recombination_map=contig.recombination_map,
                population_configurations=model.population_configurations,
                migration_matrix=model.migration_matrix,
                demographic_events=model.demographic_events,
                random_seed=seed)
        mutation_model = msprime.InfiniteSites(msprime.NUCLEOTIDES)
        return msprime.mutate(
            ts, rate=contig.mutation_rate, model=mutation_model,
            random_seed=seed)

    def get_version(self):
        return msprime.__version__


register_engine(_MsprimeEngine())


def get_default_engine():
    """
    Returns the default simulation engine (msprime).
    """
    return get_engine("msprime")
