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
    if engine.id in _registered_engines:
        raise ValueError(f"Simulation engine '{engine.id}' already registered.")
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
    class. At a minimum, the ``id``, ``description``
    and ``citations`` attributes must
    be set, and the :func:`simulate()` and :func:`get_version()` methods must
    be implemented. See msprime example in ``engines.py``.

    :ivar id: The unique identifier for the simulation engine.
    :vartype id: str
    :ivar description: A short description of this engine.
    :vartype description: str
    :ivar citations: A list of citations for the simulation engine.
    :vartype citations: list of :class:`.Citation`
    """

    def simulate(self, model=None, contig=None, samples=None, **kwargs):
        """
        Simulates the model for the specified contig and samples.

        :param model: The demographic model to simulate.
        :type model: :class:`.DemographicModel`
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
    description = "Msprime coalescent simulator"
    citations = [
            stdpopsim.Citation(
                doi="https://doi.org/10.1371/journal.pcbi.1004842",
                year="2016",
                author="Kelleher et al.",
                reasons={stdpopsim.CiteReason.ENGINE}),
            ]

    def simulate(self, demographic_model=None, contig=None, samples=None, seed=None,
                 **kwargs):
        return msprime.simulate(
                samples=samples,
                recombination_map=contig.recombination_map,
                mutation_rate=contig.mutation_rate,
                population_configurations=demographic_model.population_configurations,
                migration_matrix=demographic_model.migration_matrix,
                demographic_events=demographic_model.demographic_events,
                random_seed=seed)

    def get_version(self):
        return msprime.__version__


register_engine(_MsprimeEngine())


def get_default_engine():
    """
    Returns the default simulation engine (msprime).
    """
    return get_engine("msprime")


class _MsprimeDTWFEngine(Engine):
    id = "msprime-dtwf"
    description = "Msprime discrete-time Wright Fisher simulator"
    citations = [
            stdpopsim.Citation(
                # biorxiv; update upon publication
                doi="https://doi.org/10.1101/674440",
                year="2019",
                author="Nelson et al.",
                reasons={stdpopsim.CiteReason.ENGINE}),
            ]

    def add_arguments(self, parser):
        parser.add_argument(
                "--time-hudson", metavar="GEN", type=int, default=None,
                help="Time (generations in the past) at which the simulation "
                     "switches from msprime's DTWF simulation model to a "
                     "coalescent-with-recombination simulation model "
                     "(Hudson's algorithm)")

    def simulate(
            self, demographic_model=None, contig=None, samples=None, seed=None,
            time_hudson=None, **kwargs):

        demographic_events = demographic_model.demographic_events.copy()
        if time_hudson is not None:
            model_change = msprime.SimulationModelChange(time_hudson, "hudson")
            demographic_events.append(model_change)
            demographic_events.sort(key=lambda x: x.time)

        return msprime.simulate(
                samples=samples,
                recombination_map=contig.recombination_map,
                mutation_rate=contig.mutation_rate,
                population_configurations=demographic_model.population_configurations,
                migration_matrix=demographic_model.migration_matrix,
                demographic_events=demographic_events,
                random_seed=seed,
                model="dtwf")

    def get_version(self):
        return msprime.__version__


register_engine(_MsprimeDTWFEngine())
