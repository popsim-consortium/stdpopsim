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


@attr.s
class Engine(object):
    """
    Abstract class representing a simulation engine.

    To implement a new simulation engine, one should inherit from this
    class. At a minimum, the ``id``, ``description``
    and ``citations`` attributes must
    be set, and the :func:`simulate()` and :func:`get_version()` methods must
    be implemented. See msprime example in ``engines.py``.

    :ivar ~.id: The unique identifier for the simulation engine.
    :vartype ~.id: str
    :ivar ~.description: A short description of this engine.
    :vartype ~.description: str
    :ivar citations: A list of citations for the simulation engine.
    :vartype citations: list of :class:`.Citation`
    """

    def simulate(
            self, demographic_model=None, contig=None, samples=None, seed=None,
            dry_run=False):
        """
        Simulates the model for the specified contig and samples.

        :param demographic_model: The demographic model to simulate.
        :type demographic_model: :class:`.DemographicModel`
        :param contig: The contig, defining the length and recombination
            rate(s).
        :type contig: :class:`msprime.simulations.Contig`
        :param samples: The samples to be obtained from the simulation.
        :type samples: list of :class:`msprime.simulations.Sample`
        :param seed: The seed for the random number generator.
        :type seed: int
        :param dry_run: If True, the simulation engine will return None without
            running the simulation.
        :type dry_run: bool
        :return: A succinct tree sequence.
        :rtype: :class:`tskit.trees.TreeSequence` or None
        """
        raise NotImplementedError()

    def get_version(self):
        """
        Returns the version of the engine.

        :rtype: str
        """
        raise NotImplementedError()


class _MsprimeEngine(Engine):
    id = "msprime"  #:
    description = "Msprime coalescent simulator"  #:
    citations = [
            stdpopsim.Citation(
                doi="https://doi.org/10.1371/journal.pcbi.1004842",
                year="2016",
                author="Kelleher et al.",
                reasons={stdpopsim.CiteReason.ENGINE}),
            ]
    # We default to the first model in the list.
    supported_models = ["hudson", "dtwf", "smc", "smc_prime"]
    model_citations = {"dtwf": [
             stdpopsim.Citation(
                 # biorxiv; update upon publication
                 doi="https://doi.org/10.1101/674440",
                 year="2019",
                 author="Nelson et al.",
                 reasons={stdpopsim.CiteReason.ENGINE}),
             ]}

    def simulate(
            self, demographic_model=None, contig=None, samples=None, seed=None,
            msprime_model=None, msprime_change_model=None, dry_run=False):
        """
        Simulate the demographic model using msprime.
        See :meth:`.Engine.simulate()` for definitions of parameters defined
        for all engines.

        :param msprime_model: The msprime simulation model to be used.
            One of ``hudson``, ``dtwf``, ``smc``, or ``smc_prime``.
            See msprime API documentation for details.
        :type msprime_model: str
        :param msprime_change_model: A list of (time, model) tuples, which
            changes the simulation model to the new model at the time specified.
        :type msprime_change_model: list of (float, str) tuples
        :param dry_run: If True, ``end_time=0`` is passed to :meth:`msprime.simulate()`
            to initialise the simulation and then immediately return.
        :type dry_run: bool
        """
        if msprime_model is None:
            msprime_model = self.supported_models[0]
        else:
            if msprime_model not in self.supported_models:
                raise ValueError(f"Unrecognised model '{msprime_model}'")
            if msprime_model in self.model_citations:
                self.citations.extend(self.model_citations[msprime_model])

        demographic_events = demographic_model.demographic_events.copy()
        if msprime_change_model is not None:
            for t, model in msprime_change_model:
                if model not in self.supported_models:
                    raise ValueError(f"Unrecognised model '{model}'")
                model_change = msprime.SimulationModelChange(t, model)
                demographic_events.append(model_change)
                if model in self.model_citations:
                    self.citations.extend(self.model_citations[model])
            demographic_events.sort(key=lambda x: x.time)

        ts = msprime.simulate(
                samples=samples,
                recombination_map=contig.recombination_map,
                mutation_rate=contig.mutation_rate,
                population_configurations=demographic_model.population_configurations,
                migration_matrix=demographic_model.migration_matrix,
                demographic_events=demographic_events,
                random_seed=seed,
                model=msprime_model,
                end_time=0 if dry_run else None)
        if dry_run:
            ts = None
        return ts

    def get_version(self):
        return msprime.__version__


register_engine(_MsprimeEngine())


def get_default_engine():
    """
    Returns the default simulation engine (msprime).
    """
    return get_engine("msprime")
