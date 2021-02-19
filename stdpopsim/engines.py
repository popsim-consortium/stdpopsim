import logging
import warnings

import attr
import msprime
import stdpopsim
import numpy as np

logger = logging.getLogger(__name__)

_registered_engines = {}


def register_engine(engine):
    """
    Registers the specified simulation engine.

    :param engine: The simulation engine object to register.
    :type engine: :class:`.Engine`
    """
    if engine.id in _registered_engines:
        raise ValueError(f"Simulation engine '{engine.id}' already registered.")
    logger.debug(f"Registering simulation engine '{engine.id}'")
    _registered_engines[engine.id] = engine


def get_engine(id):
    """
    Returns the simulation engine with the specified ``id``.

    :param str id: The string identifier for the requested engine.
        The currently supported engines are "msprime" and "slim".
    :return: A simulation engine object with a ``simulate()`` method.
    :rtype: :class:`.Engine`
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
class Engine:
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
        self,
        demographic_model,
        contig,
        samples,
        *,
        seed=None,
        dry_run=False,
    ):
        """
        Simulates the model for the specified contig and samples. ``demographic_model``,
        ``contig``, and ``samples`` must be specified.

        :param demographic_model: The demographic model to simulate.
        :type demographic_model: :class:`.DemographicModel`
        :param contig: The contig, defining the length, mutation rate,
            and recombination rate(s).
        :type contig: :class:`.Contig`
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

    def _warn_zigzag(self, demographic_model):
        if demographic_model.id == "Zigzag_1S14":
            warnings.warn(
                "In stdpopsim <= 0.1.2, the Zigzag_1S14 model produced population "
                "sizes 5x lower than those in Schiffels & Durbin (2014). "
                "The population sizes have now been corrected. For details see: "
                "https://github.com/popsim-consortium/stdpopsim/issues/745"
            )


class _MsprimeEngine(Engine):
    id = "msprime"  #:
    description = "Msprime coalescent simulator"  #:
    citations = [
        stdpopsim.Citation(
            doi="https://doi.org/10.1371/journal.pcbi.1004842",
            year="2016",
            author="Kelleher et al.",
            reasons={stdpopsim.CiteReason.ENGINE},
        )
    ]
    # We default to the first model in the list.
    supported_models = ["hudson", "dtwf", "smc", "smc_prime"]
    model_citations = {
        "dtwf": [
            stdpopsim.Citation(
                doi="https://doi.org/10.1371/journal.pgen.1008619",
                year="2020",
                author="Nelson et al.",
                reasons={stdpopsim.CiteReason.ENGINE},
            )
        ]
    }

    def simulate(
        self,
        demographic_model,
        contig,
        samples,
        *,
        seed=None,
        msprime_model=None,
        msprime_change_model=None,
        dry_run=False,
        **kwargs,
    ):
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
        :param \\**kwargs: Further arguments passed to :meth:`msprime.simulate()`
        """
        if msprime_model is None:
            msprime_model = self.supported_models[0]
        else:
            if msprime_model not in self.supported_models:
                raise ValueError(f"Unrecognised model '{msprime_model}'")
            if msprime_model in self.model_citations:
                self.citations.extend(self.model_citations[msprime_model])

        if msprime_change_model is not None:
            msprime_model = [msprime_model]
            for t, model in msprime_change_model:
                if model not in self.supported_models:
                    raise ValueError(f"Unrecognised model '{model}'")
                msprime_model.append((t, model))
                if model in self.model_citations:
                    self.citations.extend(self.model_citations[model])

        if "random_seed" in kwargs.keys():
            if seed is None:
                seed = kwargs["random_seed"]
                del kwargs["random_seed"]
            else:
                raise ValueError("Cannot set both seed and random_seed")

        # TODO: remove this after a release or two. See #745.
        self._warn_zigzag(demographic_model)

        rng = np.random.default_rng(seed)
        seeds = rng.integers(1, 2 ** 31 - 1, size=2)

        ts = msprime.sim_ancestry(
            samples=samples,
            recombination_rate=contig.recombination_map,
            demography=demographic_model.model,
            ploidy=2,
            random_seed=seeds[0],
            model=msprime_model,
            end_time=0 if dry_run else None,
            **kwargs,
        )
        ts = msprime.sim_mutations(
            ts,
            end_time=0 if dry_run else None,
            random_seed=seeds[1],
            rate=contig.mutation_rate,
        )

        if contig.inclusion_mask is not None:
            ts = stdpopsim.utils.mask_tree_sequence(ts, contig.inclusion_mask, False)
        if contig.exclusion_mask is not None:
            ts = stdpopsim.utils.mask_tree_sequence(ts, contig.exclusion_mask, True)

        if dry_run:
            ts = None
        return ts

    def get_version(self):
        return msprime.__version__


register_engine(_MsprimeEngine())


def get_default_engine():
    """
    Returns the default simulation engine (msprime).

    :rtype: :class:`.Engine`
    """
    return get_engine("msprime")
