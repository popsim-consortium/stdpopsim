"""
Tests for simulation engine infrastructure.
"""

import json

import msprime
import numpy as np
import pytest

import stdpopsim


class TestEngineAPI:
    """
    Tests for the API exposed for simulation engines.
    """

    def _test_engine(self, engine):
        assert isinstance(engine, stdpopsim.Engine)
        assert engine.simulate is not None
        assert len(engine.citations) != 0
        assert len(engine.get_version()) != 0

    def test_get_default_engine(self):
        engine = stdpopsim.get_default_engine()
        self._test_engine(engine)

    def test_all_engines(self):
        for engine in stdpopsim.all_engines():
            self._test_engine(engine)

    def test_register_engine(self):
        class MyEngine(stdpopsim.Engine):
            id = "test-engine"
            name = "test"
            citations = []

        engine1 = MyEngine()
        stdpopsim.register_engine(engine1)
        engine2 = stdpopsim.get_engine(engine1.id)
        assert engine1 == engine2
        # remove engine to avoid possible problems with other tests
        del stdpopsim.engines._registered_engines[engine1.id]

    def test_register_duplicate(self):
        engine = stdpopsim.get_default_engine()
        with pytest.raises(ValueError):
            stdpopsim.register_engine(engine)

    def test_get_engine(self):
        with pytest.raises(ValueError):
            stdpopsim.get_engine("nonexistent")

    def test_abstract_base_class(self):
        e = stdpopsim.Engine()
        with pytest.raises(NotImplementedError):
            e.simulate(None, None, None)
        with pytest.raises(NotImplementedError):
            e.get_version()


@pytest.mark.filterwarnings(
    "ignore:.*model has mutation rate.*but this simulation used.*"
)
class TestBehaviour:
    def test_simulate_nonexistent_param(self):
        species = stdpopsim.get_species("HomSap")
        model = species.get_demographic_model("AshkSub_7G19")
        good_kwargs = dict(
            demographic_model=model,
            contig=species.get_contig("chr1"),
            samples={"YRI": 5, "CHB": 5, "CEU": 5},
            dry_run=True,
        )
        bad_kwargs = good_kwargs.copy().update(nonexistent_param=None)
        for engine in stdpopsim.all_engines():
            engine.simulate(**good_kwargs)
            with pytest.raises(TypeError):
                engine.simulate(**bad_kwargs)

    def test_required_params(self):
        species = stdpopsim.get_species("HomSap")
        model = species.get_demographic_model("AshkSub_7G19")
        contig = (species.get_contig("chr1"),)
        for engine in stdpopsim.all_engines():
            with pytest.raises(TypeError):
                engine.simulate(model, contig)

    def test_msprime_kwargs(self):
        species = stdpopsim.get_species("HomSap")
        model = species.get_demographic_model("AshkSub_7G19")
        contig = species.get_contig("chr22", right=5e5)
        samples = {"YRI": 5}
        engine = stdpopsim.get_engine("msprime")
        sim_arg = engine.simulate(
            model, contig, samples, record_full_arg=True, random_seed=1
        )
        assert any(msprime.NODE_IS_RE_EVENT == sim_arg.tables.nodes.flags)

    def test_msprime_seed(self):
        species = stdpopsim.get_species("HomSap")
        model = species.get_demographic_model("AshkSub_7G19")
        contig = species.get_contig("chr22", right=5e5)
        samples = {"YRI": 5}
        engine = stdpopsim.get_engine("msprime")
        with pytest.raises(ValueError):
            engine.simulate(model, contig, samples, seed=1, random_seed=1)
        sim_seed = engine.simulate(model, contig, samples, seed=1)
        sim_random_seed = engine.simulate(model, contig, samples, random_seed=1)
        assert sim_seed.tables.edges == sim_random_seed.tables.edges

    def test_non_neutral_contig(self):
        species = stdpopsim.get_species("HomSap")
        model = species.get_demographic_model("AshkSub_7G19")
        samples = {"YRI": 5}
        contig = stdpopsim.Contig.basic_contig(length=100)
        contig.clear_dfes()
        props = [1]
        mt = [stdpopsim.MutationType(distribution_type="f", distribution_args=[1])]
        dfes = [
            stdpopsim.DFE(
                id=str(0),
                description="test",
                long_description="test test",
                proportions=props,
                mutation_types=mt,
            )
        ]
        contig.add_dme(intervals=np.array([[0, 50]]), DME=dfes[0])
        engine = stdpopsim.get_engine("msprime")
        with pytest.raises(ValueError):
            engine.simulate(model, contig, samples, seed=1)
        # okay now change selection coefficient to neutral
        contig.clear_dfes()
        props = [1]
        mt = [stdpopsim.MutationType(distribution_type="f", distribution_args=[0])]
        dfes = [
            stdpopsim.DFE(
                id=str(0),
                description="test",
                long_description="test test",
                proportions=props,
                mutation_types=mt,
            )
        ]
        contig.add_dme(intervals=np.array([[0, 50]]), DME=dfes[0])
        engine = stdpopsim.get_engine("msprime")
        engine.simulate(model, contig, samples, seed=1)

    def test_gene_conversion(self):
        species = stdpopsim.get_species("DroMel")
        model = species.get_demographic_model("African3Epoch_1S16")
        samples = {"AFR": 5}
        contig = species.get_contig(length=1000, use_species_gene_conversion=True)
        assert contig.gene_conversion_fraction > 0
        assert contig.gene_conversion_length > 0
        engine = stdpopsim.get_engine("msprime")
        engine.simulate(model, contig, samples, seed=1)

    def test_msprime_bad_samples(self):
        engine = stdpopsim.get_engine("msprime")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr1")
        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        samples = [1, 2, ["foo"]]
        with pytest.raises(ValueError, match="Samples must be a dict"):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
            )

    def check_execution_order(self, expected, msprime_model, msprime_change_model=None):
        """
        Run a short simulation and verify that the model execution order
        recorded in the tree sequence provenance matches `expected`.
        Each entry in `expected` is a (ClassName, duration) tuple.
        """
        species = stdpopsim.get_species("AraTha")
        engine = stdpopsim.get_engine("msprime")
        kwargs = dict(
            demographic_model=stdpopsim.PiecewiseConstantSize(species.population_size),
            contig=species.get_contig(length=100),
            samples={"pop_0": 5},
            msprime_model=msprime_model,
        )
        if msprime_change_model is not None:
            kwargs["msprime_change_model"] = msprime_change_model
        ts = engine.simulate(**kwargs)
        prov = ts.provenance(0)
        record = json.loads(prov.record)
        models = record["parameters"]["model"]
        # Flexible provenance parsing
        if isinstance(models, dict):
            models = [models]
        actual = []
        for model in models:
            name = model["__class__"].replace("msprime.ancestry.", "")
            duration = model["duration"]
            actual.append((name, duration))
        assert expected == actual

    def test_msprime_ancestry_model(self):
        engine = stdpopsim.get_engine("msprime")
        species = stdpopsim.get_species("HomSap")
        kwargs = dict(
            demographic_model=species.get_demographic_model("AshkSub_7G19"),
            contig=species.get_contig("chr1"),
            samples={"YRI": 5, "CHB": 5, "CEU": 5},
            dry_run=True,
        )
        # Valid single None, string or `msprime.AncestryModel`
        self.check_execution_order([("StandardCoalescent", None)], None)
        self.check_execution_order([("StandardCoalescent", None)], "hudson")
        self.check_execution_order(
            [("StandardCoalescent", None)], msprime.StandardCoalescent()
        )
        self.check_execution_order([("SMCK", None)], msprime.SMCK(k=1))

        # Fail if provided an invalid `msprime.AncestryModel` subclass
        class MyClass:
            pass

        with pytest.raises(TypeError):
            engine.simulate(**kwargs, msprime_model=MyClass())

        # Fail if msprime_model is a list — passing a list is not supported;
        # use msprime_change_model instead.
        with pytest.raises(TypeError):
            engine.simulate(
                **kwargs,
                # Note: this is valid input for `msprime`
                msprime_model=[
                    msprime.DiscreteTimeWrightFisher(duration=10),
                    msprime.SMCK(k=1),
                ],
            )

        # Test that we do not mutate the input model
        model = msprime.SMCK(k=1, duration=10)
        engine.simulate(
            **kwargs,
            msprime_model=model,
            msprime_change_model=[(12, model), (99, model)],
        )
        assert model.duration == 10
        model = msprime.SMCK(k=1, duration=10)
        engine.simulate(**kwargs, msprime_model=model)
        assert model.duration == 10, "model should not be mutated"

    def test_msprime_change_model(self):
        expected = [
            ("DiscreteTimeWrightFisher", 20),
            ("StandardCoalescent", None),
        ]
        self.check_execution_order(expected, "dtwf", [(20, "hudson")])
        self.check_execution_order(
            expected,
            msprime.DiscreteTimeWrightFisher(),
            [(20, msprime.StandardCoalescent())],
        )

        model_changes = [
            (20, msprime.SMCK(k=1)),
            (21, msprime.SMCK(k=0)),
            (100, "hudson"),
        ]
        expected = [
            ("DiscreteTimeWrightFisher", 20),
            ("SMCK", 21 - 20),
            ("SMCK", 100 - 21),
            ("StandardCoalescent", None),
        ]
        self.check_execution_order(expected, "dtwf", model_changes)

    def test_msprime_change_model_errors(self):
        engine = stdpopsim.get_engine("msprime")
        species = stdpopsim.get_species("HomSap")
        kwargs = dict(
            demographic_model=species.get_demographic_model("AshkSub_7G19"),
            contig=species.get_contig("chr1"),
            samples={"YRI": 5, "CHB": 5, "CEU": 5},
        )
        with pytest.raises(
            ValueError,
            match=r"`model_changes` must be a list of \(time, model\) tuples",
        ):
            engine.simulate(
                **kwargs,
                msprime_model="dtwf",
                msprime_change_model=[
                    ("hudson"),
                ],
            )

        invalid_models = ["notvalid", None, 100]
        for invalid in invalid_models:
            with pytest.raises((TypeError, ValueError)):
                engine.simulate(
                    **kwargs,
                    msprime_model="dtwf",
                    msprime_change_model=[
                        (10, invalid),
                    ],
                )

        invalid_times = [-1, None, "hudson", np.inf]
        for invalid in invalid_times:
            with pytest.raises(
                ValueError,
                match=r"Specified time .* is not a valid non-negative number",
            ):
                engine.simulate(
                    **kwargs,
                    msprime_model="dtwf",
                    msprime_change_model=[
                        (invalid, "hudson"),
                    ],
                )
        with pytest.raises(
            ValueError,
            match=r"Tuples in `model_changes` must be sorted",
        ):
            engine.simulate(
                **kwargs,
                msprime_model="dtwf",
                msprime_change_model=[
                    (2, "hudson"),
                    (1, "dtwf"),
                ],
            )
