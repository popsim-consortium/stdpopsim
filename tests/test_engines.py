"""
Tests for simulation engine infrastructure.
"""
import unittest

import stdpopsim
import msprime


class TestEngineAPI(unittest.TestCase):
    """
    Tests for the API exposed for simulation engines.
    """

    def _test_engine(self, engine):
        self.assertIsInstance(engine, stdpopsim.Engine)
        self.assertIsNotNone(engine.simulate)
        self.assertNotEqual(len(engine.citations), 0)
        self.assertNotEqual(len(engine.get_version()), 0)

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
        self.assertEqual(engine1, engine2)
        # remove engine to avoid possible problems with other tests
        del stdpopsim.engines._registered_engines[engine1.id]

    def test_register_duplicate(self):
        engine = stdpopsim.get_default_engine()
        with self.assertRaises(ValueError):
            stdpopsim.register_engine(engine)

    def test_get_engine(self):
        with self.assertRaises(ValueError):
            stdpopsim.get_engine("nonexistent")

    def test_abstract_base_class(self):
        e = stdpopsim.Engine()
        with self.assertRaises(NotImplementedError):
            e.simulate(None, None, None)
        self.assertRaises(NotImplementedError, e.get_version)


class TestBehaviour(unittest.TestCase):
    def test_simulate_nonexistent_param(self):
        species = stdpopsim.get_species("HomSap")
        model = species.get_demographic_model("AshkSub_7G19")
        good_kwargs = dict(
            demographic_model=model,
            contig=species.get_contig("chr1"),
            samples=model.get_samples(10, 10, 10),
            dry_run=True,
        )
        bad_kwargs = good_kwargs.copy().update(nonexistent_param=None)
        for engine in stdpopsim.all_engines():
            engine.simulate(**good_kwargs)
            with self.assertRaises(TypeError):
                engine.simulate(**bad_kwargs)

    def test_required_params(self):
        species = stdpopsim.get_species("HomSap")
        model = species.get_demographic_model("AshkSub_7G19")
        contig = (species.get_contig("chr1"),)
        for engine in stdpopsim.all_engines():
            with self.assertRaises(TypeError):
                engine.simulate(model, contig)

    def test_msprime_kwargs(self):
        species = stdpopsim.get_species("HomSap")
        model = species.get_demographic_model("AshkSub_7G19")
        contig = species.get_contig("chr22", length_multiplier=0.01)
        samples = model.get_samples(10)
        engine = stdpopsim.get_engine("msprime")
        sim_arg = engine.simulate(
            model, contig, samples, record_full_arg=True, random_seed=1
        )
        assert any(msprime.NODE_IS_RE_EVENT == sim_arg.tables.nodes.flags)

    def test_msprime_seed(self):
        species = stdpopsim.get_species("HomSap")
        model = species.get_demographic_model("AshkSub_7G19")
        contig = species.get_contig("chr22", length_multiplier=0.01)
        samples = model.get_samples(10)
        engine = stdpopsim.get_engine("msprime")
        with self.assertRaises(ValueError):
            engine.simulate(model, contig, samples, seed=1, random_seed=1)
        sim_seed = engine.simulate(model, contig, samples, seed=1)
        sim_random_seed = engine.simulate(model, contig, samples, random_seed=1)
        self.assertEquals(sim_seed.tables.edges, sim_random_seed.tables.edges)
