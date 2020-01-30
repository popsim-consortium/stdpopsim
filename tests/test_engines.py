"""
Tests for simulation engine infrastructure.
"""
import unittest
from unittest import mock

import stdpopsim
import stdpopsim.cli


class TestEngine_API(unittest.TestCase):
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
        self.assertRaises(NotImplementedError, e.simulate)
        self.assertRaises(NotImplementedError, e.get_version)


class TestMsprimeDTWF(unittest.TestCase):

    def test_engine_can_simulate(self):
        cmd = f"-e msprime-dtwf AraTha -l 0.001 --seed 1234 -q 10"
        with mock.patch("sys.stdout"):
            stdpopsim.cli.stdpopsim_main(cmd.split())

    def test_simulation_model_switching(self):
        cmd = f"-e msprime-dtwf --time-hudson 50 AraTha -l 0.001 --seed 1234 -q 10"
        with mock.patch("sys.stdout"):
            stdpopsim.cli.stdpopsim_main(cmd.split())
