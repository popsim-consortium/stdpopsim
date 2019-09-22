"""
Tests for the e. coli data definitions.
"""
import unittest
import io
import msprime
from stdpopsim import e_coli
from qc import e_coli_qc


class TestLapierreConstant(unittest.TestCase):
    """
    Basic tests for the LapierreConstant model.
    """
    def test_debug_runs(self):
        model = e_coli.LapierreConstant()
        output = io.StringIO()
        model.debug(output)
        s = output.getvalue()
        self.assertGreater(len(s), 0)

    def test_qc_model_equal(self):
        model = e_coli.LapierreConstant()
        self.assertTrue(model.equals(e_coli_qc.LapierreConstant()))


class TestGenericConstantSize(unittest.TestCase):
    """
    Basic tests for the GenericConstantSize model.
    """

    def test_simulation_runs(self):
        model = e_coli.GenericConstantSize()
        ts = msprime.simulate(
            samples=[msprime.Sample(0, 0), msprime.Sample(0, 0)],
            **model.asdict())
        self.assertEqual(ts.num_populations, 1)

    def test_debug_runs(self):
        model = e_coli.GenericConstantSize()
        output = io.StringIO()
        model.debug(output)
        s = output.getvalue()
        self.assertGreater(len(s), 0)


class TestGenericTwoEpoch(unittest.TestCase):
    """
    Basic tests for the GenericTwoEpoch model.
    """

    def test_simulation_runs(self):
        model = e_coli.GenericTwoEpoch()
        ts = msprime.simulate(
            samples=[msprime.Sample(0, 0), msprime.Sample(0, 0)],
            **model.asdict())
        self.assertEqual(ts.num_populations, 1)

    def test_debug_runs(self):
        model = e_coli.GenericTwoEpoch()
        output = io.StringIO()
        model.debug(output)
        s = output.getvalue()
        self.assertGreater(len(s), 0)

    def test_debug_runs_v2(self):
        model = e_coli.GenericTwoEpoch(100, 4)
        output = io.StringIO()
        model.debug(output)
        s = output.getvalue()
        self.assertGreater(len(s), 0)
