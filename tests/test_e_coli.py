"""
Tests for the e. coli data definitions.
"""
import unittest
import io

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
