"""
Tests for the e. coli data definitions.
"""
import unittest
import io

import stdpopsim
from tests import test_models
from qc import e_coli_qc


class TestLapierreConstant(unittest.TestCase, test_models.QcdModelTestMixin):
    """
    Basic tests for the LapierreConstant model.
    """
    model = stdpopsim.e_coli._LapierreConstant()
    qc_model = e_coli_qc.LapierreConstant()
