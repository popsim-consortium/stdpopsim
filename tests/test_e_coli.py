"""
Tests for the e. coli data definitions.
"""
import unittest

import stdpopsim
from tests import test_models
from tests import test_species
from qc import e_coli_qc


class TestGenome(unittest.TestCase, test_species.GenomeTestMixin):
    """
    Tests for the e_coli genome.
    """
    genome = stdpopsim.get_species("e_coli").genome

    def test_basic_attributes(self):
        self.assertEqual(len(self.genome.chromosomes), 1)


class TestLapierreConstant(unittest.TestCase, test_models.QcdModelTestMixin):
    """
    Basic tests for the LapierreConstant model.
    """
    model = stdpopsim.e_coli._LapierreConstant()
    qc_model = e_coli_qc.LapierreConstant()
