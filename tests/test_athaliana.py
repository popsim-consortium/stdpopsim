import unittest
import io

import stdpopsim
from tests import test_models
from qc import arabidopsis_thaliana_qc


class TestGenome(unittest.TestCase):
    """
    Tests for the arabidopsis_thaliana genome.
    """
    genome = stdpopsim.get_species("arabidopsis_thaliana").genome

    def test_basic_attributes(self):
        self.assertEqual(len(self.genome.chromosomes), 5)

    def test_str(self):
        s = str(self.genome)
        self.assertGreater(len(s), 0)
        self.assertIsInstance(s, str)

    def test_chromosome_lengths(self):
        genome = self.genome
        self.assertEqual(genome.get_chromosome("chr1").length, 30427671)
        self.assertEqual(genome.get_chromosome("chr2").length, 19698289)
        self.assertEqual(genome.get_chromosome("chr3").length, 23459830)
        self.assertEqual(genome.get_chromosome("chr4").length, 18585056)
        self.assertEqual(genome.get_chromosome("chr5").length, 26975502)


class TestDurvasula2017MSMC(unittest.TestCase, test_models.QcdModelTestMixin):
    """
    Basic tests for the Durvasula MSMC model.
    """
    model = stdpopsim.arabidopsis_thaliana._Durvasula2017MSMC()
    qc_model = arabidopsis_thaliana_qc.Durvasula2017MSMC()
