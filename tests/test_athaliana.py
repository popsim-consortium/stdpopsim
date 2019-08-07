import unittest
import io

from stdpopsim import arabidopsis_thaliana
from qc import arabidopsis_thaliana_qc


class TestGenome(unittest.TestCase):
    """
    Tests for the arabidopsis_thaliana genome.
    """

    def test_basic_attributes(self):
        genome = arabidopsis_thaliana.genome
        self.assertEqual(genome.species, "arabidopsis_thaliana")
        self.assertEqual(genome.default_genetic_map, "Salome2012")
        self.assertEqual(len(genome.chromosomes), 5)

    def test_str(self):
        s = str(arabidopsis_thaliana.genome)
        self.assertGreater(len(s), 0)
        self.assertIsInstance(s, str)

    def test_chromosome_lengths(self):
        genome = arabidopsis_thaliana.genome
        self.assertEqual(genome.chromosomes["chr1"].length, 30427671)
        self.assertEqual(genome.chromosomes["chr2"].length, 19698289)
        self.assertEqual(genome.chromosomes["chr3"].length, 23459830)
        self.assertEqual(genome.chromosomes["chr4"].length, 18585056)
        self.assertEqual(genome.chromosomes["chr5"].length, 26975502)


class TestDurvasula2017MSMC(unittest.TestCase):
    """
    Basic tests for the Durvasula MSMC model.
    """

    def test_debug_runs(self):
        model = arabidopsis_thaliana.Durvasula2017MSMC()
        output = io.StringIO()
        model.debug(output)
        s = output.getvalue()
        self.assertGreater(len(s), 0)

    def test_qc_model_equal(self):
        model = arabidopsis_thaliana.Durvasula2017MSMC()
        self.assertTrue(model.equals(arabidopsis_thaliana_qc.Durvasula2017MSMC()))
