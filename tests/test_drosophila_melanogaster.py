"""
Tests for the drosophila_melanogaster data definitions.
"""
import unittest

import stdpopsim
from tests import test_models
from qc import drosophlia_melanogaster_qc


# TODO refactor most of this into a superclass.

@unittest.skip("Skip for now")
class TestGenome(unittest.TestCase):
    """
    Tests for the drosophila_melanogaster genome.
    """
    # def test_basic_attributes(self):
    #     genome = drosophila_melanogaster.genome
    #     self.assertEqual(genome.species, "drosophila_melanogaster")
    #     self.assertEqual(genome.default_genetic_map, "Comeron2012_dm6")
    #     self.assertEqual(len(genome.chromosomes), 8)

    # def test_str(self):
    #     s = str(drosophila_melanogaster.genome)
    #     self.assertGreater(len(s), 0)
    #     self.assertIsInstance(s, str)

    # def test_chromosome_lengths(self):
    #     genome = drosophila_melanogaster.genome
    #     # Numbers from DM6 release
    #     # `dm6 <https://www.ncbi.nlm.nih.gov/assembly/GCF_000001215.4/>`_.
    #     self.assertEqual(genome.chromosomes["chr2L"].length, 23513712)
    #     self.assertEqual(genome.chromosomes["chr2R"].length, 25286936)
    #     self.assertEqual(genome.chromosomes["chr3L"].length, 28110227)
    #     self.assertEqual(genome.chromosomes["chr3R"].length, 32079331)
    #     self.assertEqual(genome.chromosomes["chrX"].length, 23542271)
    #     self.assertEqual(genome.chromosomes["chr4"].length, 1348131)
    #     self.assertEqual(genome.chromosomes["chrY"].length, 3667352)


class TestSheehanSongThreeEpoch(unittest.TestCase, test_models.QcdModelTestMixin):
    model = stdpopsim.drosophila_melanogaster._SheehanSongThreeEpoch()
    qc_model = drosophlia_melanogaster_qc.SheehanSongThreeEpic()


class TestLiStephanTwoPopulation(unittest.TestCase, test_models.QcdModelTestMixin):
    model = stdpopsim.drosophila_melanogaster._LiStephanTwoPopulation()
    qc_model = drosophlia_melanogaster_qc.LiStephanTwoPopulation()
