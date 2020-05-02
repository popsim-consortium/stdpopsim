import unittest

import stdpopsim
from tests import test_species


class TestSpecies(unittest.TestCase, test_species.SpeciesTestMixin):
    species = stdpopsim.get_species("AraTha")

    def test_basic_attributes(self):
        self.assertEqual(self.species.population_size, 10**4)
        self.assertEqual(self.species.generation_time, 1)


class TestGenome(unittest.TestCase, test_species.GenomeTestMixin):
    """
    Tests for the arabidopsis_thaliana genome.
    """
    genome = stdpopsim.get_species("AraTha").genome

    def test_basic_attributes(self):
        self.assertEqual(len(self.genome.chromosomes), 5)

    def test_chromosome_lengths(self):
        genome = self.genome
        self.assertEqual(genome.get_chromosome("chr1").length, 30427671)
        self.assertEqual(genome.get_chromosome("chr2").length, 19698289)
        self.assertEqual(genome.get_chromosome("chr3").length, 23459830)
        self.assertEqual(genome.get_chromosome("chr4").length, 18585056)
        self.assertEqual(genome.get_chromosome("chr5").length, 26975502)
