"""
Tests for the e. coli data definitions.
"""
import unittest

import stdpopsim
from tests import test_species


class TestSpecies(unittest.TestCase, test_species.SpeciesTestMixin):
    species = stdpopsim.get_species("EscCol")

    def test_basic_attributes(self):
        # From paper https://doi.org/10.1093/molbev/msw048
        # Ne taken from Table 2
        self.assertEqual(self.species.population_size, 1.8e8)
        # 20 minutes per generation
        generation_time = 1.0 / (525600 / 20)
        self.assertAlmostEqual(self.species.generation_time, generation_time)


class TestGenome(unittest.TestCase, test_species.GenomeTestMixin):
    """
    Tests for the e_coli genome.
    """
    genome = stdpopsim.get_species("EscCol").genome

    def test_basic_attributes(self):
        self.assertEqual(len(self.genome.chromosomes), 1)
