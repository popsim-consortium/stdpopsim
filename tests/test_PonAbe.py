import unittest

import stdpopsim
from tests import test_species


class TestSpecies(unittest.TestCase, test_species.SpeciesTestMixin):
    species = stdpopsim.get_species("PonAbe")

    def test_basic_attributes(self):
        self.assertEqual(self.species.population_size, 1.79 * 10 ** 4)
        self.assertEqual(self.species.generation_time, 20)


class TestGenome(unittest.TestCase, test_species.GenomeTestMixin):
    """
    Tests for the Pongo abelii genome.
    """

    genome = stdpopsim.get_species("PonAbe").genome

    def test_basic_attributes(self):
        # 23 + X + MT
        self.assertEqual(len(self.genome.chromosomes), 25)
