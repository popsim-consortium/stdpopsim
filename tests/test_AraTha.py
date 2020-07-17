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
        # 5 autosomes + Mt + Pt
        self.assertEqual(len(self.genome.chromosomes), 7)
