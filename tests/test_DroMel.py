"""
Tests for the drosophila_melanogaster data definitions.
"""
import unittest

import stdpopsim
from tests import test_species


class TestSpecies(unittest.TestCase, test_species.SpeciesTestMixin):
    species = stdpopsim.get_species("DroMel")

    def test_basic_attributes(self):
        self.assertEqual(self.species.population_size, 1720600)
        self.assertEqual(self.species.generation_time, 0.1)


class TestGenome(unittest.TestCase, test_species.GenomeTestMixin):
    """
    Tests for the drosophila_melanogaster genome.
    """
    genome = stdpopsim.get_species("DroMel").genome

    def test_basic_attributes(self):
        self.assertEqual(len(self.genome.chromosomes), 8)

    def test_non_recombining_chrs(self):
        chrom = self.genome.get_chromosome("Y")
        self.assertEqual(chrom.recombination_rate, 0)
        chrom = self.genome.get_chromosome("mitochondrion_genome")
        self.assertEqual(chrom.recombination_rate, 0)

    def test_recombining_chrs(self):
        for name in ["2L", "2R", "3L", "3R", "4", "X"]:
            chrom = self.genome.get_chromosome(name)
            self.assertGreater(chrom.recombination_rate, 0)
