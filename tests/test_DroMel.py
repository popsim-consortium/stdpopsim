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
        self.assertEqual(len(self.genome.chromosomes), 7)

    def test_chromosome_lengths(self):
        genome = self.genome
        # Numbers from DM6 release
        # `dm6 <https://www.ncbi.nlm.nih.gov/assembly/GCF_000001215.4/>`_.
        self.assertEqual(genome.get_chromosome("chr2L").length, 23513712)
        self.assertEqual(genome.get_chromosome("chr2R").length, 25286936)
        self.assertEqual(genome.get_chromosome("chr3L").length, 28110227)
        self.assertEqual(genome.get_chromosome("chr3R").length, 32079331)
        self.assertEqual(genome.get_chromosome("chrX").length, 23542271)
        self.assertEqual(genome.get_chromosome("chr4").length, 1348131)
        self.assertEqual(genome.get_chromosome("chrY").length, 3667352)
