"""
Tests for the human data definitions.
"""
import unittest

import stdpopsim

from tests import test_species


class TestSpecies(unittest.TestCase, test_species.SpeciesTestMixin):
    species = stdpopsim.get_species("HomSap")

    def test_basic_attributes(self):
        self.assertEqual(self.species.population_size, 10 ** 4)
        self.assertEqual(self.species.generation_time, 30)


class TestGenome(unittest.TestCase, test_species.GenomeTestMixin):
    """
    Tests for the human genome.
    """

    genome = stdpopsim.get_species("HomSap").genome

    def test_basic_attributes(self):
        self.assertEqual(len(self.genome.chromosomes), 25)

    def test_recombination_rates(self):
        # recompute recombination rates from HapMapII_GRCh37 map then
        # compare the results to the current recombination rates for each chromosome
        genetic_map = "HapMapII_GRCh37"
        species = stdpopsim.get_species("HomSap")
        for chrom in self.genome.chromosomes:
            if chrom.id == "Y":
                with self.assertWarns(Warning):
                    contig = species.get_contig(chrom.id, genetic_map=genetic_map)
                print("HERE")
            else:
                contig = species.get_contig(chrom.id, genetic_map=genetic_map)
            self.assertAlmostEqual(
                chrom.recombination_rate,
                contig.recombination_map.mean_rate,
            )
