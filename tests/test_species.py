"""
Tests for the genetic species interface.
"""
import unittest
import math

# import msprime

import stdpopsim
# import tests


class TestAllSpecies(unittest.TestCase):
    """
    Tests for basic methods on all species in the catalog.
    """
    def test_str(self):
        for species in stdpopsim.all_species():
            s = str(species)
            self.assertIsInstance(s, str)
            self.assertGreater(len(s), 0)


class GenomeTestMixin(object):
    """
    Mixin class for testing individual genome properties.
    """
    genome = None  # To be defined in subclasses.

    def test_str(self):
        s = str(self.genome)
        self.assertGreater(len(s), 0)
        self.assertIsInstance(s, str)

    def test_mean_recombination_rate(self):
        # test that the mean recombination rate lies between the max and min values
        # TODO this test case is failing in a couple of places. Is this really a
        # good test?
        highest_rr = 0
        lowest_rr = 1e200
        for chrom in self.genome.chromosomes:
            rr = chrom.recombination_rate
            lowest_rr = min(lowest_rr, rr)
            highest_rr = max(highest_rr, rr)
        mean_genome_rr = self.genome.mean_mutation_rate
        if not math.isclose(mean_genome_rr, lowest_rr):
            self.assertGreaterEqual(mean_genome_rr, lowest_rr)
            self.assertGreaterEqual(highest_rr, mean_genome_rr)

    def test_mean_mutation_rate(self):
        # test that the mean mutation rate lies between the max and min values
        highest_mr = 0
        lowest_mr = 1e200
        for chrom in self.genome.chromosomes:
            mr = chrom.mutation_rate
            lowest_mr = min(lowest_mr, mr)
            highest_mr = max(highest_mr, mr)
        mean_genome_mr = self.genome.mean_mutation_rate
        if not math.isclose(mean_genome_mr, lowest_mr):
            self.assertGreaterEqual(mean_genome_mr, lowest_mr)
            self.assertGreaterEqual(highest_mr, mean_genome_mr)


class TestAllGenomes(unittest.TestCase):
    """
    Tests for basic properties aon all genomes.
    """

#     def get_mean_rr_numpy(self, chromosome):
#         chrom_length = chromosome.length
#         recombination_map = chromosome.recombination_map()
#         positions = np.array(recombination_map.get_positions())
#         positions_diff = recombination_map.get_positions()[1:]
#         positions_diff.append(chrom_length)
#         positions_diff = np.array(positions_diff)
#         window_sizes = positions_diff - positions
#         weights = window_sizes / chrom_length
#         rates = recombination_map.get_rates()
#         return np.average(rates, weights=weights)

#     def test_default_recombination_rates(self):
#         # recompute default recombination rates from maps found at
#         # "http://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/"
#         # "20110106_recombination_hotspots/ then compare the results to
#         # the current default recombination rates for each chromosome"

#         # We catch the warning that will be thrown when we iterate over
#         # Y chromosome.
#         with self.assertWarns(Warning):
#             for chrom in self.genome.chromosomes:
#                 default_rr = chrom.default_recombination_rate
#                 numpy_rr = self.get_mean_rr_numpy(chrom)
#                 # Assert that the difference in mean recombination rate
#                 # is small when computed with numpy.
#                 self.assertTrue(np.allclose(default_rr, numpy_rr))

    # def test_warning_from_no_mapped_chromosome(self):
    #     """
    #     Test that a known chromosome throws a warning
    #     if there is no recombination map associated
    #     """
    #     with self.assertWarns(Warning):
    #         chrom = homo_sapiens.genome.chromosomes["chrY"]
    #         cm = chrom.recombination_map()
    #         self.assertIsInstance(cm, msprime.RecombinationMap)

    # def test_warning_from_mapped_chromosome(self):
    #     # Test that a known chromosome can get a
    #     # recombination map associated with it
    #     chrom = homo_sapiens.genome.chromosomes["chr1"]
    #     cm = chrom.recombination_map()
    #     self.assertIsInstance(cm, msprime.RecombinationMap)

    # def test_chromosome_errors(self):
    #     # Assert that unknown chromosomes throw a KeyError
    #     with self.assertRaises(KeyError):
    #         homo_sapiens.genome.chromosomes["jibberish"]


# @unittest.skip("Skip for now")
# class TestContigFactory(unittest.TestCase):
#     """
#     Tests for the contig_factory function.

#     TODO move these into the appropriate location.
#     """
#     def test_length_multiplier(self):
#         contig1 = genomes.contig_factory("homo_sapiens", "chr22")
#         for x in [0.125, 1.0, 2.0]:
#             contig2 = genomes.contig_factory(
#                 "homo_sapiens", "chr22", length_multiplier=x)
#             self.assertEqual(
#                 contig1.recombination_map.get_positions()[-1] * x,
#                 contig2.recombination_map.get_positions()[-1])

#     def test_length_multiplier_on_empirical_map(self):
#         with self.assertRaises(ValueError):
#             genomes.contig_factory(
#                 "homo_sapiens", "chr1", "HapmapII_GRCh37", length_multiplier=2)


# @unittest.skip("Skip for now")
# class TestGeneticMap(unittest.TestCase):
#     """
#     Basic tests for the GeneticMap class
#     """

#     def test_get_genetic_map(self):
#         default_map = homo_sapiens.genome.default_genetic_map
#         HapmapII_GRCh37 = genetic_maps.get_genetic_map("homo_sapiens", default_map)
#         self.assertIsInstance(HapmapII_GRCh37, genetic_maps.GeneticMap)

#     def test_unknown_get_chromosome_map(self):
#         default_map = homo_sapiens.genome.default_genetic_map
#         HapmapII_GRCh37 = genetic_maps.get_genetic_map("homo_sapiens", default_map)
#         with self.assertRaises(ValueError):
#             HapmapII_GRCh37.get_chromosome_map("jibberish")

#     def test_known_get_chromosome_map(self):
#         default_map = homo_sapiens.genome.default_genetic_map
#         HapmapII_GRCh37 = genetic_maps.get_genetic_map("homo_sapiens", default_map)
#         recombination_map = HapmapII_GRCh37.get_chromosome_map("chr1")
#         self.assertIsInstance(recombination_map, msprime.RecombinationMap)
