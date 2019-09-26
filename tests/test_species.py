"""
Tests for the genetic species interface.
"""
# import unittest

# import msprime

# import stdpopsim
# from stdpopsim import genetic_maps
# import tests

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
