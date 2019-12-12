"""
Tests for the genetic species interface.
"""
import unittest
import math

import msprime

import stdpopsim


class TestSpecies(unittest.TestCase):
    """
    Tests for basic methods for the species.
    """
    def test_str(self):
        for species in stdpopsim.all_species():
            s = str(species)
            self.assertIsInstance(s, str)
            self.assertGreater(len(s), 0)

    def test_get_known_species(self):
        good = ["HomSap", "EscCol"]
        for species_id in good:
            species = stdpopsim.get_species(species_id)
            self.assertIsInstance(species, stdpopsim.Species)
            self.assertEqual(species.id, species_id)

    def test_get_unknown_species(self):
        bad = ["XXXX", ""]
        for species_name in bad:
            with self.assertRaises(ValueError):
                stdpopsim.get_species(species_name)

    def test_add_duplicate_species(self):
        species = stdpopsim.get_species("HomSap")
        with self.assertRaises(ValueError):
            stdpopsim.register_species(species)

    def test_get_known_genetic_map(self):
        good = ["HapmapII_GRCh37", "Decode_2010_sex_averaged"]
        species = stdpopsim.get_species("HomSap")
        for name in good:
            gmap = species.get_genetic_map(name)
            self.assertIsInstance(gmap, stdpopsim.GeneticMap)
            self.assertEqual(gmap.name, name)

    def test_get_unknown_genetic_map(self):
        bad = ["GDXXX", "", None]
        species = stdpopsim.get_species("HomSap")
        for name in bad:
            with self.assertRaises(ValueError):
                species.get_genetic_map(name)

    def test_add_duplicate_genetic_map(self):
        species = stdpopsim.get_species("HomSap")
        genetic_map = species.get_genetic_map("HapmapII_GRCh37")
        with self.assertRaises(ValueError):
            species.add_genetic_map(genetic_map)

    def test_add_duplicate_model(self):
        species = stdpopsim.get_species("HomSap")
        model = species.get_demographic_model("OutOfAfrica_3G09")
        with self.assertRaises(ValueError):
            species.add_demographic_model(model)


class SpeciesTestMixin(object):
    """
    Mixin class for testing individual species properties.
    """
    species = None  # To be defined in subclasses.

    def test_str(self):
        s = str(self.species)
        self.assertGreater(len(s), 0)
        self.assertIsInstance(s, str)


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
        mean_genome_rr = self.genome.mean_recombination_rate
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
    def test_str(self):
        for species in stdpopsim.all_species():
            s = str(species.genome)
            self.assertIsInstance(s, str)
            self.assertGreater(len(s), 0)


class TestGetContig(unittest.TestCase):
    """
    Tests for the get contig method.
    """
    species = stdpopsim.get_species("HomSap")

    def test_length_multiplier(self):
        contig1 = self.species.get_contig("chr22")
        for x in [0.125, 1.0, 2.0]:
            contig2 = self.species.get_contig("chr22", length_multiplier=x)
            self.assertEqual(
                contig1.recombination_map.get_positions()[-1] * x,
                contig2.recombination_map.get_positions()[-1])

    def test_length_multiplier_on_empirical_map(self):
        with self.assertRaises(ValueError):
            self.species.get_contig(
                "chr1", genetic_map="HapmapII_GRCh37", length_multiplier=2)

    def test_genetic_map(self):
        # TODO we should use a different map here so we're not hitting the cache.
        contig = self.species.get_contig("chr22", genetic_map="HapmapII_GRCh37")
        self.assertIsInstance(contig.recombination_map, msprime.RecombinationMap)
