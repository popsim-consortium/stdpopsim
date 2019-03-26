"""
Tests for the human data definitions.
"""
import unittest
import io

import msprime
import numpy as np

from stdpopsim import homo_sapiens
from stdpopsim import genetic_maps

from qc import homo_sapiens_qc


class TestGenome(unittest.TestCase):
    """
    Tests for the human genome.
    """
    def test_basic_attributes(self):
        genome = homo_sapiens.genome
        self.assertEqual(genome.species, "homo_sapiens")
        self.assertEqual(genome.default_genetic_map, "HapmapII_GRCh37")
        self.assertEqual(len(genome.chromosomes), 24)

    def test_str(self):
        s = str(homo_sapiens.genome)
        self.assertGreater(len(s), 0)
        self.assertIsInstance(s, str)

    def test_chromosome_lengths(self):
        genome = homo_sapiens.genome
        self.assertEqual(genome.chromosomes["chr1"].length, 249250621)
        self.assertEqual(genome.chromosomes["chr2"].length, 243199373)
        self.assertEqual(genome.chromosomes["chr3"].length, 198022430)
        self.assertEqual(genome.chromosomes["chr4"].length, 191154276)
        self.assertEqual(genome.chromosomes["chr5"].length, 180915260)
        self.assertEqual(genome.chromosomes["chr6"].length, 171115067)
        self.assertEqual(genome.chromosomes["chr7"].length, 159138663)
        self.assertEqual(genome.chromosomes["chr8"].length, 146364022)
        self.assertEqual(genome.chromosomes["chr9"].length, 141213431)
        self.assertEqual(genome.chromosomes["chr10"].length, 135534747)
        self.assertEqual(genome.chromosomes["chr11"].length, 135006516)
        self.assertEqual(genome.chromosomes["chr12"].length, 133851895)
        self.assertEqual(genome.chromosomes["chr13"].length, 115169878)
        self.assertEqual(genome.chromosomes["chr14"].length, 107349540)
        self.assertEqual(genome.chromosomes["chr15"].length, 102531392)
        self.assertEqual(genome.chromosomes["chr16"].length, 90354753)
        self.assertEqual(genome.chromosomes["chr17"].length, 81195210)
        self.assertEqual(genome.chromosomes["chr18"].length, 78077248)
        self.assertEqual(genome.chromosomes["chr19"].length, 59128983)
        self.assertEqual(genome.chromosomes["chr20"].length, 63025520)
        self.assertEqual(genome.chromosomes["chr21"].length, 48129895)
        self.assertEqual(genome.chromosomes["chr22"].length, 51304566)
        self.assertEqual(genome.chromosomes["chrX"].length, 155270560)
        self.assertEqual(genome.chromosomes["chrY"].length, 59373566)

    def get_mean_rr_numpy(self, chromosome):
        chrom_length = chromosome.length
        recombination_map = chromosome.recombination_map()
        positions = np.array(recombination_map.get_positions())
        positions_diff = recombination_map.get_positions()[1:]
        positions_diff.append(chrom_length)
        positions_diff = np.array(positions_diff)
        window_sizes = positions_diff - positions
        weights = window_sizes / chrom_length
        rates = recombination_map.get_rates()
        return np.average(rates, weights=weights)

    def test_default_recombination_rates(self):
        # recompute default recombination rates from maps found at
        # "http://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/"
        # "20110106_recombination_hotspots/ then compare the results to
        # the current default recombination rates for each chromosome"

        # We catch the warning that will be thrown when we iterate over
        # Y chromosome.
        with self.assertWarns(Warning):
            for chrom in homo_sapiens.genome.chromosomes.values():
                default_rr = chrom.default_recombination_rate
                numpy_rr = self.get_mean_rr_numpy(chrom)
                # Assert that the difference in mean recombination rate
                # is small when computed with numpy.
                self.assertTrue(np.allclose(default_rr, numpy_rr))

    def test_default_recombination_rate(self):
        genome = homo_sapiens.genome
        highest_rr = 0
        lowest_rr = 10000
        for chrom in genome.chromosomes.values():
            rr = chrom.default_recombination_rate
            lowest_rr = min(lowest_rr, rr)
            highest_rr = max(highest_rr, rr)
        mean_genome_rr = genome.mean_recombination_rate
        self.assertGreater(mean_genome_rr, lowest_rr)
        self.assertGreater(highest_rr, mean_genome_rr)

    def test_warning_from_no_mapped_chromosome(self):
        """
        Test that a known chromosome throws a warning
        if there is no recombination map associated
        """
        with self.assertWarns(Warning):
            chrom = homo_sapiens.genome.chromosomes["chrY"]
            cm = chrom.recombination_map()
            self.assertIsInstance(cm, msprime.RecombinationMap)

    def test_warning_from_mapped_chromosome(self):
        # Test that a known chromosome can get a
        # recombination map associated with it
        chrom = homo_sapiens.genome.chromosomes["chr1"]
        cm = chrom.recombination_map()
        self.assertIsInstance(cm, msprime.RecombinationMap)

    def test_chromosome_errors(self):
        # Assert that unknown chromosomes throw a KeyError
        with self.assertRaises(KeyError):
            homo_sapiens.genome.chromosomes["jibberish"]


class TestGeneticMap(unittest.TestCase):
    """
    Basic tests for the GeneticMap class
    """

    def test_get_genetic_map(self):
        default_map = homo_sapiens.genome.default_genetic_map
        HapmapII_GRCh37 = genetic_maps.get_genetic_map("homo_sapiens", default_map)
        self.assertIsInstance(HapmapII_GRCh37, genetic_maps.GeneticMap)

    def test_unknown_get_chromosome_map(self):
        default_map = homo_sapiens.genome.default_genetic_map
        HapmapII_GRCh37 = genetic_maps.get_genetic_map("homo_sapiens", default_map)
        with self.assertRaises(ValueError):
            HapmapII_GRCh37.get_chromosome_map("jibberish")

    def test_known_get_chromosome_map(self):
        default_map = homo_sapiens.genome.default_genetic_map
        HapmapII_GRCh37 = genetic_maps.get_genetic_map("homo_sapiens", default_map)
        recombination_map = HapmapII_GRCh37.get_chromosome_map("chr1")
        self.assertIsInstance(recombination_map, msprime.RecombinationMap)


class TestGutenkunstOutOfAfrica(unittest.TestCase):
    """
    Basic tests for the GutenkunstThreePopOutOfAfrica model.
    """

    def test_simulation_runs(self):
        model = homo_sapiens.GutenkunstThreePopOutOfAfrica()
        ts = msprime.simulate(
            samples=[msprime.Sample(pop, 0) for pop in range(3)],
            **model.asdict())
        self.assertEqual(ts.num_populations, 3)

    def test_debug_runs(self):
        model = homo_sapiens.GutenkunstThreePopOutOfAfrica()
        output = io.StringIO()
        model.debug(output)
        s = output.getvalue()
        self.assertGreater(len(s), 0)


class TestTennessenOutOfAfrica(unittest.TestCase):
    """
    Basic tests for the TennessenTwoPopOutOfAfrica model.
    """

    def test_simulation_runs(self):
        model = homo_sapiens.TennessenTwoPopOutOfAfrica()
        ts = msprime.simulate(
            samples=[msprime.Sample(pop, 0) for pop in range(2)],
            **model.asdict())
        self.assertEqual(ts.num_populations, 2)

    def test_debug_runs(self):
        model = homo_sapiens.TennessenTwoPopOutOfAfrica()
        output = io.StringIO()
        model.debug(output)
        s = output.getvalue()
        self.assertGreater(len(s), 0)

    def test_qc_model_equal(self):
        model = homo_sapiens.TennessenTwoPopOutOfAfrica()
        self.assertTrue(model.equals(homo_sapiens_qc.TennessenTwoPopOutOfAfrica()))


class TestBrowningAmerica(unittest.TestCase):
    """
    Basic tests for the BrowningAmerica model.
    """

    def test_simulation_runs(self):
        model = homo_sapiens.BrowningAmerica()
        ts = msprime.simulate(
            samples=[msprime.Sample(pop, 0) for pop in range(4)],
            **model.asdict())
        self.assertEqual(ts.num_populations, 4)

    def test_debug_runs(self):
        model = homo_sapiens.BrowningAmerica()
        output = io.StringIO()
        model.debug(output)
        s = output.getvalue()
        self.assertGreater(len(s), 0)

    def test_qc_model_equal(self):
        model = homo_sapiens.BrowningAmerica()
        self.assertTrue(model.equals(homo_sapiens_qc.BrowningAmerica()))


class TestChromosomeFactory(unittest.TestCase):
    """
    Simple tests for the temporary chromosome factory object.

    TODO move these into the appropriate location later.
    """
    def test_length_multiplier(self):
        chrom1 = homo_sapiens.chromosome_factory("chr22")
        for x in [0.125, 1.0, 2.0]:
            chrom2 = homo_sapiens.chromosome_factory("chr22", length_multiplier=x)
            self.assertEqual(
                chrom1.recombination_map.get_positions()[-1] * x,
                chrom2.recombination_map.get_positions()[-1])

    def test_length_multiplier_on_empirical_map(self):
        with self.assertRaises(ValueError):
            homo_sapiens.chromosome_factory(
                "chr1", "HapmapII_GRCh37", length_multiplier=2)
