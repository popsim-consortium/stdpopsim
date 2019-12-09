"""
Tests for the human data definitions.
"""
import unittest

import stdpopsim

from tests import test_models
from tests import test_species
from qc import homo_sapiens_qc


class TestSpecies(unittest.TestCase, test_species.SpeciesTestMixin):
    species = stdpopsim.get_species("homsap")

    def test_basic_attributes(self):
        self.assertEqual(self.species.population_size, 10**4)
        self.assertEqual(self.species.generation_time, 25)


class TestGenome(unittest.TestCase, test_species.GenomeTestMixin):
    """
    Tests for the human genome.
    """
    genome = stdpopsim.get_species("homsap").genome

    def test_basic_attributes(self):
        self.assertEqual(len(self.genome.chromosomes), 24)

    def test_chromosome_lengths(self):
        genome = self.genome
        self.assertEqual(genome.get_chromosome("chr1").length, 249250621)
        self.assertEqual(genome.get_chromosome("chr2").length, 243199373)
        self.assertEqual(genome.get_chromosome("chr3").length, 198022430)
        self.assertEqual(genome.get_chromosome("chr4").length, 191154276)
        self.assertEqual(genome.get_chromosome("chr5").length, 180915260)
        self.assertEqual(genome.get_chromosome("chr6").length, 171115067)
        self.assertEqual(genome.get_chromosome("chr7").length, 159138663)
        self.assertEqual(genome.get_chromosome("chr8").length, 146364022)
        self.assertEqual(genome.get_chromosome("chr9").length, 141213431)
        self.assertEqual(genome.get_chromosome("chr10").length, 135534747)
        self.assertEqual(genome.get_chromosome("chr11").length, 135006516)
        self.assertEqual(genome.get_chromosome("chr12").length, 133851895)
        self.assertEqual(genome.get_chromosome("chr13").length, 115169878)
        self.assertEqual(genome.get_chromosome("chr14").length, 107349540)
        self.assertEqual(genome.get_chromosome("chr15").length, 102531392)
        self.assertEqual(genome.get_chromosome("chr16").length, 90354753)
        self.assertEqual(genome.get_chromosome("chr17").length, 81195210)
        self.assertEqual(genome.get_chromosome("chr18").length, 78077248)
        self.assertEqual(genome.get_chromosome("chr19").length, 59128983)
        self.assertEqual(genome.get_chromosome("chr20").length, 63025520)
        self.assertEqual(genome.get_chromosome("chr21").length, 48129895)
        self.assertEqual(genome.get_chromosome("chr22").length, 51304566)
        self.assertEqual(genome.get_chromosome("chrX").length, 155270560)
        self.assertEqual(genome.get_chromosome("chrY").length, 59373566)

    def test_recombination_rates(self):
        # recompute recombination rates from HapmapII_GRCh37 map then
        # compare the results to the current recombination rates for each chromosome
        genetic_map = "HapmapII_GRCh37"
        species = stdpopsim.get_species("homsap")
        for chrom in self.genome.chromosomes:
            if chrom.id == "chrY":
                with self.assertWarns(Warning):
                    contig = species.get_contig(chrom.id, genetic_map=genetic_map)
            else:
                contig = species.get_contig(chrom.id, genetic_map=genetic_map)
            self.assertAlmostEqual(
                chrom.recombination_rate,
                contig.recombination_map.mean_recombination_rate)


species = stdpopsim.get_species("homsap")


class TestTennessenTwoPopOutOfAfrica(
        unittest.TestCase, test_models.QcdCatalogDemographicModelTestMixin):
    model = species.get_demographic_model("OutOfAfrica_2T12")
    qc_model = homo_sapiens_qc.TennessenTwoPopOutOfAfrica()


class TestTennessenOnePopAfrica(
        unittest.TestCase, test_models.QcdCatalogDemographicModelTestMixin):
    model = species.get_demographic_model("Africa_1T12")
    qc_model = homo_sapiens_qc.TennessenOnePopAfrica()


class TestBrowningAmerica(
        unittest.TestCase, test_models.QcdCatalogDemographicModelTestMixin):
    model = species.get_demographic_model("AmericanAdmixture_4B11")
    qc_model = homo_sapiens_qc.BrowningAmerica()


class TestRagsdaleArchaic(
        unittest.TestCase, test_models.QcdCatalogDemographicModelTestMixin):
    model = species.get_demographic_model("OutOfAfricaArchaicAdmixture_5R19")
    qc_model = homo_sapiens_qc.RagsdaleArchaic()


class TestKammAncientEurasia(
        unittest.TestCase, test_models.QcdCatalogDemographicModelTestMixin):
    model = species.get_demographic_model("AncientEurasia_9K19")
    qc_model = homo_sapiens_qc.KammAncientSamples()


# Models that have not been QC'd:

class TestGutenkunstThreePopOutOfAfrica(
        unittest.TestCase, test_models.CatalogDemographicModelTestMixin):
    model = species.get_demographic_model("OutOfAfrica_3G09")


class TestSchiffelsZigzag(
        unittest.TestCase, test_models.CatalogDemographicModelTestMixin):
    model = species.get_demographic_model("Zigzag_1S14")
