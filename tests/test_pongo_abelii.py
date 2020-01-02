import unittest
import stdpopsim
from stdpopsim import pongo_abelii
from tests import test_models
from tests import test_species
from qc import pongo_abelii_qc


class TestSpecies(unittest.TestCase, test_species.SpeciesTestMixin):
    species = stdpopsim.get_species("PonAbe")

    def test_basic_attributes(self):
        self.assertEqual(self.species.population_size, 1.79*10**4)
        self.assertEqual(self.species.generation_time, 20)


class TestGenome(unittest.TestCase, test_species.GenomeTestMixin):
    """
    Tests for the Pongo abelii genome.
    """
    genome = stdpopsim.get_species("PonAbe").genome

    def test_basic_attributes(self):
        self.assertEqual(len(self.genome.chromosomes), 24)

    def test_chromosome_lengths(self):
        genome = self.genome
        self.assertEqual(genome.get_chromosome("chr1").length, 229942017)
        self.assertEqual(genome.get_chromosome("chr2a").length, 113028656)
        self.assertEqual(genome.get_chromosome("chr2b").length, 135000294)
        self.assertEqual(genome.get_chromosome("chr3").length, 202140232)
        self.assertEqual(genome.get_chromosome("chr4").length, 198332218)
        self.assertEqual(genome.get_chromosome("chr5").length, 183952662)
        self.assertEqual(genome.get_chromosome("chr6").length, 174210431)
        self.assertEqual(genome.get_chromosome("chr7").length, 157549271)
        self.assertEqual(genome.get_chromosome("chr8").length, 153482349)
        self.assertEqual(genome.get_chromosome("chr9").length, 135191526)
        self.assertEqual(genome.get_chromosome("chr10").length, 133410057)
        self.assertEqual(genome.get_chromosome("chr11").length, 132107971)
        self.assertEqual(genome.get_chromosome("chr12").length, 136387465)
        self.assertEqual(genome.get_chromosome("chr13").length, 117095149)
        self.assertEqual(genome.get_chromosome("chr14").length, 108868599)
        self.assertEqual(genome.get_chromosome("chr15").length, 99152023)
        self.assertEqual(genome.get_chromosome("chr16").length, 77800216)
        self.assertEqual(genome.get_chromosome("chr17").length, 73212453)
        self.assertEqual(genome.get_chromosome("chr18").length, 94050890)
        self.assertEqual(genome.get_chromosome("chr19").length, 60714840)
        self.assertEqual(genome.get_chromosome("chr20").length, 62736349)
        self.assertEqual(genome.get_chromosome("chr21").length, 48394510)
        self.assertEqual(genome.get_chromosome("chr22").length, 46535552)
        self.assertEqual(genome.get_chromosome("chrX").length, 156195299)


class TestPongo(unittest.TestCase, test_models.QcdCatalogDemographicModelTestMixin):
    model = pongo_abelii._orangutan()
    qc_model = pongo_abelii_qc.LockePongo()
