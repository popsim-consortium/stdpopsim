import unittest


import stdpopsim
from stdpopsim import pongo
from tests import test_models
from tests import test_species
from qc import pongo_qc


class TestSpecies(unittest.TestCase, test_species.SpeciesTestMixin):
    species = stdpopsim.get_species("pongo")

    def test_basic_attributes(self):
        self.assertEqual(self.species.population_size, 1.79*10**4)
        self.assertEqual(self.species.generation_time, 20)


class TestGenome(unittest.TestCase, test_species.GenomeTestMixin):
    """
    Tests for the pongo genome.
    """
    genome = stdpopsim.get_species("pongo").genome

    def test_basic_attributes(self):
        self.assertEqual(len(self.genome.chromosomes), 24)

    def test_chromosome_lengths(self):
        genome = self.genome
        self.assertEqual(genome.get_chromosome("chr1").length, 227913704)
        self.assertEqual(genome.get_chromosome("chr2A").length, 109511694)
        self.assertEqual(genome.get_chromosome("chr2B").length, 129937803)
        self.assertEqual(genome.get_chromosome("chr3").length, 193656255)
        self.assertEqual(genome.get_chromosome("chr4").length, 189387572)
        self.assertEqual(genome.get_chromosome("chr5").length, 179185813)
        self.assertEqual(genome.get_chromosome("chr6").length, 169501136)
        self.assertEqual(genome.get_chromosome("chr7").length, 145408105)
        self.assertEqual(genome.get_chromosome("chr8").length, 144036388)
        self.assertEqual(genome.get_chromosome("chr9").length, 112206110)
        self.assertEqual(genome.get_chromosome("chr10").length, 132178492)
        self.assertEqual(genome.get_chromosome("chr11").length, 128122151)
        self.assertEqual(genome.get_chromosome("chr12").length, 132184051)
        self.assertEqual(genome.get_chromosome("chr13").length, 98475126)
        self.assertEqual(genome.get_chromosome("chr14").length, 88963417)
        self.assertEqual(genome.get_chromosome("chr15").length, 82547911)
        self.assertEqual(genome.get_chromosome("chr16").length, 68237989)
        self.assertEqual(genome.get_chromosome("chr17").length, 75914007)
        self.assertEqual(genome.get_chromosome("chr18").length, 75923960)
        self.assertEqual(genome.get_chromosome("chr19").length, 57575784)
        self.assertEqual(genome.get_chromosome("chr20").length, 60841859)
        self.assertEqual(genome.get_chromosome("chr21").length, 34683425)
        self.assertEqual(genome.get_chromosome("chr22").length, 35308119)
        self.assertEqual(genome.get_chromosome("chrX").length, 151242693)


class TestPongo(unittest.TestCase, test_models.QcdModelTestMixin):
    model = pongo._orangutan()
    qc_model = pongo_qc.LockePongo()
