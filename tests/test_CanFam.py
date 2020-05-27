import os
import unittest
from unittest import mock

import stdpopsim
import stdpopsim.cli
from tests import test_species


class TestSpecies(unittest.TestCase, test_species.SpeciesTestMixin):
    species = stdpopsim.get_species("CanFam")

    def test_basic_attributes(self):
        self.assertEqual(self.species.population_size, 13000)
        self.assertEqual(self.species.generation_time, 3)


class TestGenome(unittest.TestCase, test_species.GenomeTestMixin):
    genome = stdpopsim.get_species("CanFam").genome

    def test_basic_attributes(self):
        nchrom = 39  # 38 + X
        self.assertEqual(len(self.genome.chromosomes), nchrom)


class TestThatDogsCanBeSimulated(unittest.TestCase):
    def test_basic_cli_usage(self):
        cmd = "-q CanFam -c chr38 -l 0.001 --seed 1234 10"
        with mock.patch("sys.stdout", autospec=True) as stdout:
            stdout.buffer = open(os.devnull, "wb")
            stdpopsim.cli.stdpopsim_main(cmd.split())
