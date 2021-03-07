"""
Tests for the e. coli data definitions.
"""
import stdpopsim
from tests import test_species


class TestSpecies(test_species.SpeciesTestBase):
    species = stdpopsim.get_species("EscCol")

    def test_basic_attributes(self):
        # From paper https://doi.org/10.1093/molbev/msw048
        # Ne taken from Table 2
        assert self.species.population_size == 1.8e8
        # 20 minutes per generation
        generation_time = 1.0 / (525600 / 20)
        assert round(abs(self.species.generation_time - generation_time), 7) == 0


class TestGenome(test_species.GenomeTestBase):
    """
    Tests for the e_coli genome.
    """

    genome = stdpopsim.get_species("EscCol").genome

    def test_basic_attributes(self):
        assert len(self.genome.chromosomes) == 1
