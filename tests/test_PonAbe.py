import stdpopsim
from tests import test_species


class TestSpecies(test_species.SpeciesTestBase):
    species = stdpopsim.get_species("PonAbe")

    def test_basic_attributes(self):
        assert self.species.population_size == 1.79 * 10**4
        assert self.species.generation_time == 20


class TestGenome(test_species.GenomeTestBase):
    """
    Tests for the Pongo abelii genome.
    """

    genome = stdpopsim.get_species("PonAbe").genome

    def test_basic_attributes(self):
        # 23 + X + MT
        assert len(self.genome.chromosomes) == 25
