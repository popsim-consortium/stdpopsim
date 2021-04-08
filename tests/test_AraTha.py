import stdpopsim
from tests import test_species


class TestSpecies(test_species.SpeciesTestBase):
    species = stdpopsim.get_species("AraTha")

    def test_basic_attributes(self):
        assert self.species.population_size == 10 ** 4
        assert self.species.generation_time == 1


class TestGenome(test_species.GenomeTestBase):
    """
    Tests for the arabidopsis_thaliana genome.
    """

    genome = stdpopsim.get_species("AraTha").genome

    def test_basic_attributes(self):
        # 5 autosomes + Mt + Pt
        assert len(self.genome.chromosomes) == 7
