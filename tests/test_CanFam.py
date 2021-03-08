import stdpopsim
from tests import test_species


class TestSpecies(test_species.SpeciesTestBase):
    species = stdpopsim.get_species("CanFam")

    def test_basic_attributes(self):
        assert self.species.population_size == 13000
        assert self.species.generation_time == 3


class TestGenome(test_species.GenomeTestBase):
    genome = stdpopsim.get_species("CanFam").genome

    def test_basic_attributes(self):
        nchrom = 40  # 38 + X + MT
        assert len(self.genome.chromosomes) == nchrom
