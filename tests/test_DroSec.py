import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("DroSec")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "drosophila_sechellia"

    def test_name(self):
        assert self.species.name == "Drosophila sechellia"

    def test_common_name(self):
        assert self.species.common_name == "Drosophila sechellia"

    # QC Tests. These tests are performed by another contributor
    # independently referring to the citations provided in the
    # species definition, filling in the appropriate values
    # and deleting the pytest "skip" annotations.
    def test_qc_population_size(self):
        assert self.species.population_size == 100000

    @pytest.mark.skip("Generation time QC not done yet")
    def test_qc_generation_time(self):
        assert self.species.generation_time == 20


class TestGenomeData(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("DroSec").genome

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "2L": 2.39,
            "2R": 2.66,
            "3L": 1.79,
            "3R": 1.96,
            "X": 2.95,
            "4": 0,
        }.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(
            self.genome.get_chromosome(name).recombination_rate
        )

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "2L": 1.5e-09,
            "2R": 1.5e-09,
            "3L": 1.5e-09,
            "3R": 1.5e-09,
            "X": 1.5e-09,
            "4": 1.5e-09,
        }.items(),
    )
    def test_mutation_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).mutation_rate)
