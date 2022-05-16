import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("GasAcu")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "9307941"

    def test_name(self):
        assert self.species.name == "Gasterosteus aculeatus"

    def test_common_name(self):
        assert self.species.common_name == "Three-spined stickleback"

    # QC Tests. These tests are performed by another contributor
    # independently referring to the citations provided in the
    # species definition, filling in the appropriate values
    # and deleting the pytest "skip" annotations.
    def test_qc_population_size(self):
        assert self.species.population_size == 1e4

    def test_qc_generation_time(self):
        assert self.species.generation_time == 2


class TestGenomeData(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("GasAcu").genome

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": 3.11e-08,
            "2": 3.11e-08,
            "3": 3.11e-08,
            "4": 3.11e-08,
            "5": 3.11e-08,
            "6": 3.11e-08,
            "7": 3.11e-08,
            "8": 3.11e-08,
            "9": 3.11e-08,
            "10": 3.11e-08,
            "11": 3.11e-08,
            "12": 3.11e-08,
            "13": 3.11e-08,
            "14": 3.11e-08,
            "15": 3.11e-08,
            "16": 3.11e-08,
            "17": 3.11e-08,
            "18": 3.11e-08,
            "19": 3.11e-08,
            "20": 3.11e-08,
            "21": 3.11e-08,
            "Y": 0,
            "MT": 0,
        }.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(
            self.genome.get_chromosome(name).recombination_rate
        )

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": 3.7e-8,
            "2": 3.7e-8,
            "3": 3.7e-8,
            "4": 3.7e-8,
            "5": 3.7e-8,
            "6": 3.7e-8,
            "7": 3.7e-8,
            "8": 3.7e-8,
            "9": 3.7e-8,
            "10": 3.7e-8,
            "11": 3.7e-8,
            "12": 3.7e-8,
            "13": 3.7e-8,
            "14": 3.7e-8,
            "15": 3.7e-8,
            "16": 3.7e-8,
            "17": 3.7e-8,
            "18": 3.7e-8,
            "19": 3.7e-8,
            "20": 3.7e-8,
            "21": 3.7e-8,
            "Y": 3.7e-8,
            "MT": 3.7e-8,
        }.items(),
    )
    def test_mutation_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).mutation_rate)
