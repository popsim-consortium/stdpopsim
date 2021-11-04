import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("HelMel")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "heliconius_melpomene"

    def test_name(self):
        assert self.species.name == "Heliconius melpomene"

    def test_common_name(self):
        assert self.species.common_name == "Heliconius melpomene"

    # QC Tests. These tests are performed by another contributor
    # independently referring to the citations provided in the
    # species definition, filling in the appropriate values
    # and deleting the pytest "skip" annotations.
    def test_qc_population_size(self):
        assert self.species.population_size == 2.1e06

    def test_qc_generation_time(self):
        assert self.species.generation_time == 1 / 10


class TestGenomeData(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("HelMel").genome

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": 3.17e-08,
            "2": 5.61e-08,
            "3": 5.10e-08,
            "4": 4.97e-08,
            "5": 5.15e-08,
            "6": 3.40e-08,
            "7": 3.76e-08,
            "8": 5.28e-08,
            "9": 5.31e-08,
            "10": 3.16e-08,
            "11": 4.47e-08,
            "12": 3.13e-08,
            "13": 3.08e-08,
            "14": 5.47e-08,
            "15": 4.78e-08,
            "16": 4.71e-08,
            "17": 3.94e-08,
            "18": 3.16e-08,
            "19": 3.11e-08,
            "20": 3.45e-08,
            "21": 3.71e-08,
        }.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(
            self.genome.get_chromosome(name).recombination_rate
        )

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": 2.9e-09,
            "2": 2.9e-09,
            "3": 2.9e-09,
            "4": 2.9e-09,
            "5": 2.9e-09,
            "6": 2.9e-09,
            "7": 2.9e-09,
            "8": 2.9e-09,
            "9": 2.9e-09,
            "10": 2.9e-09,
            "11": 2.9e-09,
            "12": 2.9e-09,
            "13": 2.9e-09,
            "14": 2.9e-09,
            "15": 2.9e-09,
            "16": 2.9e-09,
            "17": 2.9e-09,
            "18": 2.9e-09,
            "19": 2.9e-09,
            "20": 2.9e-09,
            "21": 2.9e-09,
        }.items(),
    )
    def test_mutation_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).mutation_rate)
