import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("ChlRei")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "chlamydomonas_reinhardtii"

    def test_name(self):
        assert self.species.name == "Chlamydomonas reinhardtii"

    def test_common_name(self):
        assert self.species.common_name == "Chlamydomonas reinhardtii"

    # QC Tests. These tests are performed by another contributor
    # independently referring to the citations provided in the
    # species definition, filling in the appropriate values
    # and deleting the pytest "skip" annotations.
    def test_qc_population_size(self):
        assert self.species.population_size == 1.4 * 1e-7

    def test_qc_generation_time(self):
        assert self.species.generation_time == 1 / 876


class TestGenomeData(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("ChlRei").genome

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": 1.21e-10,
            "2": 1.49e-10,
            "3": 1.52e-10,
            "4": 1.47e-10,
            "5": 1.70e-10,
            "6": 1.17e-10,
            "7": 8.66e-11,
            "8": 1.39e-10,
            "9": 1.12e-10,
            "10": 1.97e-10,
            "11": 1.63e-10,
            "12": 9.15e-11,
            "13": 1.43e-10,
            "14": 1.90e-10,
            "15": 3.93e-10,
            "16": 1.71e-10,
            "17": 1.83e-10,
        }.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(
            self.genome.get_chromosome(name).recombination_rate
        )

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": 9.74e-10,
            "2": 8.62e-10,
            "3": 9.50e-10,
            "4": 9.66e-10,
            "5": 1.17e-09,
            "6": 9.12e-10,
            "7": 9.14e-10,
            "8": 8.98e-10,
            "9": 9.17e-10,
            "10": 9.27e-10,
            "11": 1.03e-09,
            "12": 9.55e-10,
            "13": 7.56e-10,
            "14": 8.96e-10,
            "15": 6.91e-10,
            "16": 9.59e-10,
            "17": 1.05e-09,
        }.items(),
    )
    def test_mutation_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).mutation_rate)
