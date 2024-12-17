import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("PhoSin")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "phocoena_sinus"

    def test_name(self):
        assert self.species.name == "Phocoena sinus"

    def test_common_name(self):
        assert self.species.common_name == "Vaquita"

    # QC Tests. These tests are performed by another contributor
    # independently referring to the citations provided in the
    # species definition, filling in the appropriate values
    # and deleting the pytest "skip" annotations.
    def test_qc_population_size(self):
        assert self.species.population_size == 3500

    def test_qc_generation_time(self):
        assert self.species.generation_time == 11.9


class TestGenomeData(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("PhoSin").genome

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": 1e-8,
            "2": 1e-8,
            "3": 1e-8,
            "4": 1e-8,
            "5": 1e-8,
            "6": 1e-8,
            "7": 1e-8,
            "8": 1e-8,
            "9": 1e-8,
            "10": 1e-8,
            "11": 1e-8,
            "12": 1e-8,
            "13": 1e-8,
            "14": 1e-8,
            "15": 1e-8,
            "16": 1e-8,
            "17": 1e-8,
            "18": 1e-8,
            "19": 1e-8,
            "20": 1e-8,
            "21": 1e-8,
            "X": 1e-8,
        }.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(
            self.genome.get_chromosome(name).recombination_rate
        )

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": 5.83e-9,
            "2": 5.83e-9,
            "3": 5.83e-9,
            "4": 5.83e-9,
            "5": 5.83e-9,
            "6": 5.83e-9,
            "7": 5.83e-9,
            "8": 5.83e-9,
            "9": 5.83e-9,
            "10": 5.83e-9,
            "11": 5.83e-9,
            "12": 5.83e-9,
            "13": 5.83e-9,
            "14": 5.83e-9,
            "15": 5.83e-9,
            "16": 5.83e-9,
            "17": 5.83e-9,
            "18": 5.83e-9,
            "19": 5.83e-9,
            "20": 5.83e-9,
            "21": 5.83e-9,
            "X": 5.83e-9,
        }.items(),
    )
    def test_mutation_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).mutation_rate)

    @pytest.mark.parametrize(
        ["name", "ploidy"],
        {
            "1": 2,
            "2": 2,
            "3": 2,
            "4": 2,
            "5": 2,
            "6": 2,
            "7": 2,
            "8": 2,
            "9": 2,
            "10": 2,
            "11": 2,
            "12": 2,
            "13": 2,
            "14": 2,
            "15": 2,
            "16": 2,
            "17": 2,
            "18": 2,
            "19": 2,
            "20": 2,
            "21": 2,
            "X": 2,
        }.items(),
    )
    def test_chromosome_ploidy(self, name, ploidy):
        assert ploidy == self.genome.get_chromosome(name).ploidy
