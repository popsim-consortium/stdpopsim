import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("MusMus")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "mus_musculus"

    def test_name(self):
        assert self.species.name == "Mus musculus"

    def test_common_name(self):
        assert self.species.common_name == "Mouse"

    # QC Tests. These tests are performed by another contributor
    # independently referring to the citations provided in the
    # species definition, filling in the appropriate values
    # and deleting the pytest "skip" annotations.
    # @pytest.mark.skip("Population size QC not done yet")
    def test_qc_population_size(self):
        assert self.species.population_size == 500000

    # @pytest.mark.skip("Generation time QC not done yet")
    def test_qc_generation_time(self):
        assert self.species.generation_time == 0.75


class TestGenomeData(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("MusMus").genome

    # @pytest.mark.skip("Recombination rate QC not done yet")
    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": 0.50e-8,
            "2": 0.57e-8,
            "3": 0.52e-8,
            "4": 0.56e-8,
            "5": 0.59e-8,
            "6": 0.53e-8,
            "7": 0.58e-8,
            "8": 0.58e-8,
            "9": 0.61e-8,
            "10": 0.61e-8,
            "11": 0.70e-8,
            "12": 0.53e-8,
            "13": 0.56e-8,
            "14": 0.53e-8,
            "15": 0.56e-8,
            "16": 0.59e-8,
            "17": 0.65e-8,
            "18": 0.66e-8,
            "19": 0.94e-8,
            "X": 0.48e-8,
            "Y": 0.0,
            "MT": 0.0,
        }.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(
            self.genome.get_chromosome(name).recombination_rate
        )

    # @pytest.mark.skip("Mutation rate QC not done yet")
    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": 5.4e-9,
            "2": 5.4e-9,
            "3": 5.4e-9,
            "4": 5.4e-9,
            "5": 5.4e-9,
            "6": 5.4e-9,
            "7": 5.4e-9,
            "8": 5.4e-9,
            "9": 5.4e-9,
            "10": 5.4e-9,
            "11": 5.4e-9,
            "12": 5.4e-9,
            "13": 5.4e-9,
            "14": 5.4e-9,
            "15": 5.4e-9,
            "16": 5.4e-9,
            "17": 5.4e-9,
            "18": 5.4e-9,
            "19": 5.4e-9,
            "X": 5.4e-9,
            "Y": 5.4e-9,
            "MT": 3.7e-8 * 0.75,
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
            "X": 2,
            "Y": 1,
            "MT": 1,
        }.items(),
    )
    def test_chromosome_ploidy(self, name, ploidy):
        assert ploidy == self.genome.get_chromosome(name).ploidy
