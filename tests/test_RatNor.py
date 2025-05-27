import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("RatNor")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "rattus_norvegicus"

    def test_name(self):
        assert self.species.name == "Rattus norvegicus"

    def test_common_name(self):
        assert self.species.common_name == "Rat"

    def test_qc_population_size(self):
        assert self.species.population_size == 1.24e5  # 4D sites

    def test_qc_generation_time(self):
        assert self.species.generation_time == 0.5


class TestGenomeData(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("RatNor").genome

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": 5.0e-09,
            "2": 4.8e-09,
            "3": 5.9e-09,
            "4": 5.7e-09,
            "5": 6.1e-09,
            "6": 6.0e-09,
            "7": 6.6e-09,
            "8": 6.8e-09,
            "9": 6.4e-09,
            "10": 8.1e-09,
            "11": 6.8e-09,
            "12": 10.1e-09,
            "13": 5.7e-09,
            "14": 6.3e-09,
            "15": 6.5e-09,
            "16": 6.8e-09,
            "17": 7.5e-09,
            "18": 6.7e-09,
            "19": 8.8e-09,
            "20": 9.1e-09,
            "X": 3.7e-09,
            "Y": 0,
            "MT": 0,
        }.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(
            self.genome.get_chromosome(name).recombination_rate
        )

    qc_mut_rate = 2.96e-09  # Deinum et al. 2015

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": qc_mut_rate,
            "2": qc_mut_rate,
            "3": qc_mut_rate,
            "4": qc_mut_rate,
            "5": qc_mut_rate,
            "6": qc_mut_rate,
            "7": qc_mut_rate,
            "8": qc_mut_rate,
            "9": qc_mut_rate,
            "10": qc_mut_rate,
            "11": qc_mut_rate,
            "12": qc_mut_rate,
            "13": qc_mut_rate,
            "14": qc_mut_rate,
            "15": qc_mut_rate,
            "16": qc_mut_rate,
            "17": qc_mut_rate,
            "18": qc_mut_rate,
            "19": qc_mut_rate,
            "20": qc_mut_rate,
            "X": qc_mut_rate,
            "Y": qc_mut_rate,
            "MT": qc_mut_rate * 15,  # mitochondrial mutation rate
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
            "X": 2,
            "Y": 1,
            "MT": 1,
        }.items(),
    )
    def test_chromosome_ploidy(self, name, ploidy):
        assert ploidy == self.genome.get_chromosome(name).ploidy
