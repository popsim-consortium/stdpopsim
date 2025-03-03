import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("GorGor")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "gorilla_gorilla"

    def test_name(self):
        assert self.species.name == "Gorilla gorilla"

    def test_common_name(self):
        assert self.species.common_name == "Gorilla"

    # QC Tests. These tests are performed by another contributor
    def test_qc_population_size(self):
        assert self.species.population_size == 25200

    def test_qc_generation_time(self):
        assert self.species.generation_time == 19.0


class TestGenomeData(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("GorGor").genome
    qc_rec_rate = 1.193e-08  # 0.944 /1000 / 4 / 19,785

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": qc_rec_rate,
            "2A": qc_rec_rate,
            "2B": qc_rec_rate,
            "3": qc_rec_rate,
            "4": qc_rec_rate,
            "5": qc_rec_rate,
            "6": qc_rec_rate,
            "7": qc_rec_rate,
            "8": qc_rec_rate,
            "9": qc_rec_rate,
            "10": qc_rec_rate,
            "11": qc_rec_rate,
            "12": qc_rec_rate,
            "13": qc_rec_rate,
            "14": qc_rec_rate,
            "15": qc_rec_rate,
            "16": qc_rec_rate,
            "17": qc_rec_rate,
            "18": qc_rec_rate,
            "19": qc_rec_rate,
            "20": qc_rec_rate,
            "21": qc_rec_rate,
            "22": qc_rec_rate,
            "X": qc_rec_rate,
            "MT": 0,
        }.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(
            self.genome.get_chromosome(name).recombination_rate
        )

    qc_mut_rate = 1.235e-8  # 0.65e-9 mutations/year times gen. time of 19 years

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": qc_mut_rate,
            "2A": qc_mut_rate,
            "2B": qc_mut_rate,
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
            "21": qc_mut_rate,
            "22": qc_mut_rate,
            "X": qc_mut_rate,
            "MT": qc_mut_rate,
        }.items(),
    )
    def test_mutation_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).mutation_rate)
