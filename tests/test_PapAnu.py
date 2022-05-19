import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("PapAnu")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "papio_anubis"

    def test_name(self):
        assert self.species.name == "Papio anubis"

    def test_common_name(self):
        assert self.species.common_name == "Olive baboon"

    def test_qc_population_size(self):
        assert self.species.population_size == 335505

    def test_qc_generation_time(self):
        assert self.species.generation_time == 11


class TestGenomeData(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("PapAnu").genome

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": 9.926379e-09,
            "2": 9.605435e-09,
            "3": 9.022377e-09,
            "4": 9.825128e-09,
            "5": 9.579804e-09,
            "6": 1.049788e-08,
            "7": 1.118884e-08,
            "8": 1.108988e-08,
            "9": 1.132883e-08,
            "10": 1.175322e-08,
            "11": 1.184026e-08,
            "12": 1.082400e-08,
            "13": 1.246772e-08,
            "14": 1.274188e-08,
            "15": 1.260836e-08,
            "16": 1.476158e-08,
            "17": 1.524101e-08,
            "18": 1.368410e-08,
            "19": 1.303735e-08,
            "20": 1.677201e-08,
            "X": 1.18898e-08,
            "Y": 0.0,
        }.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(
            self.genome.get_chromosome(name).recombination_rate
        )

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": 5.7e-9,
            "2": 5.7e-9,
            "3": 5.7e-9,
            "4": 5.7e-9,
            "5": 5.7e-9,
            "6": 5.7e-9,
            "7": 5.7e-9,
            "8": 5.7e-9,
            "9": 5.7e-9,
            "10": 5.7e-9,
            "11": 5.7e-9,
            "12": 5.7e-9,
            "13": 5.7e-9,
            "14": 5.7e-9,
            "15": 5.7e-9,
            "16": 5.7e-9,
            "17": 5.7e-9,
            "18": 5.7e-9,
            "19": 5.7e-9,
            "20": 5.7e-9,
            "X": 5.7e-9,
            "Y": 5.7e-9,
        }.items(),
    )
    def test_mutation_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).mutation_rate)
