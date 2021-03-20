import pytest

import stdpopsim
from tests import test_species


class TestSpecies(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("BosTau")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "bos_taurus"

    def test_name(self):
        assert self.species.name == "Bos taurus"

    def test_common_name(self):
        assert self.species.common_name == "Cattle"

    def test_qc_population_size(self):
        assert self.species.population_size == 90  # most recent Ne in _MacLeodEtAl

    def test_qc_generation_time(self):
        assert self.species.generation_time == 5


class TestGenome(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("BosTau").genome

    def test_basic_attributes(self):
        nchrom = 31  # 29 + X + MT
        assert len(self.genome.chromosomes) == nchrom

    @pytest.mark.skip("Recombination rate QC not done yet")
    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": -1,
            "2": -1,
            "3": -1,
            "4": -1,
            "5": -1,
            "6": -1,
            "7": -1,
            "8": -1,
            "9": -1,
            "10": -1,
            "11": -1,
            "12": -1,
            "13": -1,
            "14": -1,
            "15": -1,
            "16": -1,
            "17": -1,
            "18": -1,
            "19": -1,
            "20": -1,
            "21": -1,
            "22": -1,
            "23": -1,
            "24": -1,
            "25": -1,
            "26": -1,
            "27": -1,
            "28": -1,
            "29": -1,
            "X": -1,
            "MT": -1,
        }.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert pytest.approx(rate, self.genome.get_chromosome(name).recombination_rate)

    @pytest.mark.skip("Mutation rate QC not done yet")
    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": -1,
            "2": -1,
            "3": -1,
            "4": -1,
            "5": -1,
            "6": -1,
            "7": -1,
            "8": -1,
            "9": -1,
            "10": -1,
            "11": -1,
            "12": -1,
            "13": -1,
            "14": -1,
            "15": -1,
            "16": -1,
            "17": -1,
            "18": -1,
            "19": -1,
            "20": -1,
            "21": -1,
            "22": -1,
            "23": -1,
            "24": -1,
            "25": -1,
            "26": -1,
            "27": -1,
            "28": -1,
            "29": -1,
            "X": -1,
            "MT": 0,
        }.items(),
    )
    def test_mutation_rate(self, name, rate):
        assert pytest.approx(rate, self.genome.get_chromosome(name).mutation_rate)
