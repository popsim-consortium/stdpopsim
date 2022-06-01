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
        assert self.species.population_size == 62000  # ancestral Ne in _MacLeodEtAl

    def test_qc_generation_time(self):
        assert self.species.generation_time == 5


class TestGenome(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("BosTau").genome

    def test_basic_attributes(self):
        nchrom = 31  # 29 + X + MT
        assert len(self.genome.chromosomes) == nchrom

    @pytest.mark.parametrize(
        ["name", "length"],
        {
            "1": 158534110,
            "2": 136231102,
            "3": 121005158,
            "4": 120000601,
            "5": 120089316,
            "6": 117806340,
            "7": 110682743,
            "8": 113319770,
            "9": 105454467,
            "10": 103308737,
            "11": 106982474,
            "12": 87216183,
            "13": 83472345,
            "14": 82403003,
            "15": 85007780,
            "16": 81013979,
            "17": 73167244,
            "18": 65820629,
            "19": 63449741,
            "20": 71974595,
            "21": 69862954,
            "22": 60773035,
            "23": 52498615,
            "24": 62317253,
            "25": 42350435,
            "26": 51992305,
            "27": 45612108,
            "28": 45940150,
            "29": 51098607,
            "X": 139009144,
            "MT": 16338,
        }.items(),
    )
    def test_chromosome_lengths(self, name, length):
        assert length == self.genome.get_chromosome(name).length

    # Performs the calculation described in the species file comments
    lengths = [c.length for c in genome.chromosomes]
    total_length = sum(lengths) - genome.get_chromosome("MT").length
    avg_crossovers = (25.5 + 23.2) / 2  # average of male and female numbers
    avg_rate = avg_crossovers / total_length

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": avg_rate,
            "2": avg_rate,
            "3": avg_rate,
            "4": avg_rate,
            "5": avg_rate,
            "6": avg_rate,
            "7": avg_rate,
            "8": avg_rate,
            "9": avg_rate,
            "10": avg_rate,
            "11": avg_rate,
            "12": avg_rate,
            "13": avg_rate,
            "14": avg_rate,
            "15": avg_rate,
            "16": avg_rate,
            "17": avg_rate,
            "18": avg_rate,
            "19": avg_rate,
            "20": avg_rate,
            "21": avg_rate,
            "22": avg_rate,
            "23": avg_rate,
            "24": avg_rate,
            "25": avg_rate,
            "26": avg_rate,
            "27": avg_rate,
            "28": avg_rate,
            "29": avg_rate,
            "X": avg_rate,
            "MT": 0.0,
        }.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(
            self.genome.get_chromosome(name).recombination_rate, rel=0.001
        )

    _mutation_rate = 1.2e-8

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": _mutation_rate,
            "2": _mutation_rate,
            "3": _mutation_rate,
            "4": _mutation_rate,
            "5": _mutation_rate,
            "6": _mutation_rate,
            "7": _mutation_rate,
            "8": _mutation_rate,
            "9": _mutation_rate,
            "10": _mutation_rate,
            "11": _mutation_rate,
            "12": _mutation_rate,
            "13": _mutation_rate,
            "14": _mutation_rate,
            "15": _mutation_rate,
            "16": _mutation_rate,
            "17": _mutation_rate,
            "18": _mutation_rate,
            "19": _mutation_rate,
            "20": _mutation_rate,
            "21": _mutation_rate,
            "22": _mutation_rate,
            "23": _mutation_rate,
            "24": _mutation_rate,
            "25": _mutation_rate,
            "26": _mutation_rate,
            "27": _mutation_rate,
            "28": _mutation_rate,
            "29": _mutation_rate,
            "X": _mutation_rate,
            "MT": _mutation_rate,
        }.items(),
    )
    def test_mutation_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).mutation_rate)
