import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("AnaPla")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "anas_platyrhynchos"

    def test_name(self):
        assert self.species.name == "Anas platyrhynchos"

    def test_common_name(self):
        assert self.species.common_name == "Mallard"

    def test_qc_population_size(self):
        assert self.species.population_size == 156000

    def test_qc_generation_time(self):
        assert self.species.generation_time == 4


class TestGenomeData(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("AnaPla").genome

    # total map lengths for chromosomes 1-18 taken from from Huang et al 2006
    # Table 1. See: https://academic.oup.com/view-large/280279985
    # The average rec rate for chroms 1-18 is this map length
    # divided by the physical length of the chromosome
    # The average rec rate of remaining chroms should be the total
    # average rate across chroms 1-18
    huang_chroms = {
        "1": {"M": 3.170},
        "2": {"M": 2.260},
        "3": {"M": 1.123},
        "4": {"M": 0.934},
        "5": {"M": 0.786},
        "6": {"M": 1.201},
        "7": {"M": 0.981},
        "8": {"M": 0.427},
        "9": {"M": 0.359},
        "10": {"M": 0.574},
        "11": {"M": 0.333},
        "12": {"M": 0.703},
        "13": {"M": 0.307},
        "14": {"M": 0.175},
        "15": {"M": 0.069},
        "16": {"M": 0.062},
        "17": {"M": 0.048},
        "18": {"M": 0.021},
    }
    for c in huang_chroms:
        h = huang_chroms[c]
        h["bp"] = genome.get_chromosome(c).length
        h["rate"] = round(h["M"] / h["bp"], 11)
    total_M = sum([x["M"] for x in huang_chroms.values()])
    total_bp = sum([x["bp"] for x in huang_chroms.values()])
    mean_rate = round(total_M / total_bp, 10)

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": huang_chroms["1"]["rate"],
            "2": huang_chroms["2"]["rate"],
            "3": huang_chroms["3"]["rate"],
            "4": huang_chroms["4"]["rate"],
            "5": huang_chroms["5"]["rate"],
            "6": huang_chroms["6"]["rate"],
            "7": huang_chroms["7"]["rate"],
            "8": huang_chroms["8"]["rate"],
            "9": huang_chroms["9"]["rate"],
            "10": huang_chroms["10"]["rate"],
            "11": huang_chroms["11"]["rate"],
            "12": huang_chroms["12"]["rate"],
            "13": huang_chroms["13"]["rate"],
            "14": huang_chroms["14"]["rate"],
            "15": huang_chroms["15"]["rate"],
            "16": huang_chroms["16"]["rate"],
            "17": huang_chroms["17"]["rate"],
            "18": huang_chroms["18"]["rate"],
            "19": mean_rate,
            "20": mean_rate,
            "21": mean_rate,
            "22": mean_rate,
            "23": mean_rate,
            "24": mean_rate,
            "25": mean_rate,
            "26": mean_rate,
            "27": mean_rate,
            "28": mean_rate,
            "29": mean_rate,
            "30": mean_rate,
            "Z": mean_rate,
            "31": mean_rate,
            "32": mean_rate,
            "33": mean_rate,
            "34": mean_rate,
            "35": mean_rate,
            "36": mean_rate,
            "37": mean_rate,
            "38": mean_rate,
            "39": mean_rate,
            "40": mean_rate,
        }.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(
            self.genome.get_chromosome(name).recombination_rate, 0.01
        )

    _genome_mutation_rate = 4.83e-9

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": _genome_mutation_rate,
            "2": _genome_mutation_rate,
            "3": _genome_mutation_rate,
            "4": _genome_mutation_rate,
            "5": _genome_mutation_rate,
            "6": _genome_mutation_rate,
            "7": _genome_mutation_rate,
            "8": _genome_mutation_rate,
            "9": _genome_mutation_rate,
            "10": _genome_mutation_rate,
            "11": _genome_mutation_rate,
            "12": _genome_mutation_rate,
            "13": _genome_mutation_rate,
            "14": _genome_mutation_rate,
            "15": _genome_mutation_rate,
            "16": _genome_mutation_rate,
            "17": _genome_mutation_rate,
            "18": _genome_mutation_rate,
            "19": _genome_mutation_rate,
            "20": _genome_mutation_rate,
            "21": _genome_mutation_rate,
            "22": _genome_mutation_rate,
            "23": _genome_mutation_rate,
            "24": _genome_mutation_rate,
            "25": _genome_mutation_rate,
            "26": _genome_mutation_rate,
            "27": _genome_mutation_rate,
            "28": _genome_mutation_rate,
            "29": _genome_mutation_rate,
            "30": _genome_mutation_rate,
            "Z": _genome_mutation_rate,
            "31": _genome_mutation_rate,
            "32": _genome_mutation_rate,
            "33": _genome_mutation_rate,
            "34": _genome_mutation_rate,
            "35": _genome_mutation_rate,
            "36": _genome_mutation_rate,
            "37": _genome_mutation_rate,
            "38": _genome_mutation_rate,
            "39": _genome_mutation_rate,
            "40": _genome_mutation_rate,
        }.items(),
    )
    def test_mutation_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).mutation_rate)

    def test_bacterial_recombination(self):
        assert self.genome.bacterial_recombination is False
