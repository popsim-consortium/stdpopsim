import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("HelAnn")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "helianthus_annuus"

    def test_name(self):
        assert self.species.name == "Helianthus annuus"

    def test_common_name(self):
        assert self.species.common_name == "Helianthus annuus"

    def test_qc_population_size(self):
        assert self.species.population_size == 673968

    def test_qc_generation_time(self):
        assert self.species.generation_time == 1


class TestGenomeData(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("HelAnn").genome

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": 0.4e-8,
            "2": 0.4e-8,
            "3": 0.4e-8,
            "4": 0.4e-8,
            "5": 0.4e-8,
            "6": 0.4e-8,
            "7": 0.4e-8,
            "8": 0.4e-8,
            "9": 0.4e-8,
            "10": 0.4e-8,
            "11": 0.4e-8,
            "12": 0.4e-8,
            "13": 0.4e-8,
            "14": 0.4e-8,
            "15": 0.4e-8,
            "16": 0.4e-8,
            "17": 0.4e-8,
        }.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(
            self.genome.get_chromosome(name).recombination_rate
        )

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": 6.1e-9,
            "2": 6.1e-9,
            "3": 6.1e-9,
            "4": 6.1e-9,
            "5": 6.1e-9,
            "6": 6.1e-9,
            "7": 6.1e-9,
            "8": 6.1e-9,
            "9": 6.1e-9,
            "10": 6.1e-9,
            "11": 6.1e-9,
            "12": 6.1e-9,
            "13": 6.1e-9,
            "14": 6.1e-9,
            "15": 6.1e-9,
            "16": 6.1e-9,
            "17": 6.1e-9,
        }.items(),
    )
    def test_mutation_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).mutation_rate)

    @pytest.mark.parametrize("chrom", [chrom for chrom in genome.chromosomes])
    @pytest.mark.filterwarnings("ignore::stdpopsim.NonAutosomalWarning")
    def test_chromosome_ploidy(self, chrom):
        assert chrom.ploidy == 2
