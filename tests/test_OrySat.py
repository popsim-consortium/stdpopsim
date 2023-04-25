import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):
    species = stdpopsim.get_species("OrySat")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "oryza_sativa"

    def test_name(self):
        assert self.species.name == "Oryza sativa"

    def test_common_name(self):
        assert self.species.common_name == "Asian rice"

    def test_qc_population_size(self):
        assert self.species.population_size == 46875

    def test_qc_generation_time(self):
        assert self.species.generation_time == 1


class TestGenomeData(test_species.GenomeTestBase):
    genome = stdpopsim.get_species("OrySat").genome

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": 8.97e-10,
            "Mt": 0,
            "Pt": 0,
            "2": 8.97e-10,
            "3": 8.97e-10,
            "4": 8.97e-10,
            "5": 8.97e-10,
            "6": 8.97e-10,
            "7": 8.97e-10,
            "8": 8.97e-10,
            "9": 8.97e-10,
            "10": 8.97e-10,
            "11": 8.97e-10,
            "12": 8.97e-10,
        }.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(
            self.genome.get_chromosome(name).recombination_rate
        )

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": 3.2e-9,
            "Mt": 3.2e-9,
            "Pt": 3.2e-9,
            "2": 3.2e-9,
            "3": 3.2e-9,
            "4": 3.2e-9,
            "5": 3.2e-9,
            "6": 3.2e-9,
            "7": 3.2e-9,
            "8": 3.2e-9,
            "9": 3.2e-9,
            "10": 3.2e-9,
            "11": 3.2e-9,
            "12": 3.2e-9,
        }.items(),
    )
    
    def test_mutation_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).mutation_rate)
    
    @pytest.mark.parametrize("chrom", [chrom for chrom in genome.chromosomes]) 
    def test_chromosome_ploidy(self, chrom): 
        if chrom.id in ["Mt", "Pt"]: 
            assert chrom.ploidy == 1 
        else: 
            assert chrom.ploidy == 2
