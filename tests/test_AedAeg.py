import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("AedAeg")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "aedes_aegypti_lvpagwg"

    def test_name(self):
        assert self.species.name == "Aedes aegypti"

    def test_common_name(self):
        assert self.species.common_name == "Yellow fever mosquito"

    def test_qc_population_size(self):
        assert self.species.population_size == 1000000

    def test_qc_generation_time(self):
        assert self.species.generation_time == 1 / 15


class TestGenomeData(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("AedAeg").genome

    @pytest.mark.parametrize(
        ["name", "rate"],
        {"1": 0.306e-8, "2": 0.249e-8, "3": 0.291e-8, "MT": 0.0}.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(
            self.genome.get_chromosome(name).recombination_rate
        )

    @pytest.mark.parametrize(
        ["name", "rate"], {"1": 3.5e-9, "2": 3.5e-9, "3": 3.5e-9, "MT": 3.5e-9}.items()
    )
    def test_mutation_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).mutation_rate)

    @pytest.mark.parametrize("chrom", [chrom for chrom in genome.chromosomes])
    def test_chromosome_ploidy(self, chrom):
        if chrom.id in ["MT"]:
            assert chrom.ploidy == 1
        else:
            assert chrom.ploidy == 2
