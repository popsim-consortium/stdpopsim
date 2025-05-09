import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("AnoCar")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "anolis_carolinensis"

    def test_name(self):
        assert self.species.name == "Anolis carolinensis"

    def test_common_name(self):
        assert self.species.common_name == "Anole lizard"

    def test_qc_population_size(self):
        assert self.species.population_size == 3.05e6

    def test_qc_generation_time(self):
        assert self.species.generation_time == 1.5


class TestGenomeData(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("AnoCar").genome

    rec_rate = 1.59e-8  # placeholder as we wait for map

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": rec_rate,
            "2": rec_rate,
            "3": rec_rate,
            "4": rec_rate,
            "5": rec_rate,
            "6": rec_rate,
            "LGa": rec_rate,
            "LGb": rec_rate,
            "LGc": rec_rate,
            "LGd": rec_rate,
            "LGf": rec_rate,
            "LGg": rec_rate,
            "LGh": rec_rate,
            "MT": 0,
        }.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(
            self.genome.get_chromosome(name).recombination_rate
        )

    mu = 2.1e-10

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": mu,
            "2": mu,
            "3": mu,
            "4": mu,
            "5": mu,
            "6": mu,
            "LGa": mu,
            "LGb": mu,
            "LGc": mu,
            "LGd": mu,
            "LGf": mu,
            "LGg": mu,
            "LGh": mu,
            "MT": mu,
        }.items(),
    )
    def test_mutation_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).mutation_rate)

    @pytest.mark.parametrize("chrom", [chrom for chrom in genome.chromosomes])
    def test_chromosome_ploidy(self, chrom):
        if chrom.id in ["MT"]:
            assert chrom.ploidy == 1
        else:
            assert chrom.ploidy == 2

    def test_assembly_source(self):
        assert self.genome.assembly_source == "ensembl"

    def test_assembly_build_version(self):
        assert self.genome.assembly_build_version is not None
