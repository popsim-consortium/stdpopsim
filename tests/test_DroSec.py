import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("DroSec")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "drosophila_sechellia"

    def test_name(self):
        assert self.species.name == "Drosophila sechellia"

    def test_common_name(self):
        assert self.species.common_name == "Drosophila sechellia"

    def test_qc_population_size(self):
        assert self.species.population_size == 1e5

    def test_qc_generation_time(self):
        assert self.species.generation_time == 1 / 20


class TestGenomeData(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("DroSec").genome

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "2L": 2.28e-8,
            "2R": 2.51e-8,
            "3L": 1.88e-8,
            "3R": 1.91e-8,
            "4": 0,
            "X": 2.85e-8,
        }.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(
            self.genome.get_chromosome(name).recombination_rate
        )

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "2L": 1.5e-9,
            "2R": 1.5e-9,
            "3L": 1.5e-9,
            "3R": 1.5e-9,
            "4": 1.5e-9,
            "X": 1.5e-9,
        }.items(),
    )
    def test_mutation_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).mutation_rate)

    @pytest.mark.parametrize("chrom", [chrom for chrom in genome.chromosomes])
    def test_chromosome_ploidy(self, chrom):
        assert chrom.ploidy == 2

    def test_assembly_source(self):
        assert self.genome.assembly_source == "manual"

    def test_assembly_build_version(self):
        assert self.genome.assembly_build_version is None
