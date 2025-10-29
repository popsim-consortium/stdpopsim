import stdpopsim

import pytest

from tests import test_species


class TestSpecies(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("PanTro")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "pan_troglodytes"

    def test_name(self):
        assert self.species.name == "Pan troglodytes"

    def test_common_name(self):
        assert self.species.common_name == "Chimpanzee"

    def test_qc_population_size(self):
        assert self.species.population_size == 16781

    def test_qc_generation_time(self):
        assert self.species.generation_time == 24.6

    def test_separate_sexes(self):
        assert self.species.separate_sexes is True


class TestGenome(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("PanTro").genome

    def test_basic_attributes(self):
        assert len(self.genome.chromosomes) == 25

    @pytest.mark.parametrize("chr_id", [chrom.id for chrom in genome.chromosomes])
    def test_recombination_rate(self, chr_id, rate=1.2e-8):
        assert rate == pytest.approx(
            self.genome.get_chromosome(chr_id).recombination_rate
        )

    @pytest.mark.parametrize("chr_id", [chrom.id for chrom in genome.chromosomes])
    def test_mutation_rate(self, chr_id, rate=1.6e-8):
        assert rate == pytest.approx(self.genome.get_chromosome(chr_id).mutation_rate)

    @pytest.mark.parametrize("chrom", [chrom for chrom in genome.chromosomes])
    def test_chromosome_ploidy(self, chrom):
        if chrom.id in ["Y"]:
            assert chrom.ploidy == 1
        else:
            assert chrom.ploidy == 2
