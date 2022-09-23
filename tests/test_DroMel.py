"""
Tests for the drosophila_melanogaster data definitions.
"""
import pytest

import stdpopsim
from tests import test_species


class TestSpecies(test_species.SpeciesTestBase):
    species = stdpopsim.get_species("DroMel")

    def test_basic_attributes(self):
        self.species.population_size == 1720600
        self.species.generation_time == 0.1


class TestGenome(test_species.GenomeTestBase):
    """
    Tests for the drosophila_melanogaster genome.
    """

    genome = stdpopsim.get_species("DroMel").genome

    def test_basic_attributes(self):
        assert len(self.genome.chromosomes) == 8

    @pytest.mark.parametrize("chr_id", ["Y", "mitochondrion_genome", "4"])
    def test_non_recombining_chrs(self, chr_id):
        chrom = self.genome.get_chromosome(chr_id)
        assert chrom.recombination_rate == 0
        assert chrom.gene_conversion_fraction == 0

    @pytest.mark.parametrize("chr_id", ["2L", "2R", "3L", "3R", "X"])
    def test_recombining_chrs(self, chr_id):
        chrom = self.genome.get_chromosome(chr_id)
        assert chrom.recombination_rate > 0
        assert chrom.gene_conversion_fraction > 0

    @pytest.mark.parametrize("chrom", [chrom for chrom in genome.chromosomes])
    def test_chromosome_ploidy(self, chrom):
        if chrom.id in ["mitochondrion_genome", "Y"]:
            assert chrom.ploidy == 1
        else:
            assert chrom.ploidy == 2
