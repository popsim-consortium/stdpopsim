"""
Tests for the human data definitions.
"""

import stdpopsim

import pytest

from tests import test_species


class TestSpecies(test_species.SpeciesTestBase):
    species = stdpopsim.get_species("HomSap")

    def test_basic_attributes(self):
        assert self.species.population_size == 10**4
        assert self.species.generation_time == 30


class TestGenome(test_species.GenomeTestBase):

    species = stdpopsim.get_species("HomSap")
    genome = species.genome

    def test_basic_attributes(self):
        assert len(self.genome.chromosomes) == 25

    @pytest.mark.parametrize("chr_id", [chrom.id for chrom in genome.chromosomes])
    def test_recombination_rates(self, chr_id):
        genetic_map = "HapMapII_GRCh38"
        chrom = self.genome.get_chromosome(chr_id)
        if chr_id in ["X", "Y", "MT"]:
            with pytest.warns(stdpopsim.NonAutosomalWarning):
                contig = self.species.get_contig(chr_id, genetic_map=genetic_map)
        else:
            contig = self.species.get_contig(chr_id, genetic_map=genetic_map)
        assert chrom.length == contig.recombination_map.sequence_length
        assert chrom.recombination_rate == pytest.approx(
            contig.recombination_map.mean_rate, rel=1e-6
        )

    def test_bacterial_recombination(self):
        assert self.genome.bacterial_recombination is False

    @pytest.mark.parametrize("chrom", [chrom for chrom in genome.chromosomes])
    def test_chromosome_ploidy(self, chrom):
        if chrom.id in ["MT", "Y"]:
            assert chrom.ploidy == 1
        else:
            assert chrom.ploidy == 2
