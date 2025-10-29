"""
Tests for the EscCol data definitions.
"""

import stdpopsim
import pytest
from tests import test_species


class TestSpecies(test_species.SpeciesTestBase):
    species = stdpopsim.get_species("EscCol")

    def test_basic_attributes(self):
        # From paper https://doi.org/10.1093/molbev/msw048
        # Ne taken from Table 2
        assert self.species.population_size == 1.8e8
        # 20 minutes per generation
        generation_time = 1.0 / (525600 / 20)
        assert round(abs(self.species.generation_time - generation_time), 7) == 0
        assert self.species.separate_sexes is False


class TestGenome(test_species.GenomeTestBase):
    """
    Tests for the EscCol genome.
    """

    genome = stdpopsim.get_species("EscCol").genome

    def test_basic_attributes(self):
        assert len(self.genome.chromosomes) == 1

    def test_mutation_rate(self):
        assert self.genome.get_chromosome("Chromosome").mutation_rate == 8.9e-11

    def test_recombination_rate(self):
        assert self.genome.get_chromosome("Chromosome").recombination_rate == 8.9e-11

    def test_gene_conversion_length(self):
        assert self.genome.get_chromosome("Chromosome").gene_conversion_length == 542

    def test_bacterial_recombination(self):
        assert self.genome.bacterial_recombination is True

    @pytest.mark.parametrize("chrom", [chrom for chrom in genome.chromosomes])
    def test_chromosome_ploidy(self, chrom):
        assert chrom.ploidy == 1

    def test_assembly_source(self):
        assert self.genome.assembly_source == "ensembl"

    def test_assembly_build_version(self):
        assert self.genome.assembly_build_version is not None
