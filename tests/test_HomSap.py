"""
Tests for the human data definitions.
"""
import stdpopsim

import pytest

from tests import test_species


class TestSpecies(test_species.SpeciesTestBase):
    species = stdpopsim.get_species("HomSap")

    def test_basic_attributes(self):
        assert self.species.population_size == 10 ** 4
        assert self.species.generation_time == 30


class TestGenome(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("HomSap").genome

    def test_basic_attributes(self):
        assert len(self.genome.chromosomes) == 25

    @pytest.mark.parametrize("chr_id", [chrom.id for chrom in genome.chromosomes])
    def test_recombination_rates(self, chr_id):
        # We should recast this test and just hard code in the values.
        # Tests should be *obvious* not clever.

        # recompute recombination rates from HapMapII_GRCh37 map then
        # compare the results to the current recombination rates for each chromosome
        genetic_map = "HapMapII_GRCh37"
        species = stdpopsim.get_species("HomSap")
        chrom = species.genome.get_chromosome(chr_id)
        if chr_id in ["X", "Y", "MT"]:
            with pytest.warns(stdpopsim.NonAutosomalWarning):
                contig = species.get_contig(chr_id, genetic_map=genetic_map)
        elif chr_id in ["3", "5", "7", "11", "16", "17", "18", "20"]:
            contig = species.get_contig(chr_id, genetic_map=genetic_map)
        else:
            # The rest of the chromosomes are currently emitting a warning about
            # the mismatch in chromosome lengths because of the fact that we're
            # on 37 for the map. This should be resolved when we start using the
            # lifted over map.
            with pytest.warns(UserWarning, match="longer than chromosome length"):
                contig = species.get_contig(chr_id, genetic_map=genetic_map)
        assert pytest.approx(
            chrom.recombination_rate,
            contig.recombination_map.mean_rate,
        )
