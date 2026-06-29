import stdpopsim
import pytest
from tests import test_species


class TestSpecies(test_species.SpeciesTestBase):
    species = stdpopsim.get_species("AraTha")

    def test_basic_attributes(self):
        assert self.species.population_size == 10**4
        assert self.species.generation_time == 1
        assert self.species.separate_sexes is False


class TestGenome(test_species.GenomeTestBase):
    """
    Tests for the arabidopsis_thaliana genome.
    """

    species = stdpopsim.get_species("AraTha")
    genome = species.genome

    def test_basic_attributes(self):
        # 5 autosomes + Mt + Pt
        assert len(self.genome.chromosomes) == 7

    @pytest.mark.parametrize("chrom", [chrom for chrom in genome.chromosomes])
    def test_chromosome_ploidy(self, chrom):
        if chrom.id in ["Mt", "Pt"]:
            assert chrom.ploidy == 1
        else:
            assert chrom.ploidy == 2

    @pytest.mark.parametrize("chr_id", [chrom.id for chrom in genome.chromosomes])
    def test_recombination_selfing_correction(self, chr_id):
        # Platt et al. 2010 https://doi.org/10.1371/journal.pgen.1000843
        selfing_rate = 0.97
        # Nordborg 2000 https://doi.org/10.1093/genetics/154.2.923
        # Caicedo et al. 2007 https://doi.org/10.1371/journal.pgen.0030163
        rate_correction = 1 - selfing_rate / (2 - selfing_rate)

        # Downloaded salome2012_liftover_maps.tar.gz
        # maps = ["salome2012_liftover_maps/TAIR10_lifted_chr1.txt",
        #        "salome2012_liftover_maps/TAIR10_lifted_chr2.txt",
        #        "salome2012_liftover_maps/TAIR10_lifted_chr3.txt",
        #        "salome2012_liftover_maps/TAIR10_lifted_chr4.txt",
        #        "salome2012_liftover_maps/TAIR10_lifted_chr5.txt",]
        #
        # Calculated mean recombination rates
        # raw_maps = [msprime.RateMap.read_hapmap(file, rate_col=2) for file in maps]
        # mean_rates = [m.mean_rate for m in raw_maps]
        # values pasted from the above output

        corrected_rate_means = {
            "1": 3.5149378684933275e-08 * rate_correction,
            "2": 3.8417724671111925e-08 * rate_correction,
            "3": 3.644268465988057e-08 * rate_correction,
            "4": 4.122114539275355e-08 * rate_correction,
            "5": 3.6486022007800115e-08 * rate_correction,
            "Pt": 0.0,
            "Mt": 0.0,
        }
        assert (
            corrected_rate_means[chr_id]
            == self.genome.get_chromosome(chr_id).recombination_rate
        )

    def test_selfing_correction_local(self):
        # test a local rate
        test_chrom = "3"
        genetic_map = "SalomeAveraged_TAIR10"
        local_rate = self.species.get_contig(
            test_chrom, genetic_map=genetic_map
        ).recombination_map.rate[2]
        # from above: raw_maps[2].rate[2] * rate_correction
        test_rate = 2.2387625242718484e-09
        assert local_rate == test_rate

    def test_assembly_source(self):
        assert self.genome.assembly_source == "ensembl"

    def test_assembly_build_version(self):
        assert self.genome.assembly_build_version is not None
