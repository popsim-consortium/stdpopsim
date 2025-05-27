import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("SusScr")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "sus_scrofa"

    def test_name(self):
        assert self.species.name == "Sus scrofa"

    def test_common_name(self):
        assert self.species.common_name == "Pig"

    def test_assembly_source(self):
        assert self.species.genome.assembly_source == "ensembl"

    def test_assembly_build_version(self):
        assert self.species.genome.assembly_build_version == "113"

    # Effective population size (Zhang et al., 2022: page 1048)
    # https://doi.org/10.1016/j.gpb.2022.02.001
    def test_qc_population_size(self):
        assert self.species.population_size == 270_000

    # Generation time/interal (Servanty et al., 2011: page 835)
    # https://doi.org/10.1111/j.1365-2664.2011.02017.x
    def test_qc_generation_time(self):
        assert self.species.generation_time == 3.6


class TestGenomeData(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("SusScr").genome

    # Pedigree-based recombination rate (Johnsson et al., 2021)
    # https://doi.org/10.1186/s12711-021-00643-0
    # 1. Downloaded Additional file 4: Table S3 sex-averaged rec rate
    #    https://static-content.springer.com/esm/art%3A10.1186%2Fs12711-021-00643-0/MediaObjects/12711_2021_643_MOESM4_ESM.csv
    # 2. Summed rec rate column for each chromosome to get length in Morgans
    # 3. Divided the sum by reported chromosome length in the publication
    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            # Morgans / n_base_pairs
            "1": 5.335247609316375e-09,  # 1.46339366 / 274287862
            "2": 8.605726100366067e-09,  # 1.307359715 / 151917421
            "3": 9.751141305671909e-09,  # 1.293980615 / 132700427
            "4": 9.794600488262457e-09,  # 1.282122117 / 130900910
            "5": 1.199263229672395e-08,  # 1.252961512 / 104477606
            "6": 8.786173522131294e-09,  # 1.500903759 / 170825645
            "7": 1.0887772281713889e-08,  # 1.325940215 / 121782508
            "8": 8.923031415618725e-09,  # 1.239683313 / 138930735
            "9": 9.337506177380429e-09,  # 1.301664487 / 139401727
            "10": 1.6542332055680058e-08,  # 1.146706799 / 69319537
            "11": 1.2145120578614156e-08,  # 0.960345302 / 79072521
            "12": 1.690506246134428e-08,  # 1.041322399 / 61598258
            "13": 6.162796630989997e-09,  # 1.283345448 / 208240759
            "14": 8.721569373637597e-09,  # 1.23601441 / 141719266
            "15": 8.114468763191382e-09,  # 1.139305203 / 140404164
            "16": 1.1407464686298615e-08,  # 0.906317069 / 79449474
            "17": 1.3779949520738576e-08,  # 0.874362601 / 63451800
            "18": 1.3679923895960038e-08,  # 0.763439439 / 55807287
            # using average acrosss autsomes
            "X": 9.415396334579284e-09,  # 21.31916806 / 2264287907
            "Y": 0.0,  # no recombination
            "MT": 0.0,  # no recombination
        }.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(
            self.genome.get_chromosome(name).recombination_rate
        )

    # De-novo mutation rate (Zhang et al., 2022: page 1042)
    # https://doi.org/10.1016/j.gpb.2022.02.001
    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": 3.6e-9,
            "2": 3.6e-9,
            "3": 3.6e-9,
            "4": 3.6e-9,
            "5": 3.6e-9,
            "6": 3.6e-9,
            "7": 3.6e-9,
            "8": 3.6e-9,
            "9": 3.6e-9,
            "10": 3.6e-9,
            "11": 3.6e-9,
            "12": 3.6e-9,
            "13": 3.6e-9,
            "14": 3.6e-9,
            "15": 3.6e-9,
            "16": 3.6e-9,
            "17": 3.6e-9,
            "18": 3.6e-9,
            "X": 3.6e-9,
            "Y": 3.6e-9,
            "MT": 3.6e-9,
        }.items(),
    )
    def test_mutation_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).mutation_rate)

    # Ploidy - pig is diploid with XY sex chromosome system
    @pytest.mark.parametrize(
        ["name", "ploidy"],
        {
            "1": 2,
            "2": 2,
            "3": 2,
            "4": 2,
            "5": 2,
            "6": 2,
            "7": 2,
            "8": 2,
            "9": 2,
            "10": 2,
            "11": 2,
            "12": 2,
            "13": 2,
            "14": 2,
            "15": 2,
            "16": 2,
            "17": 2,
            "18": 2,
            "X": 2,
            "Y": 1,
            "MT": 1,
        }.items(),
    )
    def test_chromosome_ploidy(self, name, ploidy):
        assert ploidy == self.genome.get_chromosome(name).ploidy
