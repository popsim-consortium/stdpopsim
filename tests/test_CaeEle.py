import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("CaeEle")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "caenorhabditis_elegans"

    def test_name(self):
        assert self.species.name == "Caenorhabditis elegans"

    def test_common_name(self):
        assert self.species.common_name == "C. elegans"

    # QC Tests. These tests are performed by another contributor
    # independently referring to the citations provided in the
    # species definition, filling in the appropriate values
    # and deleting the pytest "skip" annotations.
    def test_qc_population_size(self):
        assert self.species.population_size == 10000

    def test_qc_generation_time(self):
        assert self.species.generation_time == 0.01


class TestGenomeData(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("CaeEle").genome
    mu = 1.84e-9

    # # downloading genetic map, for calculating total cM per chromosome
    # comm = "wget http://sesame.uoregon.edu/~ateterina/rockman2009_maps.tgz"
    # os.system(comm)
    # try:
    #     comm = "tar -zxvf rockman2009_maps.tgz"
    #     os.system(comm)
    # except:
    #     pass
    # my_recs = {}
    # genome = stdpopsim.get_species("CaeEle").genome
    # chromlist = ["I", "II", "III", "IV", "V", "X"]
    # for chrom in chromlist:
    #     bps, rates, cM = [], [], 0
    #     with open(
    #         "genetic_map/C.elegans.Rockman.Kruglyak.2009." + chrom + ".hapmap.txt"
    #     ) as infile:
    #         infile.readline()  # header
    #         newline = infile.readline().strip().split()
    #         bp = float(newline[1])
    #         rate = float(newline[2])
    #         bps.append(bp)
    #         rates.append(rate)
    #         for line in infile:
    #             newline = line.strip().split()
    #             bp = float(newline[1])
    #             rate = float(newline[2])
    #             bps.append(bp)
    #             rates.append(rate)
    #             previous_index = bps.index(bp) - 1
    #             length = bp - bps[previous_index]
    #             cM += rates[previous_index] * length
    #     total_rate = cM / genome.get_chromosome(chrom).length
    #     total_rate /= 100  # cM to rate
    #     total_rate /= 1000000  # Mb to bp
    #     total_rate /= 1000  # 0.1% outcrossing
    #     my_recs[chrom] = total_rate
    #     print(my_recs)

    my_recs = {
        "I": 3.1216265402124167e-11,
        "II": 3.529290802315087e-11,
        "III": 3.906598767640363e-11,
        "IV": 2.712098077556377e-11,
        "V": 2.4705737572511805e-11,
        "X": 2.9472374817864404e-11,
        "MtDNA": 0,
    }

    @pytest.mark.parametrize(
        ["name", "rate"],
        my_recs.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(
            self.genome.get_chromosome(name).recombination_rate
        )

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "I": mu,
            "II": mu,
            "III": mu,
            "IV": mu,
            "V": mu,
            "X": mu,
            "MtDNA": 1.05e-7,
        }.items(),
    )
    def test_mutation_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).mutation_rate)
