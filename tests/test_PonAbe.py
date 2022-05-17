import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("PonAbe")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "pongo_abelii"

    def test_name(self):
        assert self.species.name == "Pongo abelii"

    def test_common_name(self):
        assert self.species.common_name == "Sumatran orangutan"

    def test_qc_population_size(self):
        assert self.species.population_size == 1.79 * 10**4

    def test_qc_generation_time(self):
        assert self.species.generation_time == 25


class TestGenomeData(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("PonAbe").genome

    # Average rates obtained from the NaterPA_PonAbe3 maps.
    # Loading all maps takes ~60 seconds, so we just hardcode the values.
    # This is how they were obtained:
    # # maps = {
    # #     chrom.id: gmap.get_chromosome_map(chrom.id)
    # #     for chrom in species.genome.chromosomes
    # #     if chrom.id not in ("X", "MT")
    # # }
    # # mean_rates = {k: v.mean_rate for k, v in maps.items()}
    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": 5.7209738520466636e-09,
            "2A": 6.746942995474612e-09,
            "2B": 6.088330307532069e-09,
            "3": 5.776152240522997e-09,
            "4": 6.061756714136327e-09,
            "5": 5.8783946576216106e-09,
            "6": 5.896045046261895e-09,
            "7": 6.82979118342936e-09,
            "8": 6.343129574251662e-09,
            "9": 7.047556965001838e-09,
            "10": 6.948119058906713e-09,
            "11": 5.791845481175716e-09,
            "12": 6.043159266363315e-09,
            "13": 6.4081697383466255e-09,
            "14": 6.110753350074597e-09,
            "15": 6.631500663980084e-09,
            "16": 6.7606578552415915e-09,
            "17": 8.01982524161532e-09,
            "18": 6.400648881415157e-09,
            "19": 9.096468143371218e-09,
            "20": 6.314242320359848e-09,
            "21": 7.588672333415156e-09,
            "22": 9.836284925758108e-09,
            # X isn't in the map, so this is the weighted genome-wide average.
            # # np.average(
            # #    [m.mean_rate for m in maps.values()],
            # #    weights=[m.sequence_length for m in maps.values()],
            # # )
            "X": 6.402884398829772e-09,
            "MT": 0,
        }.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(
            self.genome.get_chromosome(name).recombination_rate
        )

    _genome_mutation_rate = 1.5e-8

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": _genome_mutation_rate,
            "2a": _genome_mutation_rate,
            "2b": _genome_mutation_rate,
            "3": _genome_mutation_rate,
            "4": _genome_mutation_rate,
            "5": _genome_mutation_rate,
            "6": _genome_mutation_rate,
            "7": _genome_mutation_rate,
            "8": _genome_mutation_rate,
            "9": _genome_mutation_rate,
            "10": _genome_mutation_rate,
            "11": _genome_mutation_rate,
            "12": _genome_mutation_rate,
            "13": _genome_mutation_rate,
            "14": _genome_mutation_rate,
            "15": _genome_mutation_rate,
            "16": _genome_mutation_rate,
            "17": _genome_mutation_rate,
            "18": _genome_mutation_rate,
            "19": _genome_mutation_rate,
            "20": _genome_mutation_rate,
            "21": _genome_mutation_rate,
            "22": _genome_mutation_rate,
            "X": _genome_mutation_rate,
            "MT": _genome_mutation_rate,  # XXX: this is almost certainly wrong!
        }.items(),
    )
    def test_mutation_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).mutation_rate)
