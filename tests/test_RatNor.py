import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("RatNor")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "rattus_norvegicus"

    def test_name(self):
        assert self.species.name == "Rattus norvegicus"

    def test_common_name(self):
        assert self.species.common_name == "Rat"

    # QC Tests. These tests are performed by another contributor
    # independently referring to the citations provided in the
    # species definition, filling in the appropriate values
    # and deleting the pytest "skip" annotations.
    @pytest.mark.skip("Population size QC not done yet")
    def test_qc_population_size(self):
        assert self.species.population_size == 124000

    @pytest.mark.skip("Generation time QC not done yet")
    def test_qc_generation_time(self):
        assert self.species.generation_time == 0.5


class TestGenomeData(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("RatNor").genome

    @pytest.mark.skip("Recombination rate QC not done yet")
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(
            self.genome.get_chromosome(name).recombination_rate
        )

    @pytest.mark.skip("Mutation rate QC not done yet")
    def test_mutation_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).mutation_rate)

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
            "19": 2,
            "20": 2,
            "X": 2,
            "Y": 1,
            "MT": 1,
        }.items(),
    )
    def test_chromosome_ploidy(self, name, ploidy):
        assert ploidy == self.genome.get_chromosome(name).ploidy
