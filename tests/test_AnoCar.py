import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("AnoCar")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "anolis_carolinensis"

    def test_name(self):
        assert self.species.name == "Anolis carolinensis"

    def test_common_name(self):
        assert self.species.common_name == "Anole lizard"

    # QC Tests. These tests are performed by another contributor
    # independently referring to the citations provided in the
    # species definition, filling in the appropriate values
    # and deleting the pytest "skip" annotations.
    @pytest.mark.skip("Population size QC not done yet")
    def test_qc_population_size(self):
        assert self.species.population_size == -1

    @pytest.mark.skip("Generation time QC not done yet")
    def test_qc_generation_time(self):
        assert self.species.generation_time == -1


class TestGenomeData(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("AnoCar").genome

    @pytest.mark.skip("Recombination rate QC not done yet")
    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": -1,
            "2": -1,
            "3": -1,
            "4": -1,
            "5": -1,
            "6": -1,
            "LGa": -1,
            "LGb": -1,
            "LGc": -1,
            "LGd": -1,
            "LGf": -1,
            "LGg": -1,
            "LGh": -1,
            "MT": -1,
        }.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert pytest.approx(rate, self.genome.get_chromosome(name).recombination_rate)

    @pytest.mark.skip("Mutation rate QC not done yet")
    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": -1,
            "2": -1,
            "3": -1,
            "4": -1,
            "5": -1,
            "6": -1,
            "LGa": -1,
            "LGb": -1,
            "LGc": -1,
            "LGd": -1,
            "LGf": -1,
            "LGg": -1,
            "LGh": -1,
            "MT": -1,
        }.items(),
    )
    def test_mutation_rate(self, name, rate):
        assert pytest.approx(rate, self.genome.get_chromosome(name).mutation_rate)
