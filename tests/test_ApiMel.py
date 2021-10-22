import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("ApiMel")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "apis_mellifera"

    def test_name(self):
        assert self.species.name == "Apis mellifera"

    def test_common_name(self):
        assert self.species.common_name == "Apis mellifera (DH4)"

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

    genome = stdpopsim.get_species("ApiMel").genome

    @pytest.mark.skip("Recombination rate QC not done yet")
    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "CM009931.2": -1,
            "CM009932.2": -1,
            "CM009933.2": -1,
            "CM009934.2": -1,
            "CM009935.2": -1,
            "CM009936.2": -1,
            "CM009937.2": -1,
            "CM009938.2": -1,
            "CM009939.2": -1,
            "CM009940.2": -1,
            "CM009941.2": -1,
            "CM009942.2": -1,
            "CM009943.2": -1,
            "CM009944.2": -1,
            "CM009945.2": -1,
            "CM009946.2": -1,
            "CM009947.2": -1,
        }.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert pytest.approx(rate, self.genome.get_chromosome(name).recombination_rate)

    @pytest.mark.skip("Mutation rate QC not done yet")
    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "CM009931.2": -1,
            "CM009932.2": -1,
            "CM009933.2": -1,
            "CM009934.2": -1,
            "CM009935.2": -1,
            "CM009936.2": -1,
            "CM009937.2": -1,
            "CM009938.2": -1,
            "CM009939.2": -1,
            "CM009940.2": -1,
            "CM009941.2": -1,
            "CM009942.2": -1,
            "CM009943.2": -1,
            "CM009944.2": -1,
            "CM009945.2": -1,
            "CM009946.2": -1,
            "CM009947.2": -1,
        }.items(),
    )
    def test_mutation_rate(self, name, rate):
        assert pytest.approx(rate, self.genome.get_chromosome(name).mutation_rate)
