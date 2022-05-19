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

    def test_qc_population_size(self):
        assert self.species.population_size == 2e05

    def test_qc_generation_time(self):
        assert self.species.generation_time == 2


class TestGenomeData(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("ApiMel").genome

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "CM009931.2": 23.9e-8,
            "CM009932.2": 24.6e-8,
            "CM009933.2": 24.1e-8,
            "CM009934.2": 27.6e-8,
            "CM009935.2": 21.4e-8,
            "CM009936.2": 21.2e-8,
            "CM009937.2": 23.4e-8,
            "CM009938.2": 20.9e-8,
            "CM009939.2": 24.6e-8,
            "CM009940.2": 24.8e-8,
            "CM009941.2": 20.3e-8,
            "CM009942.2": 21.2e-8,
            "CM009943.2": 23.4e-8,
            "CM009944.2": 24.6e-8,
            "CM009945.2": 22.1e-8,
            "CM009946.2": 22.8e-8,
            "CM009947.2": 0,
        }.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(
            self.genome.get_chromosome(name).recombination_rate
        )

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "CM009931.2": 3.4e-9,
            "CM009932.2": 3.4e-9,
            "CM009933.2": 3.4e-9,
            "CM009934.2": 3.4e-9,
            "CM009935.2": 3.4e-9,
            "CM009936.2": 3.4e-9,
            "CM009937.2": 3.4e-9,
            "CM009938.2": 3.4e-9,
            "CM009939.2": 3.4e-9,
            "CM009940.2": 3.4e-9,
            "CM009941.2": 3.4e-9,
            "CM009942.2": 3.4e-9,
            "CM009943.2": 3.4e-9,
            "CM009944.2": 3.4e-9,
            "CM009945.2": 3.4e-9,
            "CM009946.2": 3.4e-9,
            "CM009947.2": 3.4e-9,
        }.items(),
    )
    def test_mutation_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).mutation_rate)
