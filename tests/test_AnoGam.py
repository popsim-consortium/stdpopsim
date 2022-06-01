import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("AnoGam")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "anopheles_gambiae"

    def test_name(self):
        assert self.species.name == "Anopheles gambiae"

    def test_common_name(self):
        assert self.species.common_name == "Anopheles gambiae"

    def test_qc_population_size(self):
        # based on theta = 4 Ne u since diversity in most pops is 1.5%
        assert self.species.population_size == pytest.approx(
            0.015 / (4 * 3.5e-9), rel=0.1
        )

    def test_qc_generation_time(self):
        assert self.species.generation_time == 1 / 11


class TestGenomeData(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("AnoGam").genome
    mut_rate = 3.5e-9

    # these values are given in the text of Pombi et al
    # except no single value is given for the X; the number
    # given is the average across chromosome from Pombi's Table 1
    # No value is given for 2L, because of the inversion,
    # so we set it equal to 2R.
    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "2L": 1.3e-8,
            "2R": 1.3e-8,
            "3L": 1.3e-8,
            "3R": 1.6e-8,
            "X": round(41.379 / (22.104923 - 1.825909) * 1e-8, 10),
            "Mt": 0,
        }.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(
            self.genome.get_chromosome(name).recombination_rate
        )

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "2L": mut_rate,
            "2R": mut_rate,
            "3L": mut_rate,
            "3R": mut_rate,
            "X": mut_rate,
            "Mt": mut_rate,
        }.items(),
    )
    def test_mutation_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).mutation_rate)
