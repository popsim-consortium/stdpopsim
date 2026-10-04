import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("PruDul")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "prunus_dulcis"

    def test_name(self):
        assert self.species.name == "Prunus dulcis"

    def test_common_name(self):
        assert self.species.common_name == "Prunus dulcis"

    def test_assembly_source(self):
        assert self.species.genome.assembly_source == "manual"

    def test_assembly_build_version(self):
        assert self.species.genome.assembly_build_version == "3"

    # Effective population size (Velasco et al., 2016: Figure S8, page 3987-8)
    # https://doi.org/10.1534/g3.116.032672
    def test_qc_population_size(self):
        assert self.species.population_size == 216627

    # Generation time/interal (Velasco et al., 2016: Figure S8, page 3987)
    # https://doi.org/10.1534/g3.116.032672
    def test_qc_generation_time(self):
        assert self.species.generation_time == 10


class TestGenomeData(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("PruDul").genome

    # Recombination rates from linkage map in Mas-Gomez et al. (2025) Table 1
    # and physical lengths from genome assembly
    # in Castanera et al. (2024) - file genome_data.py
    # https://doi.org/10.1016/j.hpj.2025.04.013
    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            # Morgans / n_base_pairs
            "1": 1.97e-8,  # 99.41 / 50452767 = 0.00000197
            "2": 1.51e-8,  # 46.88 / 30972344 = 0.00000151
            "3": 1.89e-8,  # 57.52 / 30498384 = 0.00000189
            "4": 2.22e-8,  # 62.41 / 28068711 = 0.00000222
            "5": 2.54e-8,  # 56.54 / 22241271 = 0.00000254
            "6": 1.98e-8,  # 64.19 / 32371708 = 0.00000198
            "7": 2.61e-8,  # 64.98 / 24859559 = 0.00000261
            "8": 2.12e-8,  # 57.36 / 27108964 = 0.00000212
            "MT": 0,  # no recombination
            "Pt": 0,  # no recombination
        }.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(
            self.genome.get_chromosome(name).recombination_rate
        )

    # Mutation rate per site per generation
    # Velasco et al. (2016), page 3987
    # https://doi.org/10.1534/g3.116.032672
    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": 1e-8,
            "2": 1e-8,
            "3": 1e-8,
            "4": 1e-8,
            "5": 1e-8,
            "6": 1e-8,
            "7": 1e-8,
            "8": 1e-8,
            "MT": 1e-8,
            "Pt": 1e-8,
        }.items(),
    )
    def test_mutation_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).mutation_rate)

    # Ploidy - almond is diploid
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
            "MT": 1,
            "Pt": 1,
        }.items(),
    )
    def test_chromosome_ploidy(self, name, ploidy):
        assert ploidy == self.genome.get_chromosome(name).ploidy
