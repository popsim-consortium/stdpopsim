"""
Tests for classes that hold information about genomic region to be simulated.
"""
import numpy as np
import pytest

import stdpopsim


class TestGenomicElementType:
    def test_bad_proportions(self):
        # no proportion but mut type, prop > 1, sum(prop) > 1
        proportions = ([], [2], [0.5, 0.5, 0.5])
        mut_type_ids = ([0], [0], [0])
        for props in proportions:
            with pytest.raises(ValueError):
                stdpopsim.GenomicElementType(
                    intervals=np.array([[0, 1]]),
                    mutation_type_ids=mut_type_ids,
                    proportions=props,
                )


class TestContig:
    def get_test_contig(self):
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr22", length_multiplier=0.001)
        return contig

    def test_all_intervals_array(self):
        contig = self.get_test_contig()
        contig.fully_neutral()
        assert len(contig.all_intervals_array) == 1

        truth = np.array([[20, 30, 1], [30, 50, 0]])
        contig.clear_genomic_mutation_types()
        contig.add_genomic_element_type(
            intervals=np.array([[30, 50, 0]]),
            mutation_types=[stdpopsim.ext.MutationType()],
            proportions=np.array([0]),
        )
        contig.add_genomic_element_type(
            intervals=np.array([[20, 30, 1]]),
            mutation_types=[stdpopsim.ext.MutationType()],
            proportions=np.array([0]),
        )
        assert (contig.all_intervals_array == truth).all()

    def test_fully_neutral_error(self):
        contig = self.get_test_contig()
        contig.recombination_map = None
        with pytest.raises(ValueError):
            contig.fully_neutral()

    def test_add_mutation_type_errors(self):
        contig = self.get_test_contig()
        # diff length mutation types and proportions
        with pytest.raises(ValueError):
            contig.add_mutation_types([stdpopsim.ext.MutationType()], [], 0)
        # invalid ge type id
        with pytest.raises(ValueError):
            contig.add_mutation_types([stdpopsim.ext.MutationType()], [0], 1)
        contig.fully_neutral()
        # sum props > 1
        with pytest.raises(ValueError):
            contig.add_mutation_types([stdpopsim.ext.MutationType()], [1.1], 0)

    def test_add_genomic_element_type_errors(self):
        contig = self.get_test_contig()
        contig.recombination_map = None
        # can't add ge type without rec map
        with pytest.raises(ValueError):
            contig.add_genomic_element_type(np.array([]), [], [])
        contig = self.get_test_contig()
        # bad intervals
        with pytest.raises(ValueError):
            contig.add_genomic_element_type(np.array([10, 20]), [], [])

    def test_add_genomic_element_type(self):
        contig = self.get_test_contig()
        assert contig.mutation_types == []
        assert contig.genomic_element_types == []
        intervals1 = np.array([[10, 20], [20, 40], [100, 200]])
        intervals2 = np.array([[5, 10], [45, 50], [200, 250]])
        m1 = stdpopsim.ext.MutationType()
        m2 = stdpopsim.ext.MutationType()
        contig.add_genomic_element_type(intervals1, [m1, m2], [0.5, 0.3])
        contig.add_genomic_element_type(intervals2, [m1], [0.2])
        assert len(contig.genomic_element_types) == 2
        assert len(contig.mutation_types) == 2
