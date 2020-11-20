"""
Tests for classes that hold information about genomic region to be simulated.
"""
import unittest
import stdpopsim
import numpy as np


class TestGenomicElementType(unittest.TestCase):
    def test_bad_proportions(self):
        # no proportion but mut type, prop > 1, sum(prop) > 1
        proportions = ([], [2], [0.5, 0.5, 0.5])
        mut_type_ids = ([0], [0], [0])
        for props in proportions:
            with self.assertRaises(ValueError):
                stdpopsim.GenomicElementType(
                    intervals=np.array([[0, 1]]),
                    mutation_type_ids=mut_type_ids,
                    proportions=props,
                )


class TestContig(unittest.TestCase):
    def get_test_contig(self, no_rec_map=False):
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr22", length_multiplier=0.001)
        if no_rec_map:
            contig.recombination_map = None
        return contig

    def test_length(self):
        contig = self.get_test_contig()
        self.assertTrue(contig.length is not None and contig.length > 0)
        contig_nomap = self.get_test_contig(no_rec_map=True)
        self.assertTrue(contig_nomap.length is None)

    def test_slim_fractions(self):
        contig = self.get_test_contig()
        self.assertTrue(contig.slim_fractions.size == 0)
        contig.fully_neutral()
        self.assertTrue(np.isclose(contig.slim_fractions, np.array([0])))
        contig.fully_neutral(slim_mutations=True)
        self.assertTrue(np.isclose(contig.slim_fractions, np.array([1])))

        proportions = ([0, 1.0], [0.5, 0.5], [0, 0], [0.2, 0])
        for props in proportions:
            slim_frac = np.array(sum(props))
            contig.clear_genomic_mutation_types()
            contig.add_genomic_element_type(
                intervals=np.array([[0, 1]]),
                mutation_types=[
                    stdpopsim.ext.MutationType(),
                    stdpopsim.ext.MutationType(),
                ],
                proportions=props,
            )
            self.assertTrue((slim_frac == contig.slim_fractions).all())

    def test_genomic_intervals(self):
        contig = self.get_test_contig()
        contig.fully_neutral()
        self.assertTrue(len(contig.genomic_intervals) == 1)

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
        self.assertTrue((contig.genomic_intervals == truth).all())

    def test_msp_mutation_rate_map(self):
        contig = self.get_test_contig()
        contig.fully_neutral()
        breaks, rates = contig.msp_mutation_rate_map
        self.assertTrue(breaks == [0, int(contig.length)])
        self.assertTrue(rates == [contig.mutation_rate])

    def test_slim_mutation_rate_map(self):
        contig = self.get_test_contig()
        contig.fully_neutral()
        breaks, rates = contig.slim_mutation_rate_map
        self.assertTrue(breaks == [int(contig.length) - 1])
        self.assertTrue(rates == [0.0])

    def test_complex_mutation_rate_maps(self):
        contig = self.get_test_contig()
        for prop1, prop2 in ((0.3, 0.1), (0.7, 0.2), (1.0, 0.0)):
            contig.clear_genomic_mutation_types()
            contig.add_genomic_element_type(
                intervals=np.array([[10, 100]]),
                mutation_types=[stdpopsim.ext.MutationType()],
                proportions=[prop1],
            )

            obs_msp_breaks, obs_msp_rates = contig.msp_mutation_rate_map
            obs_slim_breaks, obs_slim_rates = contig.slim_mutation_rate_map

            exp_msp_breaks = [0, 10, 100, int(contig.length)]
            exp_msp_rates = [
                contig.mutation_rate,
                contig.mutation_rate * (1 - prop1),
                contig.mutation_rate,
            ]
            self.assertEqual(obs_msp_breaks, exp_msp_breaks)
            self.assertTrue(
                (np.isclose(np.array(obs_msp_rates), np.array(exp_msp_rates))).all()
            )

            exp_slim_breaks = [9, 99, int(contig.length) - 1]
            exp_slim_rates = [0, contig.mutation_rate * prop1, 0]
            self.assertEqual(obs_slim_breaks, exp_slim_breaks)
            self.assertTrue(
                (np.isclose(np.array(obs_slim_rates), np.array(exp_slim_rates))).all()
            )

            contig.clear_genomic_mutation_types()
            contig.add_genomic_element_type(
                intervals=np.array([[0, 50]]),
                mutation_types=[stdpopsim.ext.MutationType()],
                proportions=[prop1],
            )
            contig.add_genomic_element_type(
                intervals=np.array([[50, 100]]),
                mutation_types=[stdpopsim.ext.MutationType()],
                proportions=[prop2],
            )

            obs_msp_breaks, obs_msp_rates = contig.msp_mutation_rate_map
            obs_slim_breaks, obs_slim_rates = contig.slim_mutation_rate_map

            exp_msp_breaks = [0, 50, 100, int(contig.length)]
            exp_msp_rates = [
                contig.mutation_rate * (1 - prop1),
                contig.mutation_rate * (1 - prop2),
                contig.mutation_rate,
            ]
            self.assertEqual(obs_msp_breaks, exp_msp_breaks)
            self.assertTrue(
                (np.isclose(np.array(obs_msp_rates), np.array(exp_msp_rates))).all()
            )

            exp_slim_breaks = [49, 99, int(contig.length) - 1]
            exp_slim_rates = [
                contig.mutation_rate * prop1,
                contig.mutation_rate * prop2,
                0,
            ]
            self.assertEqual(obs_slim_breaks, exp_slim_breaks)
            self.assertTrue(
                (np.isclose(np.array(obs_slim_rates), np.array(exp_slim_rates))).all()
            )

    def test_add_mutation_type_errors(self):
        contig = self.get_test_contig()
        # diff length mutation types and proportions
        with self.assertRaises(ValueError):
            contig.add_mutation_types([stdpopsim.ext.MutationType()], [], 0)
        # invalid ge type id
        with self.assertRaises(ValueError):
            contig.add_mutation_types([stdpopsim.ext.MutationType()], [0], 0)
        contig.fully_neutral()
        # sum props > 1
        with self.assertRaises(ValueError):
            contig.add_mutation_types([stdpopsim.ext.MutationType()], [1.1], 0)

    def test_add_genomic_element_type_errors(self):
        contig = self.get_test_contig()
        contig.recombination_map = None
        # can't add ge type without rec map
        with self.assertRaises(ValueError):
            contig.add_genomic_element_type(np.array([]), [], [])
        contig = self.get_test_contig()
        # bad intervals
        with self.assertRaises(ValueError):
            contig.add_genomic_element_type(np.array([10, 20]), [], [])

    def test_add_genomic_element_type(self):
        contig = self.get_test_contig()
        self.assertTrue(contig.mutation_types == [])
        self.assertTrue(contig.genomic_element_types == [])
        intervals1 = np.array([[10, 20], [20, 40], [100, 200]])
        intervals2 = np.array([[5, 10], [45, 50], [200, 250]])
        m1 = stdpopsim.ext.MutationType()
        m2 = stdpopsim.ext.MutationType()
        contig.add_genomic_element_type(intervals1, [m1, m2], [0.5, 0.3])
        contig.add_genomic_element_type(intervals2, [m1], [0.2])
        self.assertTrue(len(contig.genomic_element_types) == 2)
        self.assertTrue(len(contig.mutation_types) == 2)
