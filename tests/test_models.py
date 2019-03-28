"""
Tests for simulation model infrastructure.
"""
import unittest
import itertools

import numpy as np
import msprime

from stdpopsim import models
from stdpopsim import homo_sapiens
from stdpopsim import drosophila_melanogaster
from stdpopsim import pongo
from stdpopsim import e_coli
# FIXME getting 'cannot import name 'arabadopsis_thaliana' errors here
# from stdpopsim import arabadopsis_thaliana


class TestPopulationConfigsEqual(unittest.TestCase):
    """
    Tests the equality comparison of different msprime population configurations.
    """
    def test_empty(self):
        self.assertTrue(models.population_configurations_equal([], []))

    def test_sample_size_error(self):
        pc1 = msprime.PopulationConfiguration(sample_size=2, initial_size=1)
        pc2 = msprime.PopulationConfiguration(initial_size=1)
        with self.assertRaises(ValueError):
            models.population_configurations_equal([pc1], [pc1])
        with self.assertRaises(ValueError):
            models.population_configurations_equal([pc1], [pc2])
        with self.assertRaises(ValueError):
            models.population_configurations_equal([pc2], [pc1])

    def test_no_initial_size_error(self):
        pc1 = msprime.PopulationConfiguration()
        pc2 = msprime.PopulationConfiguration(initial_size=1)
        with self.assertRaises(ValueError):
            self.assertTrue(models.population_configurations_equal([pc1], [pc1]))
        with self.assertRaises(ValueError):
            self.assertTrue(models.population_configurations_equal([pc1], [pc2]))
        with self.assertRaises(ValueError):
            self.assertTrue(models.population_configurations_equal([pc2], [pc1]))

    def test_different_lengths(self):
        pc = msprime.PopulationConfiguration(initial_size=1)
        self.assertFalse(models.population_configurations_equal([], [pc]))
        self.assertFalse(models.population_configurations_equal([pc], []))
        self.assertFalse(models.population_configurations_equal([pc, pc], [pc]))
        self.assertFalse(models.population_configurations_equal([pc], [pc, pc]))
        self.assertTrue(models.population_configurations_equal([pc], [pc]))

    def test_initial_sizes(self):
        test_sizes = [
            ([1], [1.001]),
            ([1.001], [1]),
            ([1, 2, 3], [2, 3, 4]),
            (np.arange(1, 100), np.arange(1, 100) - 0.001),
        ]

        for sizes1, sizes2 in test_sizes:
            pc_list1 = [
                msprime.PopulationConfiguration(initial_size=size) for size in sizes1]
            pc_list2 = [
                msprime.PopulationConfiguration(initial_size=size) for size in sizes2]
            self.assertFalse(models.population_configurations_equal(pc_list1, pc_list2))
            self.assertFalse(models.population_configurations_equal(pc_list2, pc_list1))
            self.assertTrue(models.population_configurations_equal(pc_list1, pc_list1))
            self.assertTrue(models.population_configurations_equal(pc_list2, pc_list2))

    def test_growth_rates(self):
        test_rates = [
            ([1], [1.001]),
            ([-1.001], [-1]),
            ([1, 2, 3], [2, 3, 4]),
            (np.arange(1, 100), np.arange(1, 100) - 0.001),
        ]
        for rates1, rates2 in test_rates:
            pc_list1 = [
                msprime.PopulationConfiguration(
                    initial_size=1, growth_rate=rate) for rate in rates1]
            pc_list2 = [
                msprime.PopulationConfiguration(
                    initial_size=1, growth_rate=rate) for rate in rates2]
            self.assertFalse(models.population_configurations_equal(pc_list1, pc_list2))
            self.assertFalse(models.population_configurations_equal(pc_list2, pc_list1))
            self.assertTrue(models.population_configurations_equal(pc_list1, pc_list1))
            self.assertTrue(models.population_configurations_equal(pc_list2, pc_list2))


class TestDemographicEventsEqual(unittest.TestCase):
    """
    Tests the equality comparison of different msprime demographic events.
    """
    def test_empty(self):
        self.assertTrue(models.demographic_events_equal([], [], 1))

    def test_different_lengths(self):
        events = [msprime.PopulationParametersChange(time=1, initial_size=1)] * 10
        self.assertFalse(models.demographic_events_equal(events[:1], [], 1))
        self.assertFalse(models.demographic_events_equal([], events[:1], 1))
        self.assertFalse(models.demographic_events_equal(events, [], 1))
        self.assertFalse(models.demographic_events_equal([], events, 1))

    def test_different_times(self):
        n = 10
        e1 = [
            msprime.PopulationParametersChange(time=j, initial_size=1)
            for j in range(n)]
        e2 = [
            msprime.PopulationParametersChange(time=j + 1, initial_size=1)
            for j in range(n)]
        for j in range(1, n):
            self.assertFalse(models.demographic_events_equal(e1[:j], e2[:j], 1))
            self.assertFalse(models.demographic_events_equal(e2[:j], e1[:j], 1))

    def test_different_types(self):
        events = [
            msprime.PopulationParametersChange(time=1, initial_size=1),
            msprime.MigrationRateChange(time=1, rate=1),
            msprime.MassMigration(time=1, source=1),
            msprime.SimpleBottleneck(time=1)]
        for a, b in itertools.combinations(events, 2):
            self.assertFalse(models.demographic_events_equal([a], [b], 1))
            self.assertFalse(models.demographic_events_equal([b], [a], 1))
            self.assertTrue(models.demographic_events_equal([a], [a], 1))
            self.assertTrue(models.demographic_events_equal([b], [b], 1))

    def test_population_parameters_change(self):

        def f(time=1, initial_size=1, growth_rate=None, population=None):
            return msprime.PopulationParametersChange(
                time=time, initial_size=initial_size, growth_rate=growth_rate,
                population=population)

        test_events = [
            (f(time=1), f(time=2)),
            (f(initial_size=1), f(initial_size=2)),
            (f(growth_rate=1), f(growth_rate=2)),
            (f(population=1), f(population=2)),
        ]
        for a, b in test_events:
            self.assertFalse(models.demographic_events_equal([a], [b], 1))
            self.assertFalse(models.demographic_events_equal([b], [a], 1))
            self.assertTrue(models.demographic_events_equal([a], [a], 1))
            self.assertTrue(models.demographic_events_equal([b], [b], 1))

    def test_migration_rate_change(self):

        def f(time=1, rate=1, matrix_index=None):
            return msprime.MigrationRateChange(
                time=time, rate=rate, matrix_index=matrix_index)

        test_events = [
            (f(time=1), f(time=2)),
            (f(rate=1), f(rate=2)),
            (f(matrix_index=[0, 1]), f(matrix_index=[0, 2])),
            (f(matrix_index=np.array([0, 1])), f(matrix_index=[0, 2])),
        ]
        for a, b in test_events:
            self.assertFalse(models.demographic_events_equal([a], [b], 1))
            self.assertFalse(models.demographic_events_equal([b], [a], 1))
            self.assertTrue(models.demographic_events_equal([a], [a], 1))
            self.assertTrue(models.demographic_events_equal([b], [b], 1))

    def test_mass_migration(self):

        def f(time=1, source=1, dest=1, proportion=1):
            return msprime.MassMigration(
                time=time, source=source, dest=dest, proportion=proportion)

        test_events = [
            (f(time=1), f(time=2)),
            (f(source=1), f(source=2)),
            (f(dest=1), f(dest=2)),
            (f(proportion=1), f(proportion=0.2)),
        ]
        for a, b in test_events:
            self.assertFalse(models.demographic_events_equal([a], [b], 1))
            self.assertFalse(models.demographic_events_equal([b], [a], 1))
            self.assertTrue(models.demographic_events_equal([a], [a], 1))
            self.assertTrue(models.demographic_events_equal([b], [b], 1))

    def test_simple_bottleneck(self):

        def f(time=1, population=1, proportion=1):
            return msprime.SimpleBottleneck(
                time=time, population=population, proportion=proportion)

        test_events = [
            (f(time=1), f(time=2)),
            (f(population=1), f(population=2)),
            (f(proportion=1), f(proportion=0.2)),
        ]
        for a, b in test_events:
            self.assertFalse(models.demographic_events_equal([a], [b], 1))
            self.assertFalse(models.demographic_events_equal([b], [a], 1))
            self.assertTrue(models.demographic_events_equal([a], [a], 1))
            self.assertTrue(models.demographic_events_equal([b], [b], 1))


class TestModelsEqual(unittest.TestCase):
    """
    Tests Model object equality comparison.
    """
    def test_known_models(self):
        # TODO this should be available via a function at the top level,
        # in __init__.py
        known_models = [
            homo_sapiens.GutenkunstThreePopOutOfAfrica(),
            homo_sapiens.TennessenEuropean(),
            drosophila_melanogaster.SheehanSongThreeEpoch(),
            drosophila_melanogaster.LiStephanTwoPopulation(),
            # arabadopsis_thaliana.Durvasula2017MSMC(),
            pongo.LockeEtAlPongoIM(),
            e_coli.LapierreConstant()]
        n = len(known_models)
        for j in range(n):
            for k in range(n):
                self.assertEqual(j == k, known_models[j].equals(known_models[k]))

    def test_different_objects(self):
        m1 = models.Model()
        self.assertFalse(m1.equals(self))
        self.assertFalse(m1.equals({}))
        self.assertFalse(m1.equals(None))

    def test_default_models(self):
        m1 = models.Model()
        m2 = models.Model()
        self.assertTrue(m1.equals(m2))

    def test_migration_matrices(self):
        m1 = models.Model()
        m2 = models.Model()
        m1.migration_matrix = [[]]
        self.assertFalse(m1.equals(m2))
        m2.migration_matrix = [[]]
        self.assertTrue(m1.equals(m2))
        mm1 = np.arange(100, dtype=float).reshape(10, 10)
        m1.migration_matrix = mm1
        self.assertFalse(m1.equals(m2))
        m2.migration_matrix = mm1
        self.assertTrue(m1.equals(m2))
        # Perturb the matrix by a tiny bit and see if we're still equal
        mm2 = mm1.copy() + 1e-9
        m2.migration_matrix = mm2
        self.assertFalse(np.all(m1.migration_matrix == m2.migration_matrix))
        self.assertTrue(m1.equals(m2))
        # If we have higher tolerances we catch the differences
        self.assertFalse(m1.equals(m2, atol=1e-10, rtol=1e-9))
