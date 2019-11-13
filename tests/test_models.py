"""
Tests for simulation model infrastructure.
"""
import unittest
import itertools
import io

import numpy as np
import msprime

import stdpopsim
from stdpopsim import models


class ModelTestMixin(object):
    """
    Mixin for testing specific models. Subclasses should extend
    unittest.TestCase and this mixin, and define the self.model (as the
    model instance).
    """
    # To be defined in subclasses.
    model = None

    def test_debug_runs(self):
        output = io.StringIO()
        self.model.debug(output)
        s = output.getvalue()
        self.assertGreater(len(s), 0)

    def test_simulation_runs(self):
        # With a recombination_map of None, we simulate a coalescent without
        # recombination in msprime, with no mutation.
        contig = stdpopsim.Contig()
        samples = self.model.get_samples(*([2] * self.model.num_populations))
        ts = self.model.simulate(contig, samples)
        self.assertEqual(ts.num_populations, self.model.num_populations)


class QcdModelTestMixin(ModelTestMixin):
    """
    Extends the tests to also check that the qc model is equal to
    the production model.
    """
    # To be defined in subclass.
    qc_model = None

    def test_qc_model_equal(self):
        self.assertTrue(self.model.equals(self.qc_model))


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
        with self.assertRaises(models.UnequalModelsError):
            models.verify_population_configurations_equal([pc], [pc, pc])

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
            with self.assertRaises(models.UnequalModelsError):
                models.verify_population_configurations_equal(pc_list2, pc_list1)

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
            with self.assertRaises(models.UnequalModelsError):
                models.verify_population_configurations_equal(pc_list2, pc_list1)


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
        with self.assertRaises(models.UnequalModelsError):
            models.verify_demographic_events_equal([], events, 1)

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
            with self.assertRaises(models.UnequalModelsError):
                models.verify_demographic_events_equal(e1[:j], e2[:j], 1)
            with self.assertRaises(models.UnequalModelsError):
                models.verify_demographic_events_equal(e2[:j], e1[:j], 1)

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
            with self.assertRaises(models.UnequalModelsError):
                models.verify_demographic_events_equal([b], [a], 1)
            with self.assertRaises(models.UnequalModelsError):
                models.verify_demographic_events_equal([a], [b], 1)

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
            with self.assertRaises(models.UnequalModelsError):
                models.verify_demographic_events_equal([b], [a], 1)
            with self.assertRaises(models.UnequalModelsError):
                models.verify_demographic_events_equal([a], [b], 1)

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
            with self.assertRaises(models.UnequalModelsError):
                models.verify_demographic_events_equal([b], [a], 1)
            with self.assertRaises(models.UnequalModelsError):
                models.verify_demographic_events_equal([a], [b], 1)

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
            with self.assertRaises(models.UnequalModelsError):
                models.verify_demographic_events_equal([b], [a], 1)
            with self.assertRaises(models.UnequalModelsError):
                models.verify_demographic_events_equal([a], [b], 1)

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
            with self.assertRaises(models.UnequalModelsError):
                models.verify_demographic_events_equal([b], [a], 1)
            with self.assertRaises(models.UnequalModelsError):
                models.verify_demographic_events_equal([a], [b], 1)


class DummyModel(models.Model):
    """
    Dummy subclass to make sure we're filtering models correctly.
    """


class TestAllModels(unittest.TestCase):
    """
    Tests that we can get all known simulation models.
    """
    def test_non_empty(self):
        self.assertGreater(len(list(stdpopsim.all_models())), 0)

    def test_all_instances(self):
        for model in stdpopsim.all_models():
            self.assertIsInstance(model, models.Model)

    def test_filtering_outside_classes(self):
        for model in stdpopsim.all_models():
            self.assertNotIsInstance(model, DummyModel)

    def test_generation_times_non_empty(self):
        self.assertGreater(len([model.generation_time for model in
                                stdpopsim.all_models()]), 0)


class TestModelsEqual(unittest.TestCase):
    """
    Tests Model object equality comparison.
    """
    def test_known_models(self):
        # All models should be equal to themselves.
        other_model = models.Model()
        for model in stdpopsim.all_models():
            self.assertTrue(model.equals(model))
            self.assertFalse(model.equals(other_model))

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


class TestModelProperties(unittest.TestCase):
    def test_model_generation_time(self):
        self.assertTrue(models.Model().generation_time == -1)
        known_models = list(stdpopsim.all_models())
        n = len(known_models)
        for j in range(n):
            self.assertTrue(known_models[j].generation_time > -2)


class TestConstantSizeModel(unittest.TestCase, ModelTestMixin):
    model = models.PiecewiseConstantSize(100)


class TestTwoEpochModel(unittest.TestCase, ModelTestMixin):
    model = models.PiecewiseConstantSize(100, (10, 10))


class TestPiecewiseConstantSize(unittest.TestCase):
    """
    Specific tests for the piecewise constant model.
    """
    def test_single_epoch(self):
        model = models.PiecewiseConstantSize(100)
        self.assertEqual(model.population_configurations[0].initial_size, 100)
        self.assertEqual(len(model.population_configurations), 1)
        self.assertEqual(len(model.demographic_events), 0)

    def test_two_epoch(self):
        model = models.PiecewiseConstantSize(50, (10, 100))
        self.assertEqual(model.population_configurations[0].initial_size, 50)
        self.assertEqual(len(model.demographic_events), 1)
        event = model.demographic_events[0]
        self.assertEqual(event.time, 10)
        self.assertEqual(event.initial_size, 100)
        self.assertEqual(event.growth_rate, 0)

    def test_three_epoch(self):
        model = models.PiecewiseConstantSize(0.1, (0.1, 10), (0.2, 100))
        self.assertEqual(model.population_configurations[0].initial_size, 0.1)
        self.assertEqual(len(model.demographic_events), 2)
        event = model.demographic_events[0]
        self.assertEqual(event.time, 0.1)
        self.assertEqual(event.initial_size, 10)
        self.assertEqual(event.growth_rate, 0)
        event = model.demographic_events[1]
        self.assertEqual(event.time, 0.2)
        self.assertEqual(event.initial_size, 100)
        self.assertEqual(event.growth_rate, 0)


class TestGenericIM(unittest.TestCase):
    """
    Tests for the generic IM model.
    """
    def test_pop_configs(self):
        model = models.GenericIM(100, 200, 300, 50, 0, 0)
        self.assertEqual(len(model.population_configurations), 3)
        self.assertEqual(model.population_configurations[0].initial_size, 100)
        self.assertEqual(model.population_configurations[1].initial_size, 200)
        self.assertEqual(model.population_configurations[2].initial_size, 300)

    def test_split(self):
        model = models.GenericIM(100, 200, 300, 50, 0, 0)
        self.assertEqual(len(model.demographic_events), 2)
        self.assertEqual(model.demographic_events[0].time, 50)
        self.assertEqual(model.demographic_events[1].time, 50)

    def test_migration_rates(self):
        model = models.GenericIM(100, 200, 300, 50, 0.002, 0.003)
        self.assertEqual(np.shape(model.migration_matrix), (3, 3))
        self.assertEqual(model.migration_matrix[0][1], 0.002)
        self.assertEqual(model.migration_matrix[1][0], 0.003)
        # check these are the only two nonzero entries in the migration matrix
        self.assertEqual(np.sum(np.array(model.migration_matrix) != 0), 2)
