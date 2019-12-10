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
        # Generate vector with 2 samples for each pop with sampling enabled
        sample_count = []
        for p in self.model.populations:
            if p.allow_samples:
                sample_count.append(2)
            else:
                sample_count.append(0)
        samples = self.model.get_samples(*sample_count)
        engine = stdpopsim.get_default_engine()
        ts = engine.simulate(self.model, contig, samples)
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


class TestAllModels(unittest.TestCase):
    """
    Tests for registered simulation models.
    """
    def test_non_empty(self):
        self.assertGreater(len(list(stdpopsim.all_demographic_models())), 0)

    def test_all_instances(self):
        for model in stdpopsim.all_demographic_models():
            self.assertIsInstance(model, models.DemographicModel)

            self.assertGreater(len(model.id), 0)
            self.assertGreater(len(model.name), 0)
            self.assertGreater(len(model.description), 0)
            self.assertGreater(len(model.citations), 0)
            self.assertGreater(model.generation_time, 0)

            npops = len(model.populations)
            self.assertGreater(npops, 0)
            self.assertEqual(len(model.population_configurations), npops)
            self.assertEqual(len(model.migration_matrix), npops)
            self.assertEqual(len(model.migration_matrix[0]), npops)
            self.assertIsInstance(model.demographic_events, list)


class TestModelsEqual(unittest.TestCase):
    """
    Tests DemographicModel object equality comparison.
    """
    def test_known_models(self):
        # All models should be equal to themselves.
        other_model = models.DemographicModel.empty()
        for model in stdpopsim.all_demographic_models():
            self.assertTrue(model.equals(model))
            self.assertFalse(model.equals(other_model))

    def test_different_objects(self):
        m1 = models.DemographicModel.empty()
        self.assertFalse(m1.equals(self))
        self.assertFalse(m1.equals({}))
        self.assertFalse(m1.equals(None))

    def test_default_models(self):
        m1 = models.DemographicModel.empty()
        m2 = models.DemographicModel.empty()
        self.assertTrue(m1.equals(m2))

    def test_migration_matrices(self):
        m1 = models.DemographicModel.empty()
        m2 = models.DemographicModel.empty()
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


class TestIsolationWithMigration(unittest.TestCase):
    """
    Tests for the generic IM model.
    """
    def test_pop_configs(self):
        model = models.IsolationWithMigration(100, 200, 300, 50, 0, 0)
        self.assertEqual(len(model.population_configurations), 3)
        self.assertEqual(model.population_configurations[0].initial_size, 200)
        self.assertEqual(model.population_configurations[1].initial_size, 300)
        self.assertEqual(model.population_configurations[2].initial_size, 100)

    def test_split(self):
        model = models.IsolationWithMigration(100, 200, 300, 50, 0, 0)
        self.assertEqual(len(model.demographic_events), 2)
        self.assertEqual(model.demographic_events[0].time, 50)
        self.assertEqual(model.demographic_events[1].time, 50)

    def test_migration_rates(self):
        model = models.IsolationWithMigration(100, 200, 300, 50, 0.002, 0.003)
        self.assertEqual(np.shape(model.migration_matrix), (3, 3))
        self.assertEqual(model.migration_matrix[0][1], 0.002)
        self.assertEqual(model.migration_matrix[1][0], 0.003)
        # check these are the only two nonzero entries in the migration matrix
        self.assertEqual(np.sum(np.array(model.migration_matrix) != 0), 2)


class TestPopulationSampling(unittest.TestCase):
    # Create populations to test on
    _pop1 = stdpopsim.Population("pop0", "Test pop. 0")
    _pop2 = stdpopsim.Population("pop1", "Test pop. 1",
                                 sampling_time=10)
    _pop3 = stdpopsim.Population("pop2", "Test pop. 2",
                                 sampling_time=None)

    # Create an empty model to hold populations
    base_mod = models.DemographicModel.empty(populations=[_pop1, _pop2, _pop3])

    def test_num_sampling_populations(self):
        self.assertEqual(self.base_mod.num_sampling_populations, 2)

    def test_get_samples(self):
        test_samples = self.base_mod.get_samples(2, 1)
        self.assertEqual(len(test_samples), 3)
        # Check for error when prohibited sampling asked for
        with self.assertRaises(ValueError):
            self.base_mod.get_samples(2, 2, 1)
        # Get the population corresponding to each sample
        sample_populations = [i.population for i in test_samples]
        # Check sample populations
        self.assertEqual(sample_populations, [0, 0, 1])
        # Test sampling times
        sample_times = [i.time for i in test_samples]
        self.assertEqual(sample_times, [0, 0, 10])

    # Test that all sampling populations are specified before non-sampling populations
    # in the model.populations list
    def test_population_order(self):
        for model in stdpopsim.all_demographic_models():
            allow_sample_status = [int(p.allow_samples) for p in model.populations]
            num_sampling = sum(allow_sample_status)
            # All sampling populations must be at the start of the list
            self.assertEqual(sum(allow_sample_status[num_sampling:]), 0)

    # Test that populations are listed in the same order in model.populations and
    # model.population_configurations
    def test_population_config_order_equal(self):
        for model in stdpopsim.all_demographic_models():
            pop_names = [pop.name for pop in model.populations]
            config_names = [
                config.metadata["name"] for config in model.population_configurations]
            for p, c in zip(pop_names, config_names):
                self.assertEqual(p, c)
