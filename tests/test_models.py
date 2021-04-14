"""
Tests for simulation model infrastructure.
"""
import unittest
import sys
import textwrap
import copy

import numpy as np
import msprime
import pytest

import stdpopsim
from stdpopsim import models
from stdpopsim import utils


class DemographicModelTestMixin(object):
    """
    Mixin for testing specific models. Subclasses should extend
    unittest.TestCase and this mixin, and define the self.model (as the
    model instance).
    """

    # To be defined in subclasses.
    model = None

    def test_debug_runs(self):
        dbg = self.model.model.debug()
        assert dbg.num_epochs > 0

    def test_simulation_runs(self):
        # With a recombination_map of None, we simulate a coalescent without
        # recombination in msprime, with mutation rate equal to rate from model.
        contig = stdpopsim.Contig(mutation_rate=self.model.mutation_rate)
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


class CatalogDemographicModelTestMixin(DemographicModelTestMixin):
    """
    Mixin for demographic models in the catalog.
    """

    def test_id_valid(self):
        self.assertTrue(utils.is_valid_demographic_model_id(self.model.id))


class QcdCatalogDemographicModelTestMixin(CatalogDemographicModelTestMixin):
    """
    Extends the tests to also check that the qc model is equal to
    the production model.
    """

    def test_qc_model_equal(self):
        d1 = self.model.model
        d2 = self.model.qc_model.model
        d1.assert_equivalent(d2, rel_tol=1e-5)
        assert d1 != d2

    def test_generation_time_match(self):
        g1 = self.model.generation_time
        g2 = self.model.qc_model.generation_time
        self.assertEqual(g1, g2)

    def test_mutation_rate_match(self):
        u1 = self.model.mutation_rate
        u2 = self.model.qc_model.mutation_rate
        self.assertEqual(u1, u2)


# Add model specific test classes, derived from one of the above.
qc_test_classes = []
for species in stdpopsim.all_species():
    for model in species.demographic_models:
        superclasses = [unittest.TestCase]
        if model.qc_model is not None:
            superclasses.append(QcdCatalogDemographicModelTestMixin)
        else:
            superclasses.append(CatalogDemographicModelTestMixin)
        classname = f"Test{species.id}{model.id}"
        cls = type(classname, tuple(superclasses), dict(model=model))
        qc_test_classes.append(cls)
# Basic sanity checks to double check that no errors get introduced
# that lead to these qc tests being skipped silently.
assert len(qc_test_classes) > 0
for cls in qc_test_classes:
    assert issubclass(cls, unittest.TestCase)
    # Insert the class into the current test module's namespace.
    setattr(sys.modules[__name__], cls.__name__, cls)


class TestRegisterQCModel(unittest.TestCase):
    def make_model(self, name):
        return models.DemographicModel(
            id=name,
            description=name,
            long_description=name,
            generation_time=1,
            model=msprime.Demography.isolated_model([1]),
        )

    def test_register_qc(self):
        model = self.make_model("test")
        model.register_qc(model)

    def test_already_registered(self):
        model = self.make_model("test")
        model.register_qc(model)
        with self.assertRaises(ValueError):
            model.register_qc(model)

    def test_bad_qc_models(self):
        model = self.make_model("test")
        for not_a_model in [None, 15, "Zigzag_1S14"]:
            with self.assertRaises(ValueError):
                model.register_qc(not_a_model)


class TestAllModels:
    """
    Tests for registered simulation models.
    """

    def test_non_empty(self):
        assert len(list(stdpopsim.all_demographic_models())) > 0

    @pytest.mark.parametrize("model", stdpopsim.all_demographic_models())
    def test_all_instances(self, model):
        assert isinstance(model, models.DemographicModel)
        assert len(model.id) > 0
        assert len(model.description) > 0
        assert len(model.long_description) > 0
        assert len(model.citations) > 0
        assert model.generation_time > 0
        if model.mutation_rate is not None:
            assert model.mutation_rate > 0
        model.model.validate()


class TestModelOutput:
    def test_str(self):
        model = stdpopsim.DemographicModel(
            id="xyz",
            description="abc",
            long_description="ABC",
            generation_time=1234,
            model=msprime.Demography.isolated_model([1]),
        )
        s = str(model)
        assert "xyz" in s
        assert "abc" in s
        assert "ABC" in s
        assert "1234" in s

    def test_wrap_long_lines(self):
        model = stdpopsim.DemographicModel(
            id="xyz",
            description="abc",
            long_description="ABC " * 50,
            generation_time=1234,
            model=msprime.Demography.isolated_model([1]),
        )
        s = str(model)
        expected = """\
        Demographic model:
        ║  id               = xyz
        ║  description      = abc
        ║  long_description = ABC ABC ABC ABC ABC ABC ABC ABC ABC ABC ABC ABC ABC ABC ABC ABC ABC
        ║                     ABC ABC ABC ABC ABC ABC ABC ABC ABC ABC ABC ABC ABC ABC ABC ABC ABC
        ║                     ABC ABC ABC ABC ABC ABC ABC ABC ABC ABC ABC ABC ABC ABC ABC ABC
        ║  generation_time  = 1234
        ║  citations        = []
        """  # noqa 501
        assert textwrap.dedent(expected) in s


class TestModelsEqual(unittest.TestCase):
    """
    Tests DemographicModel object equality comparison.
    """

    def test_known_models(self):
        other_model = msprime.Demography.isolated_model([1])
        for model in stdpopsim.all_demographic_models():
            model.model.assert_equivalent(model.model)
            model.model.assert_equal(model.model)
            assert not model.model.is_equivalent(other_model)


class TestConstantSizeModel(unittest.TestCase, DemographicModelTestMixin):
    model = models.PiecewiseConstantSize(100)


class TestTwoEpochModel(unittest.TestCase, DemographicModelTestMixin):
    model = models.PiecewiseConstantSize(100, (10, 10))


class TestPiecewiseConstantSize(unittest.TestCase):
    """
    Specific tests for the piecewise constant model.
    """

    def test_single_epoch(self):
        model = models.PiecewiseConstantSize(100)
        self.assertEqual(model.model.populations[0].initial_size, 100)
        self.assertEqual(model.model.num_populations, 1)
        self.assertEqual(model.model.num_events, 0)

    def test_two_epoch(self):
        model = models.PiecewiseConstantSize(50, (10, 100))
        self.assertEqual(model.model.populations[0].initial_size, 50)
        self.assertEqual(model.model.num_populations, 1)
        self.assertEqual(model.model.num_events, 1)
        event = model.model.events[0]
        self.assertEqual(event.time, 10)
        self.assertEqual(event.initial_size, 100)
        self.assertIsNone(event.growth_rate)

    def test_three_epoch(self):
        model = models.PiecewiseConstantSize(0.1, (0.1, 10), (0.2, 100))
        self.assertEqual(model.model.populations[0].initial_size, 0.1)
        self.assertEqual(model.model.num_populations, 1)
        self.assertEqual(model.model.num_events, 2)
        event = model.model.events[0]
        self.assertEqual(event.time, 0.1)
        self.assertEqual(event.initial_size, 10)
        self.assertIsNone(event.growth_rate)
        event = model.model.events[1]
        self.assertEqual(event.time, 0.2)
        self.assertEqual(event.initial_size, 100)
        self.assertIsNone(event.growth_rate)


class TestIsolationWithMigration(unittest.TestCase):
    """
    Tests for the generic IM model.
    """

    def test_pop_configs(self):
        model = models.IsolationWithMigration(100, 200, 300, 50, 0, 0).model
        self.assertEqual(model.num_populations, 3)
        self.assertEqual(model.populations[0].initial_size, 200)
        self.assertEqual(model.populations[1].initial_size, 300)
        self.assertEqual(model.populations[2].initial_size, 100)

    def test_split(self):
        model = models.IsolationWithMigration(100, 200, 300, 50, 0, 0).model
        self.assertEqual(len(model.events), 1)
        self.assertEqual(model.events[0].time, 50)

    def test_migration_rates(self):
        model = models.IsolationWithMigration(100, 200, 300, 50, 0.002, 0.003).model
        self.assertEqual(np.shape(model.migration_matrix), (3, 3))
        self.assertEqual(model.migration_matrix[0][1], 0.002)
        self.assertEqual(model.migration_matrix[1][0], 0.003)
        # check these are the only two nonzero entries in the migration matrix
        self.assertEqual(np.sum(np.array(model.migration_matrix) != 0), 2)


class TestPopulationSampling(unittest.TestCase):
    def make_model(self):
        # Create populations to test on
        _pop1 = stdpopsim.Population("pop0", "Test pop. 0")
        _pop2 = stdpopsim.Population("pop1", "Test pop. 1", sampling_time=10)
        _pop3 = stdpopsim.Population("pop2", "Test pop. 2", sampling_time=None)

        # Create an model to hold populations
        base_mod = models.DemographicModel(
            id="x",
            description="y",
            long_description="z",
            populations=[_pop1, _pop2, _pop3],
            population_configurations=[
                msprime.PopulationConfiguration(initial_size=1),
                msprime.PopulationConfiguration(initial_size=1),
                msprime.PopulationConfiguration(initial_size=1),
            ],
        )
        return base_mod

    def test_num_sampling_populations(self):
        base_mod = self.make_model()
        self.assertEqual(base_mod.num_sampling_populations, 2)

    def test_get_samples(self):
        base_mod = self.make_model()
        test_samples = base_mod.get_samples(2, 1)
        self.assertEqual(len(test_samples), 2)
        self.assertEqual(sum([ss.num_samples for ss in test_samples]), 3)

        # Check for error when prohibited sampling asked for
        with self.assertRaises(ValueError):
            base_mod.get_samples(2, 2, 1)
        # Get the population corresponding to each sample
        sample_populations = [i.population for i in test_samples]
        # Check sample populations
        self.assertEqual(sample_populations, [0, 1])
        # Test sampling times
        sample_times = [i.time for i in test_samples]
        self.assertEqual(sample_times, [0, 10])

    # Test that all sampling populations are specified before non-sampling populations
    # in the model.populations list
    def test_population_order(self):
        for model in stdpopsim.all_demographic_models():
            allow_sample_status = [int(p.allow_samples) for p in model.populations]
            num_sampling = sum(allow_sample_status)
            # All sampling populations must be at the start of the list
            self.assertEqual(sum(allow_sample_status[num_sampling:]), 0)


class TestZigZagWarning(unittest.TestCase):
    def test_zigzag_produces_warning(self):
        species = stdpopsim.get_species("HomSap")
        model = species.get_demographic_model("Zigzag_1S14")
        contig = species.get_contig("chr22")
        samples = model.get_samples(10)
        for engine in stdpopsim.all_engines():
            with self.assertWarnsRegex(UserWarning, "Zigzag_1S14"):
                engine.simulate(model, contig, samples, dry_run=True)


class TestMutationRates(unittest.TestCase):
    def test_mutation_rate_warning(self):
        species = stdpopsim.get_species("HomSap")
        model = copy.deepcopy(species.get_demographic_model("OutOfAfrica_3G09"))
        ###
        # unneeded after mutation rates are included in catalog!!
        model.mutation_rate = 2.35e-8
        ###
        contig = species.get_contig("chr22")
        samples = model.get_samples(10, 10, 10)
        for engine in stdpopsim.all_engines():
            with self.assertWarnsRegex(UserWarning, "mutation"):
                engine.simulate(model, contig, samples, dry_run=True)

    def test_mutation_rate_match(self):
        species = stdpopsim.get_species("HomSap")
        model = copy.deepcopy(species.get_demographic_model("OutOfAfrica_3G09"))
        ###
        # unneeded after mutation rates are included in catalog!!
        model.mutation_rate = 2.35e-8
        ###
        contig = species.get_contig("chr22")
        self.assertNotEqual(model.mutation_rate, contig.mutation_rate)
        contig = species.get_contig("chr22", mutation_rate=model.mutation_rate)
        self.assertEqual(model.mutation_rate, contig.mutation_rate)
        contig = species.get_contig(length=100)
        self.assertNotEqual(model.mutation_rate, contig.mutation_rate)
        contig = species.get_contig(length=100, mutation_rate=model.mutation_rate)
        self.assertEqual(model.mutation_rate, contig.mutation_rate)
