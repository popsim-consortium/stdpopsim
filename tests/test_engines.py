"""
Tests for simulation engine infrastructure.
"""
import stdpopsim
import msprime
import pytest


class TestEngineAPI:
    """
    Tests for the API exposed for simulation engines.
    """

    def _test_engine(self, engine):
        assert isinstance(engine, stdpopsim.Engine)
        assert engine.simulate is not None
        assert len(engine.citations) != 0
        assert len(engine.get_version()) != 0

    def test_get_default_engine(self):
        engine = stdpopsim.get_default_engine()
        self._test_engine(engine)

    def test_all_engines(self):
        for engine in stdpopsim.all_engines():
            self._test_engine(engine)

    def test_register_engine(self):
        class MyEngine(stdpopsim.Engine):
            id = "test-engine"
            name = "test"
            citations = []

        engine1 = MyEngine()
        stdpopsim.register_engine(engine1)
        engine2 = stdpopsim.get_engine(engine1.id)
        assert engine1 == engine2
        # remove engine to avoid possible problems with other tests
        del stdpopsim.engines._registered_engines[engine1.id]

    def test_register_duplicate(self):
        engine = stdpopsim.get_default_engine()
        with pytest.raises(ValueError):
            stdpopsim.register_engine(engine)

    def test_get_engine(self):
        with pytest.raises(ValueError):
            stdpopsim.get_engine("nonexistent")

    def test_abstract_base_class(self):
        e = stdpopsim.Engine()
        with pytest.raises(NotImplementedError):
            e.simulate(None, None, None)
        with pytest.raises(NotImplementedError):
            e.get_version()


@pytest.mark.filterwarnings(
    "ignore:.*model has mutation rate.*but this simulation used.*"
)
class TestBehaviour:
    def test_simulate_nonexistent_param(self):
        species = stdpopsim.get_species("HomSap")
        model = species.get_demographic_model("AshkSub_7G19")
        good_kwargs = dict(
            demographic_model=model,
            contig=species.get_contig("chr1"),
            samples=model.get_samples(10, 10, 10),
            dry_run=True,
        )
        bad_kwargs = good_kwargs.copy().update(nonexistent_param=None)
        for engine in stdpopsim.all_engines():
            engine.simulate(**good_kwargs)
            with pytest.raises(TypeError):
                engine.simulate(**bad_kwargs)

    def test_required_params(self):
        species = stdpopsim.get_species("HomSap")
        model = species.get_demographic_model("AshkSub_7G19")
        contig = (species.get_contig("chr1"),)
        for engine in stdpopsim.all_engines():
            with pytest.raises(TypeError):
                engine.simulate(model, contig)

    def test_msprime_kwargs(self):
        species = stdpopsim.get_species("HomSap")
        model = species.get_demographic_model("AshkSub_7G19")
        contig = species.get_contig("chr22", length_multiplier=0.01)
        samples = model.get_samples(10)
        engine = stdpopsim.get_engine("msprime")
        sim_arg = engine.simulate(
            model, contig, samples, record_full_arg=True, random_seed=1
        )
        assert any(msprime.NODE_IS_RE_EVENT == sim_arg.tables.nodes.flags)

    def test_msprime_seed(self):
        species = stdpopsim.get_species("HomSap")
        model = species.get_demographic_model("AshkSub_7G19")
        contig = species.get_contig("chr22", length_multiplier=0.01)
        samples = model.get_samples(10)
        engine = stdpopsim.get_engine("msprime")
        with pytest.raises(ValueError):
            engine.simulate(model, contig, samples, seed=1, random_seed=1)
        sim_seed = engine.simulate(model, contig, samples, seed=1)
        sim_random_seed = engine.simulate(model, contig, samples, random_seed=1)
        assert sim_seed.tables.edges == sim_random_seed.tables.edges
