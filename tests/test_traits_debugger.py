"""
Tests for the traits debugger and for aligning traits models
with demographic models.
"""

import io

import msprime
import numpy as np
import pytest

import stdpopsim
from stdpopsim import traits


def split_model():
    # A and B split from anc 3000 generations ago.
    d = msprime.Demography()
    d.add_population(name="A", initial_size=1000)
    d.add_population(name="B", initial_size=500)
    d.add_population(name="anc", initial_size=800)
    d.add_population_split(time=3000, derived=["A", "B"], ancestral="anc")
    return stdpopsim.DemographicModel(
        id="test_model",
        description="test",
        long_description="test",
        generation_time=1,
        model=d,
    )


def traits_model(**condition_kwargs):
    tm = stdpopsim.TraitsModel(traits=[stdpopsim.Trait(id="add1", type="additive")])
    tm.add_fitness_function(
        id="fit1",
        trait_ids=["add1"],
        function_type="gaussian",
        function_args=[np.zeros(1), np.eye(1)],
        **condition_kwargs,
    )
    return tm


class TestCollectValidPopulationIntervals:
    def test_split_model(self):
        got = traits._collect_valid_population_intervals(split_model())
        assert got == {
            0: [[0, 3000]],
            1: [[0, 3000]],
            2: [[3000, float("inf")]],
        }

    def test_populations_inactive_at_present(self):
        # Ancestral populations must not be reported as active at recent
        # times, and ancient default sampling times must not hide the
        # recent activity of present-day populations.
        species = stdpopsim.get_species("HomSap")
        dm = species.get_demographic_model("AncientEurope_4A21")
        got = traits._collect_valid_population_intervals(dm)
        by_name = {p.name: p.id for p in dm.model.populations}
        assert got[by_name["OOA"]] == [[1500, float("inf")]]
        assert got[by_name["NE"]] == [[600, 1500]]
        assert got[by_name["Bronze"]] == [[0, 140]]


# All of the other parts of _standardize_condition are covered by other tests
def test_standardize_condition():
    valid_intervals = {
        0: [[0, 3000]],
        1: [[0, 3000]],
        2: [[3000, float("inf")]],
    }

    condition = traits_model().fitness_functions[0]
    condition.population_list = [4]
    with pytest.raises(ValueError, match="Population index out of bounds."):
        stdpopsim.traits._standardize_condition(condition, valid_intervals)

    condition = traits_model(
        population_list=[0], time_intervals=[[5, float("inf")]]
    ).fitness_functions[0]
    new_conditions = stdpopsim.traits._standardize_condition(condition, valid_intervals)
    assert new_conditions[0].time_intervals[0] == [5, 3000]

    condition = traits_model(
        population_list=[0], time_intervals=[[3001, float("inf")]]
    ).fitness_functions[0]
    with pytest.raises(ValueError, match="An environment or a fitness"):
        stdpopsim.traits._standardize_condition(condition, valid_intervals)


def test_resolve_population_ids():
    tm = traits_model(population_list=[0, 2, 1, 3])
    species = stdpopsim.get_species("HomSap")
    dm = species.get_demographic_model("AncientEurope_4A21")
    stdpopsim.traits._resolve_population_ids(tm, dm)
    assert tm.fitness_functions[0].population_list == [0, 2, 1, 3]

    tm = traits_model(population_list=["OOA", "WA", "NE", "CHG"])
    stdpopsim.traits._resolve_population_ids(tm, dm)
    assert tm.fitness_functions[0].population_list == [0, 2, 1, 3]

    tm = traits_model(population_list=[-3])
    with pytest.raises(ValueError, match="Population index out of bounds"):
        stdpopsim.traits._resolve_population_ids(tm, dm)

    tm = traits_model()
    # We have to set this manually, because intializing a fitness function
    # already checks for repeated indices.
    tm.fitness_functions[0].population_list = [0, 1, 2, 1]
    with pytest.raises(ValueError, match="Repeated population indices."):
        stdpopsim.traits._resolve_population_ids(tm, dm)


class TestAlignTraitsModelDemography:
    def align(self, tm):
        return traits._align_traits_model_demography(tm, split_model())

    def test_resolves_population_names(self):
        tm = self.align(traits_model(time_intervals=[(0, 1000)], population_list=["A"]))
        assert tm.fitness_functions[0].population_list == [0]

    def test_no_population_list_expanded(self):
        # A condition without populations applies to every population
        # active during its times, clipped to each one's activity.
        tm = self.align(traits_model(time_intervals=[(0, 4000)]))
        by_pop = {ff.population_list[0]: ff for ff in tm.fitness_functions}
        assert set(by_pop) == {0, 1, 2}
        assert by_pop[0].time_intervals == [[0, 3000]]
        assert by_pop[1].time_intervals == [[0, 3000]]
        assert by_pop[2].time_intervals == [[3000, 4000]]

    def test_no_population_list_drops_inactive(self):
        # anc is not active on [0, 1000), so it gets no copy.
        tm = self.align(traits_model(time_intervals=[(0, 1000)]))
        by_pop = {ff.population_list[0]: ff for ff in tm.fitness_functions}
        assert set(by_pop) == {0, 1}

    def test_no_population_list_no_times(self):
        tm = self.align(traits_model(time_intervals=None))
        by_pop = {ff.population_list[0]: ff for ff in tm.fitness_functions}
        assert set(by_pop) == {0, 1, 2}
        assert by_pop[0].time_intervals == [[0, 3000]]
        assert by_pop[2].time_intervals == [[3000, float("inf")]]

    def test_infinite_interval_clipped(self):
        tm = self.align(
            traits_model(time_intervals=[(0, float("inf"))], population_list=["B"])
        )
        assert tm.fitness_functions[0].time_intervals == [[0, 3000]]

    def test_finite_interval_outside_population(self):
        with pytest.raises(ValueError, match="does not exist"):
            self.align(traits_model(time_intervals=[(0, 4000)], population_list=["B"]))

    def test_bad_population_name(self):
        with pytest.raises(ValueError, match="not in demographic model"):
            self.align(traits_model(population_list=["nope"]))

    def test_bad_population_index(self):
        with pytest.raises(ValueError, match="out of bounds"):
            self.align(traits_model(population_list=[3]))


class TestTraitsDebugger:
    def epochs_text(self, dbg):
        # Split the printed output into one chunk per epoch box.
        chunks = str(dbg).split("Epoch[")
        return chunks[0], chunks[1:]

    def test_no_traits_model(self):
        dm = split_model()
        dbg = stdpopsim.TraitsDebugger(dm)
        head, epochs = self.epochs_text(dbg)
        assert len(epochs) == len(dm.model.debug().epochs)
        assert "Burn-in" not in head

    def test_epochs_split_at_condition_boundaries(self):
        tm = traits_model(time_intervals=[(0, 1000)])
        tm.add_environment(
            id="env1",
            trait_ids=["add1"],
            distribution_type="mvn",
            distribution_args=[np.zeros(1), np.eye(1)],
            time_intervals=[(0, float("inf"))],
            population_list=["A"],
        )
        dbg = stdpopsim.TraitsDebugger(split_model(), tm)
        head, epochs = self.epochs_text(dbg)
        # [0, 1000), [1000, 3000), [3000, inf)
        assert len(epochs) == 3
        assert "fit1" in epochs[0]
        assert "fit1" not in epochs[1]
        assert "fit1" not in epochs[2]
        # env1 is clipped to the times A is active
        assert "env1" in epochs[0]
        assert "env1" in epochs[1]
        assert "env1" not in epochs[2]

    def test_infinite_interval_includes_burn_in(self):
        tm = traits_model(time_intervals=[(0, float("inf"))])
        _, epochs = self.epochs_text(stdpopsim.TraitsDebugger(split_model(), tm))
        assert all("fit1" in epoch for epoch in epochs)
        assert "includes burn-in" in epochs[-1]
        assert sum("includes burn-in" in epoch for epoch in epochs) == 1

    def test_burn_in_note(self):
        tm = traits_model(time_intervals=[(3000, 4000)], population_list=["anc"])
        dbg = stdpopsim.TraitsDebugger(split_model(), tm)
        head, epochs = self.epochs_text(dbg)
        assert "Burn-in ends 4e+03" in head
        assert "'fit1'" in head
        # the extra breakpoint at 4000 splits the oldest epoch
        assert len(epochs) == 3
        assert "fit1" not in epochs[0]
        assert "fit1" in epochs[1]
        assert "fit1" not in epochs[2]

    def test_input_model_not_mutated(self):
        tm = traits_model(time_intervals=[(0, float("inf"))], population_list=["A"])
        stdpopsim.TraitsDebugger(split_model(), tm)
        assert tm.fitness_functions[0].population_list == ["A"]
        assert tm.fitness_functions[0].time_intervals == [(0, float("inf"))]

    def test_print_history(self):
        dbg = stdpopsim.TraitsDebugger(split_model(), traits_model())
        buf = io.StringIO()
        dbg.print_history(buf)
        assert buf.getvalue() == str(dbg)

    def test_catalog_model(self):
        species = stdpopsim.get_species("HomSap")
        dm = species.get_demographic_model("AncientEurope_4A21")
        out = str(stdpopsim.TraitsDebugger(dm))
        assert "Bronze" in out
        assert "includes burn-in" in out
