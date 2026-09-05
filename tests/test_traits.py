"""
Tests for infrastructure around traits.
"""

import pytest
import stdpopsim
import numpy as np
import copy
import textwrap


def check_trait_ids_errors(fun, args):
    for bad_list in ("foo", 123, ("a", "bc"), []):
        with pytest.raises(ValueError, match="Trait IDs must be"):
            fun(trait_ids=bad_list, **args)
    for bad_list in (["a", "a"], ["c", "abc", "d", "abc"]):
        with pytest.raises(ValueError, match="nonempty list of unique"):
            fun(trait_ids=bad_list, **args)

    tids = ["fitness", "height", "number_of_nostrils"]
    for bad_value in (123, None, ["xyz", "abc"], ""):
        tids[1] = bad_value
        with pytest.raises(ValueError, match="Each trait ID must be"):
            fun(trait_ids=tids, **args)


def check_particular_nonoverlapping_intervals_error(intervals, match):
    with pytest.raises(ValueError, match=match):
        stdpopsim.traits._check_nonoverlapping_intervals(intervals)
    tm = stdpopsim.TraitsModel(
        traits=[
            stdpopsim.Trait(id="height", type="additive"),
            stdpopsim.Trait(id="not_height", type="additive"),
        ]
    )
    with pytest.raises(ValueError, match=match):
        tm.add_environment(
            id="abc",
            trait_ids=["height", "not_height"],
            distribution_type="mvn",
            distribution_args=[np.zeros(2), np.eye(2)],
            time_intervals=intervals,
        )
    with pytest.raises(ValueError, match=match):
        tm.add_fitness_function(
            id="abc",
            trait_ids=["height", "not_height"],
            function_type="gaussian",
            function_args=[np.zeros(2), np.eye(2)],
            time_intervals=intervals,
        )


def test_check_nonoverlapping_intervals_errors():
    check_particular_nonoverlapping_intervals_error(5, "Intervals must be a list.")
    check_particular_nonoverlapping_intervals_error(
        [], "Cannot supply an empty list of intervals."
    )
    check_particular_nonoverlapping_intervals_error(
        [(1, 2, 3)], "Each interval in time_interval must be of the form"
    )
    check_particular_nonoverlapping_intervals_error(
        [("blee", 1)], "Intervals must be numeric."
    )
    check_particular_nonoverlapping_intervals_error(
        [(2.0, "blee")], "Intervals must be numeric."
    )
    check_particular_nonoverlapping_intervals_error(
        [(-3, 5)], "Intervals must start at the present"
    )
    check_particular_nonoverlapping_intervals_error(
        [(6, 2)], "Lower end of interval must be"
    )
    check_particular_nonoverlapping_intervals_error(
        [(0, 3), (2, 5)], "Intervals must be non-overlapping."
    )


def check_arg_copies(fun, args, check_arg):
    # To avoid gotchas, several of these methods should take a copy of their
    # arguments; test we aren't still referring to the original list
    args = copy.deepcopy(args)
    orig = args[check_arg].copy()
    x = fun(**args).__dict__
    assert args[check_arg] == orig
    assert x[check_arg] == orig
    assert args[check_arg][0] != 0, "test needs nonzero first entry"
    args[check_arg][0] *= 2
    assert x[check_arg] == orig
    assert x[check_arg] != args[check_arg]


class TestTraitsModel:

    def test_make_traits_model(self):
        tm = stdpopsim.TraitsModel([])
        assert len(tm.traits) == 1
        assert tm.traits[0].id == "fitness"
        assert tm.traits[0].type == "multiplicative"
        assert tm.environments == []
        assert tm.fitness_functions == []
        traits = [stdpopsim.Trait(id=u, type="additive") for u in "abc"]
        tm = stdpopsim.TraitsModel(traits=traits)
        assert len(tm.traits) == len(traits) + 1
        for i in range(len(traits)):
            assert tm.traits[i + 1] == traits[i]
        assert tm.environments == []
        assert tm.fitness_functions == []

    def test_unique_trait_ids(self):
        traits = [stdpopsim.Trait(id=u, type="additive") for u in "aba"]
        with pytest.raises(ValueError, match="must be unique"):
            stdpopsim.TraitsModel(traits=traits)

    def test_add_fitness_function(self):
        traits = [stdpopsim.Trait(id=u, type="additive") for u in "abc"]
        tm = stdpopsim.TraitsModel(traits=traits)
        fl = [
            {
                "id": u,
                "trait_ids": [u],
                "function_type": "gaussian",
                "function_args": [np.array([0]), np.array([[1]])],
                "population_list": [u + "_pop"],
                "time_intervals": [(0, end)],
            }
            for u, end in zip("abc", [1, 2, 3])
        ]
        for f in fl:
            tm.add_fitness_function(**f)
        for (u, end), ff in zip(zip("abc", [1, 2, 3]), tm.fitness_functions):
            assert ff.id == u
            assert ff.trait_ids == [u]
            assert ff.function_type == "gaussian"
            assert ff.function_args == [0, 1]
            assert ff.time_intervals == [(0, end)]
            assert ff.population_list == [u + "_pop"]

    def test_add_fitness_function_errors(self):
        traits = [stdpopsim.Trait(id=u, type="additive") for u in "abc"]
        tm = stdpopsim.TraitsModel(traits=traits)
        f = {
            "trait_ids": ["a"],
            "function_type": "gaussian",
            "function_args": [np.array([0]), np.array([[1]])],
        }
        tm.add_fitness_function(id="a", **f)
        with pytest.raises(ValueError, match="must be unique"):
            tm.add_fitness_function(id="a", **f)
        f["trait_ids"] = ["z"]
        with pytest.raises(ValueError, match="Unknown trait ID"):
            tm.add_fitness_function(id="z", **f)

    def test_add_environment(self):
        traits = [stdpopsim.Trait(id=u, type="additive") for u in "abc"]
        tm = stdpopsim.TraitsModel(traits=traits)
        el = [
            {
                "id": u,
                "trait_ids": [u],
                "distribution_type": "n",
                "distribution_args": [0, 1],
                "population_list": [u + "_pop"],
                "time_intervals": [(0, end)],
            }
            for u, end in zip("abc", [1, 2, 3])
        ]
        for e in el:
            tm.add_environment(**e)
        for (u, end), env in zip(zip("abc", [1, 2, 3]), tm.environments):
            assert env.id == u
            assert env.trait_ids == [u]
            assert env.distribution_type == "n"
            assert env.distribution_args == [0, 1]
            assert env.time_intervals == [(0, end)]
            assert env.population_list == [u + "_pop"]

    def test_add_environments_errors(self):
        traits = [stdpopsim.Trait(id=u, type="additive") for u in "abc"]
        tm = stdpopsim.TraitsModel(traits=traits)
        e = {"trait_ids": ["a"], "distribution_type": "n", "distribution_args": [0, 1]}
        tm.add_environment(id="a", **e)
        with pytest.raises(ValueError, match="already exists"):
            tm.add_environment(id="a", **e)
        e = {"trait_ids": ["z"], "distribution_type": "n", "distribution_args": [0, 1]}
        with pytest.raises(ValueError, match="Unknown trait ID"):
            tm.add_environment(id="z", **e)


class TestEnvironment:

    def test_make_environment(self):
        # one trait
        env = stdpopsim.Environment(
            id="abc",
            trait_ids=["height"],
            distribution_type="g",
            distribution_args=[1, 2],
        )
        assert env.id == "abc"
        assert env.trait_ids == ["height"]
        assert env.distribution_type == "g"
        assert env.distribution_args == [1, 2]
        assert env.population_list is None
        assert env.time_intervals is None
        # three traits
        tids = ["height", "boop", "num_nostrils"]
        env = stdpopsim.Environment(
            id="abc",
            trait_ids=tids,
            distribution_type="f",
            distribution_args=[1, 2, 3],
        )
        assert env.id == "abc"
        assert env.trait_ids == tids
        assert env.distribution_type == "f"
        assert env.distribution_args == [1, 2, 3]
        assert env.population_list is None
        assert env.time_intervals is None
        # population_list
        env = stdpopsim.Environment(
            id="abc",
            trait_ids=["height"],
            distribution_type="n",
            distribution_args=[1, 2],
            population_list=[0, 1, "pop_2"],
        )
        assert env.id == "abc"
        assert env.trait_ids == ["height"]
        assert env.distribution_type == "n"
        assert env.distribution_args == [1, 2]
        assert env.population_list == [0, 1, "pop_2"]
        assert env.time_intervals is None
        # time_intervals
        env = stdpopsim.Environment(
            id="abc",
            trait_ids=["height"],
            distribution_type="n",
            distribution_args=[1, 2],
            population_list=[0, 1, "pop_2"],
            time_intervals=[(0, 1), (1, 2.2), (4, float("inf"))],
        )
        assert env.id == "abc"
        assert env.trait_ids == ["height"]
        assert env.distribution_type == "n"
        assert env.distribution_args == [1, 2]
        assert env.population_list == [0, 1, "pop_2"]
        assert env.time_intervals == [(0, 1), (1, 2.2), (4, float("inf"))]

    def test_make_environment_copies(self):
        # making an environment should take a copy of its arguments;
        # test we aren't still referring to the original list
        args = {
            "id": "abc",
            "trait_ids": ["height", "boop", "num_nostrils"],
            "distribution_type": "f",
            "distribution_args": [1, 2, 3],
            "population_list": [1, 2, "pop_3"],
            "time_intervals": [(0, 1), (1, 2.2), (4, float("inf"))],
        }
        check_arg_copies(stdpopsim.Environment, args, "trait_ids")
        check_arg_copies(stdpopsim.Environment, args, "distribution_args")
        check_arg_copies(stdpopsim.Environment, args, "population_list")
        check_arg_copies(stdpopsim.Environment, args, "time_intervals")

    def test_make_environment_errors(self):
        with pytest.raises(TypeError, match="required keyword-only"):
            stdpopsim.Environment()
        for bad_id in (123, [], None):
            with pytest.raises(ValueError, match="id must be a nonempty string"):
                stdpopsim.Environment(
                    id=bad_id,
                    trait_ids=["foo"],
                    distribution_type="f",
                    distribution_args=[0],
                )
        with pytest.raises(TypeError, match="required keyword-only"):
            stdpopsim.Environment(
                id="abc",
            )
        with pytest.raises(TypeError, match="required keyword-only"):
            stdpopsim.Environment(
                id="abc",
                trait_ids=["abc"],
            )
        with pytest.raises(TypeError, match="required keyword-only"):
            stdpopsim.Environment(
                id="abc",
                distribution_type="f",
                distribution_args=[0],
            )
        with pytest.raises(ValueError, match="population_list contains"):
            stdpopsim.Environment(
                id="abc",
                trait_ids=["foo"],
                distribution_type="f",
                distribution_args=[0],
                population_list=[0, 1, 0],
            )
        for bad_pop in [2.2, [], set([]), ["pop"]]:
            with pytest.raises(ValueError, match="entries must be integers"):
                stdpopsim.Environment(
                    id="abc",
                    trait_ids=["foo"],
                    distribution_type="f",
                    distribution_args=[0],
                    population_list=[0, 1, bad_pop],
                )

    def test_trait_ids_errors(self):
        args = {
            "id": "abc",
            "distribution_type": "f",
            "distribution_args": [0],
        }
        check_trait_ids_errors(stdpopsim.Environment, args)

    def test_distribution_errors(self):
        args = {
            "id": "abc",
            "trait_ids": ["height", "boop", "num_nostrils"],
            "distribution_type": "f",
            "distribution_args": [1, 2, 3],
        }

        bad_args = args.copy()
        bad_args["distribution_type"] = 123
        with pytest.raises(ValueError, match="distribution_type must be str"):
            stdpopsim.Environment(**bad_args)

        bad_args = args.copy()
        bad_args["distribution_type"] = "e"
        with pytest.raises(ValueError, match="not implemented as a multivariate"):
            stdpopsim.Environment(**bad_args)

        bad_args = args.copy()
        for b in ([1, 2, 3, 4], []):
            bad_args["distribution_args"] = b
            with pytest.raises(ValueError, match="must be a list of length 3"):
                stdpopsim.Environment(**bad_args)

    def test_mvn_errors(self):
        args = {
            "id": "abc",
            "trait_ids": ["trait1", "trait2", "trait3"],
            "distribution_type": "mvn",
            "distribution_args": [
                np.array([0, 1, 2]),
                np.array([[1, 0.5, 0], [0.5, 2, 0], [0, 0, 3]]),
            ],
        }
        f = stdpopsim.Environment(**args)
        assert f.id == "abc"
        assert f.trait_ids == ["trait1", "trait2", "trait3"]
        assert f.distribution_type == "mvn"

        bad_args = args.copy()
        bad_args["distribution_args"] = [np.array([0]), np.array([1]), np.array([2])]
        with pytest.raises(ValueError, match="requires two parameters"):
            stdpopsim.Environment(**bad_args)

        bad_args = args.copy()
        bad_args["distribution_args"] = [np.array([0, 1]), args["distribution_args"][1]]
        with pytest.raises(ValueError, match="must be 1 dimensional of length 3"):
            stdpopsim.Environment(**bad_args)

        bad_args = args.copy()
        bad_args["distribution_args"] = [
            args["distribution_args"][0],
            np.array([0, 1, 2]),
        ]
        with pytest.raises(ValueError, match="must be 2 dimensional"):
            stdpopsim.Environment(**bad_args)

        bad_args = args.copy()
        bad_args["distribution_args"] = [
            args["distribution_args"][0],
            np.array([[1, 1, 0], [1, 1, 1]]),
        ]
        with pytest.raises(ValueError, match="must be square.* with dimensions .3, 3."):
            stdpopsim.Environment(**bad_args)

        bad_args = args.copy()
        S = args["distribution_args"][1]
        S[1, 2] = 10
        bad_args["distribution_args"] = [args["distribution_args"][0], S]
        with pytest.raises(ValueError, match="must be symmetric"):
            stdpopsim.Environment(**bad_args)

        bad_args = args.copy()
        S = args["distribution_args"][1]
        S[1, 2] = S[2, 1] = 10
        bad_args["distribution_args"] = [args["distribution_args"][0], S]
        with pytest.raises(ValueError, match="is not positive definite"):
            stdpopsim.Environment(**bad_args)


class TestTrait:

    def test_make_trait(self):
        for t, ta in [("identity", []), ("threshold", [-1]), ("liability", [-2, 3])]:
            tr = stdpopsim.Trait(
                id="num_bristles", type="additive", transform=t, transform_args=ta
            )
            assert tr.id == "num_bristles"
            assert tr.type == "additive"
            assert tr.transform == t
            assert tr.transform_args == ta

        tr = stdpopsim.Trait(
            id="num_bristles", type="multiplicative", transform=t, transform_args=ta
        )
        assert tr.id == "num_bristles"
        assert tr.type == "multiplicative"
        assert tr.transform == t
        assert tr.transform_args == ta

    def test_make_trait_copies(self):
        args = {
            "id": "num_bristles",
            "type": "additive",
            "transform": "threshold",
            "transform_args": [2],
        }
        check_arg_copies(stdpopsim.Trait, args, "transform_args")

    def test_make_trait_errors(self):
        with pytest.raises(TypeError, match="required .* argument"):
            stdpopsim.Trait()
        for bad_id in (123, [], None):
            with pytest.raises(ValueError, match="id must be a nonempty string"):
                stdpopsim.Trait(id=bad_id, type="additive")
        for bad_type in (123, [], None):
            with pytest.raises(ValueError, match="Unknown trait type"):
                stdpopsim.Trait(id="foo", type=bad_type)
        for bad_transform in (123, [], None):
            with pytest.raises(ValueError, match="Transform .* unknown"):
                stdpopsim.Trait(
                    id="foo",
                    type="additive",
                    transform=bad_transform,
                    transform_args=[],
                )
        for bad_args in ("abc", 123):
            with pytest.raises(ValueError, match="transform_args must be"):
                stdpopsim.Trait(
                    id="foo",
                    type="additive",
                    transform="identity",
                    transform_args=bad_args,
                )

    def test_type_errors(self):
        for bad_type in ("abc", 123, None, ""):
            with pytest.raises(ValueError, match="Unknown trait type"):
                stdpopsim.Trait(
                    id="foo", type=bad_type, transform="identity", transform_args=[]
                )

    def test_identity_errors(self):
        for bad_args in ([0], [-1]):
            with pytest.raises(ValueError, match="identity transform"):
                stdpopsim.Trait(
                    id="foo",
                    type="additive",
                    transform="identity",
                    transform_args=bad_args,
                )

    def test_threshold_errors(self):
        for bad_args in (None, [], [-1, 2]):
            with pytest.raises(ValueError, match="threshold transform"):
                stdpopsim.Trait(
                    id="foo",
                    type="additive",
                    transform="threshold",
                    transform_args=bad_args,
                )

    def test_liability_errors(self):
        for bad_args in (None, [], [1], [-1, 2, 3]):
            with pytest.raises(ValueError, match="liability transform"):
                stdpopsim.Trait(
                    id="foo",
                    type="additive",
                    transform="liability",
                    transform_args=bad_args,
                )
        with pytest.raises(ValueError, match="must be positive"):
            stdpopsim.Trait(
                id="foo", type="additive", transform="liability", transform_args=[0, -2]
            )

    def test_unknown_transform(self):
        with pytest.raises(ValueError, match="Transform .* unknown"):
            stdpopsim.Trait(
                id="foo", type="additive", transform="xyz", transform_args=[]
            )


class TestFitnessFunction:

    def test_make_fitness_function(self):
        f = stdpopsim.FitnessFunction(
            id="foo",
            trait_ids=["fitness"],
            function_type="gaussian",
            function_args=[np.array([0]), np.array([[1]])],
        )
        assert f.id == "foo"
        assert f.trait_ids == ["fitness"]
        assert f.function_type == "gaussian"
        assert f.function_args == [0, 1]
        # population_list
        f = stdpopsim.FitnessFunction(
            id="foo",
            trait_ids=["height"],
            function_type="gaussian",
            function_args=[np.array([0]), np.array([[1]])],
            population_list=[0, 1, "pop_2"],
        )
        assert f.id == "foo"
        assert f.trait_ids == ["height"]
        assert f.function_type == "gaussian"
        assert f.function_args == [0, 1]
        assert f.population_list == [0, 1, "pop_2"]
        assert f.time_intervals is None
        # time_intervals
        f = stdpopsim.FitnessFunction(
            id="foo",
            trait_ids=["height"],
            function_type="gaussian",
            function_args=[np.array([0]), np.array([[1]])],
            population_list=[0, 1, "pop_2"],
            time_intervals=[(0, 1), (1, 2.2), (4, float("inf"))],
        )
        assert f.id == "foo"
        assert f.trait_ids == ["height"]
        assert f.function_type == "gaussian"
        assert f.function_args == [0, 1]
        assert f.population_list == [0, 1, "pop_2"]
        assert f.time_intervals == [(0, 1), (1, 2.2), (4, float("inf"))]

    def test_make_fitness_function_errors(self):
        with pytest.raises(TypeError, match="required keyword-only arg"):
            stdpopsim.FitnessFunction()
        with pytest.raises(ValueError, match="Unknown function type foo."):
            stdpopsim.FitnessFunction(
                id="ff", trait_ids=["fitness"], function_type="foo", function_args=[0]
            )
        args = {
            "id": "foo",
            "trait_ids": ["fitness"],
            "function_type": "gaussian",
            "function_args": [np.array([0]), np.array([[1]])],
        }
        for bad_arg in [0, [], None, id]:
            args["function_type"] = bad_arg
            with pytest.raises(ValueError, match="must be a str"):
                stdpopsim.FitnessFunction(**args)
        args["function_type"] = "gaussian"
        for bad_id in (123, [], None):
            args["id"] = bad_id
            with pytest.raises(ValueError, match="id must be a nonempty string"):
                stdpopsim.FitnessFunction(**args)
        with pytest.raises(ValueError, match="population_list contains"):
            stdpopsim.FitnessFunction(
                id="abc",
                trait_ids=["foo"],
                function_type="gaussian",
                function_args=[np.ones(1), np.ones(1)],
                population_list=[0, 1, 0],
            )
        for bad_pop in [2.2, [], set([]), ["pop"]]:
            with pytest.raises(ValueError, match="entries must be integers"):
                stdpopsim.FitnessFunction(
                    id="abc",
                    trait_ids=["foo"],
                    function_type="gaussian",
                    function_args=[np.ones(1), np.ones(1)],
                    population_list=[0, 1, bad_pop],
                )

    def test_trait_ids_errors(self):
        args = {
            "id": "foo",
            "function_type": "gaussian",
            "function_args": [np.array([0]), np.array([[1]])],
        }
        check_trait_ids_errors(stdpopsim.FitnessFunction, args)

    def test_make_fitness_function_copies(self):
        args = {
            "id": "foo",
            "trait_ids": ["abc"],
            "function_type": "threshold",
            "function_args": [1, 2, 3],
            "population_list": [1, 2, "pop_3"],
            "time_intervals": [(0, 1), (1, 2.2), (4, float("inf"))],
        }
        check_arg_copies(stdpopsim.FitnessFunction, args, "trait_ids")
        check_arg_copies(stdpopsim.FitnessFunction, args, "function_args")
        check_arg_copies(stdpopsim.FitnessFunction, args, "population_list")
        check_arg_copies(stdpopsim.FitnessFunction, args, "time_intervals")

    def test_threshold_errors(self):
        args = {
            "id": "foo",
            "trait_ids": ["trait1"],
            "function_type": "threshold",
            "function_args": [0.1, 0, 2],
        }
        f = stdpopsim.FitnessFunction(**args)
        assert f.id == "foo"
        assert f.trait_ids == ["trait1"]
        assert f.function_type == "threshold"
        args["function_args"] = [2, 3]
        with pytest.raises(ValueError, match="takes three arguments"):
            stdpopsim.FitnessFunction(**args)
        for bad_val in [-1, 10]:
            args["function_args"] = [bad_val, 5, 3]
            with pytest.raises(ValueError, match="must be between 0 and 1"):
                stdpopsim.FitnessFunction(**args)
        for bad_val in [
            -1,
        ]:
            args["function_args"] = [0.5, bad_val, 3]
            with pytest.raises(ValueError, match="must be nonnegative"):
                stdpopsim.FitnessFunction(**args)
            args["function_args"] = [0.5, 2, bad_val]
            with pytest.raises(ValueError, match="must be nonnegative"):
                stdpopsim.FitnessFunction(**args)

    def test_gaussian_errors(self):
        args = {
            "id": "foo",
            "trait_ids": ["trait1", "trait2", "trait3"],
            "function_type": "gaussian",
            "function_args": [
                np.array([0, 1, 2]),
                np.array([[1, 0.5, 0], [0.5, 2, 0], [0, 0, 3]]),
            ],
        }
        f = stdpopsim.FitnessFunction(**args)
        assert f.id == "foo"
        assert f.trait_ids == ["trait1", "trait2", "trait3"]
        assert f.function_type == "gaussian"

        bad_args = args.copy()
        bad_args["function_args"] = [np.array([0]), np.array([1]), np.array([2])]
        with pytest.raises(ValueError, match="requires two parameters"):
            stdpopsim.FitnessFunction(**bad_args)

        bad_args = args.copy()
        bad_args["function_args"] = [np.array([0, 1]), args["function_args"][1]]
        with pytest.raises(ValueError, match="must be 1 dimensional of length 3"):
            stdpopsim.FitnessFunction(**bad_args)

        bad_args = args.copy()
        bad_args["function_args"] = [args["function_args"][0], np.array([0, 1, 2])]
        with pytest.raises(ValueError, match="must be 2 dimensional"):
            stdpopsim.FitnessFunction(**bad_args)

        # do we really want to require this?
        bad_args = args.copy()
        bad_args["function_args"] = [
            args["function_args"][0],
            list(args["function_args"][1]),
        ]
        with pytest.raises(ValueError, match="is not a numpy array"):
            stdpopsim.FitnessFunction(**bad_args)

        bad_args = args.copy()
        bad_args["function_args"] = [
            args["function_args"][0],
            np.array([[1, 1, 0], [1, 1, 1]]),
        ]
        with pytest.raises(ValueError, match="must be square.* with dimensions .3, 3."):
            stdpopsim.FitnessFunction(**bad_args)

        bad_args = args.copy()
        S = args["function_args"][1]
        S[1, 2] = 10
        bad_args["function_args"] = [args["function_args"][0], S]
        with pytest.raises(ValueError, match="must be symmetric"):
            stdpopsim.FitnessFunction(**bad_args)

        bad_args = args.copy()
        S = args["function_args"][1]
        S[1, 2] = S[2, 1] = 10
        bad_args["function_args"] = [args["function_args"][0], S]
        with pytest.raises(ValueError, match="is not positive definite"):
            stdpopsim.FitnessFunction(**bad_args)


class TestCreateMutationType:
    """
    Tests for creating a MutationType instance.
    """

    def test_default_mutation_type(self):
        mt = stdpopsim.MutationType()
        assert mt.dominance_coeff == 0.5
        assert mt.distribution_type == "f"
        assert mt.distribution_args == [0]
        assert mt.convert_to_substitution is True
        assert mt.dominance_coeff_list is None
        assert mt.dominance_coeff_breaks is None

    def test_default_trait_ids(self):
        mt = stdpopsim.MutationType()
        assert mt.trait_ids == ["fitness"]

    def test_bad_trait_ids(self):
        check_trait_ids_errors(stdpopsim.MutationType, {})

    def test_Q_scaled_index(self):
        mut_params = {
            "f": ([], [0]),
            "e": ([1], [0]),
            "g": ([0.014, 0.19], [0]),
            "n": ([0.5, 1], [0, 1]),
            "ln": ([0.5, 1], []),
            "lp": ([0.5, 1], []),
            "u": ([-0.5, 1], []),
        }
        for t in mut_params:
            if t == "f":
                mt = stdpopsim.MutationType(
                    distribution_type=t,
                )
            else:
                mt = stdpopsim.MutationType(
                    distribution_type=t, distribution_args=mut_params[t][0]
                )
            assert mt.Q_scaled_index == mut_params[t][1]

    def test_create_bad_mutation_type_message(self):
        # dominance_coeff must be a number
        with pytest.raises(ValueError, match="dominance_coeff must be a number."):
            stdpopsim.MutationType(
                dominance_coeff="abc",
            )

        # distribution_type must be str
        with pytest.raises(ValueError, match="distribution_type must be str."):
            stdpopsim.MutationType(
                distribution_type=1,
            )

        # distribution_args must be list
        with pytest.raises(ValueError, match="distribution_args must be list."):
            stdpopsim.MutationType(distribution_args=dict())

        # elements in distribution_args must be numbers
        with pytest.raises(ValueError, match="is not a number."):
            stdpopsim.MutationType(
                distribution_type="g",
                distribution_args=[0.5, "1"],
            )

        # elements in distribution_args must be valid.
        with pytest.raises(ValueError, match="is an invalid parameter."):
            stdpopsim.MutationType(
                distribution_type="n",
                distribution_args=[1, np.inf],
            )

        # convert_to_substitution must be bool
        with pytest.raises(ValueError, match="convert_to_substitution must be bool."):
            stdpopsim.MutationType(
                convert_to_substitution=1,
            )

        for dc in [np.inf, np.nan, -np.inf]:
            with pytest.raises(
                ValueError, match=f"Invalid dominance coefficient {dc}."
            ):
                stdpopsim.MutationType(
                    dominance_coeff=dc,
                )

        # unsupported distribution type
        with pytest.raises(
            ValueError, match="abc is not a supported distribution type."
        ):
            stdpopsim.MutationType(
                distribution_type="abc",
            )

        # fixed-value selection coefficient
        with pytest.raises(ValueError, match="must be a list of length"):
            stdpopsim.MutationType(
                distribution_type="f",
                distribution_args=[1, 2],
            )

        # gamma-distributed selection coefficient
        with pytest.raises(ValueError, match="uses a .mean, shape."):
            stdpopsim.MutationType(
                distribution_type="g",
                distribution_args=[1],
            )

        with pytest.raises(ValueError, match="The shape parameter must be positive."):
            stdpopsim.MutationType(
                distribution_type="g",
                distribution_args=[1, -1],
            )

        # exponentially-distributed selection coefficients
        with pytest.raises(ValueError, match="uses a .mean"):
            stdpopsim.MutationType(
                distribution_type="e",
                distribution_args=[1, 2],
            )

        # normally-distributed selection coefficients
        with pytest.raises(ValueError, match="uses a .mean, sd."):
            stdpopsim.MutationType(
                distribution_type="n",
                distribution_args=[1, 2, 3],
            )

        with pytest.raises(ValueError, match="The sd parameter must be nonnegative."):
            stdpopsim.MutationType(
                distribution_type="n",
                distribution_args=[1, -1],
            )

        # Weibull-distributed selection coefficients
        with pytest.raises(ValueError, match="uses a .scale, shape. parameterisation."):
            stdpopsim.MutationType(
                distribution_type="w",
                distribution_args=[1, 2, 3, 4],
            )

        with pytest.raises(ValueError, match="The scale parameter must be positive."):
            stdpopsim.MutationType(
                distribution_type="w",
                distribution_args=[-1, 2],
            )

        with pytest.raises(ValueError, match="The shape parameter must be positive."):
            stdpopsim.MutationType(
                distribution_type="w",
                distribution_args=[1, -2],
            )

        # Uniformly-distributed selection coefficients
        for bad_args in ([0], [1, 2, 3], [3, -2]):
            with pytest.raises(ValueError, match="uses a .min, max"):
                stdpopsim.MutationType(
                    distribution_type="u",
                    distribution_args=bad_args,
                )

        # Lognormally-distributed selection coefficients
        for dt in ["lp", "ln"]:
            with pytest.raises(
                ValueError, match="uses a .meanlog, sdlog. parameterisation"
            ):
                stdpopsim.MutationType(
                    distribution_type=dt,
                    distribution_args=[1, 2, 3, 4],
                )

            with pytest.raises(
                ValueError, match="The sdlog parameter must be nonnegative."
            ):
                stdpopsim.MutationType(
                    distribution_type=dt,
                    distribution_args=[1, -2],
                )

    def test_mutation_type_is_neutral(self):
        mt = stdpopsim.MutationType()
        assert mt.is_neutral is True

        mt = stdpopsim.MutationType(
            distribution_type="g",
            distribution_args=[0.014, 0.19],
        )
        assert mt.is_neutral is False

        mt = stdpopsim.MutationType(
            distribution_type="f",
            distribution_args=[1],
        )
        assert not mt.is_neutral

        mt = stdpopsim.MutationType(
            distribution_type="e",
            distribution_args=[1],
        )
        assert not mt.is_neutral

    def test_mutation_types(self):
        mut_params = {
            "f": ([-0.1], [0], [0.1], [50]),
            "g": ([-0.1, 0.1], [0.1, 0.1], [50, 50]),
            "e": ([0.1], [10], [5000], [0]),
            "n": ([-0.1, 0.2], [0.1, 0.1], [50, 50]),
            "w": ([0.1, 0.2], [0.1, 0.1], [50, 50]),
            "u": ([0.1, 0.2], [0.1, 0.1], [-5, 50]),
            "lp": ([-0.1, 0.2], [0.1, 0.1], [50, 50]),
            "ln": ([-0.1, 0.2], [0.1, 0.1], [50, 50]),
        }
        for t in mut_params:
            for p in mut_params[t]:
                mt = stdpopsim.MutationType(distribution_type=t, distribution_args=p)
                if t in ("lp", "ln", "u"):
                    assert mt.distribution_type == "s"
                else:
                    assert mt.distribution_type == t
                    assert len(mt.distribution_args) == len(p)
                    for a, b in zip(mt.distribution_args, p):
                        assert a == b

    def test_bad_mutation_types(self):
        bad_mut_params = {
            "f": ([0.1, 0.2], [], [np.inf]),
            "g": ([], [0.1, 0], [0.1, -0.1], [0.1, 0.4, 0.5], [0.1, np.inf]),
            "e": ([], [0, 1], [0.1, 0.4, 0.5], [np.inf]),
            "n": ([], [0.1, -1], [0.1, 0.4, 0.5], [0.1], [0.3, np.inf]),
            "w": ([], [-0.1, 1], [0.1, -1], [0.1, 0.4, 0.5], [0.1], [np.inf, 2.3]),
            "u": ([], [-0.1, -0.5], [0.1, 0.4, 0.5], [0.1], [2.3, np.inf]),
            "lp": ([], [0.1, -1], [0.1, 0.4, 0.5], [0.1], [0.1, np.inf]),
            "ln": ([], [0.1, -1], [0.1, 0.4, 0.5], [0.1], [0.1, np.inf]),
        }
        for t in bad_mut_params:
            for p in bad_mut_params[t]:
                print(t, p)
                with pytest.raises(ValueError):
                    stdpopsim.MutationType(distribution_type=t, distribution_args=p)

    def test_convert_to_substitution(self):
        mt = stdpopsim.MutationType()
        assert mt.convert_to_substitution is True
        for c in (True, False):
            mt = stdpopsim.MutationType(convert_to_substitution=c)
            assert mt.convert_to_substitution == c

    def test_dominance_coeff(self):
        mt = stdpopsim.MutationType()
        assert mt.dominance_coeff == 0.5
        for dominance_coeff in (-10, 0, 0.5, 1, 50):
            mt = stdpopsim.MutationType(dominance_coeff=dominance_coeff)
            assert mt.dominance_coeff == dominance_coeff

    def test_dominance_coeff_list(self):
        for dcl, dcb in (
            ([-0.1, 0.7, 1.2], [-2.1, 1.0]),
            ([-0.1, -0.7], [-2.1]),
        ):
            mt = stdpopsim.MutationType(
                dominance_coeff_list=dcl, dominance_coeff_breaks=dcb
            )
            assert np.allclose(dcl, mt.dominance_coeff_list)
            assert np.allclose(dcb, mt.dominance_coeff_breaks)

    def test_pass_by_value(self):
        # make sure that for the arguments that are lists
        # we can't post-hoc modify them (and thus bypass validation)
        val = 0.5
        x = [val]
        mt = stdpopsim.MutationType(distribution_args=x)
        x[0] = 2 * val + 1
        assert mt.distribution_args[0] == val
        x = [val, val]
        mt = stdpopsim.MutationType(
            dominance_coeff_list=x, dominance_coeff_breaks=[0.0]
        )
        x[0] = 2 * val + 1
        assert mt.dominance_coeff_list[0] == val
        x = [val]
        mt = stdpopsim.MutationType(
            dominance_coeff_list=[0.0, 0.0], dominance_coeff_breaks=x
        )
        x[0] = 2 * val + 1
        assert mt.dominance_coeff_breaks[0] == val

    def test_bad_dominance_coeff(self):
        for dominance_coeff in (np.inf, np.nan, "abc", [], {}):
            with pytest.raises(ValueError, match="dominance.coeff"):
                stdpopsim.MutationType(dominance_coeff=dominance_coeff)

    def test_bad_distribution_type(self):
        for distribution_type in (1, {}, None, "~", "!", "F"):
            with pytest.raises(ValueError):
                stdpopsim.MutationType(distribution_type=distribution_type)

    def test_bad_dominance_coeff_list(self):
        dcl = [-0.1, 0.7, 1.2]
        dcb = [-2.1, 1.0]
        # can't specify both dominance_coeff and list
        with pytest.raises(ValueError, match="both dominance_coeff and"):
            stdpopsim.MutationType(
                dominance_coeff=0.5,
                dominance_coeff_list=dcl,
                dominance_coeff_breaks=dcb,
            )
        with pytest.raises(ValueError, match="both dominance_coeff and"):
            stdpopsim.MutationType(
                dominance_coeff=0.5,
                dominance_coeff_list=dcl,
            )
        with pytest.raises(ValueError, match="both dominance_coeff and"):
            stdpopsim.MutationType(
                dominance_coeff=0.5,
                dominance_coeff_breaks=dcb,
            )
        # must have both coeffs and breaks
        with pytest.raises(ValueError, match="dominance.*no breaks"):
            stdpopsim.MutationType(dominance_coeff_list=dcl)
        # must have at least 2 bins
        with pytest.raises(ValueError, match="dominance.*at least 2"):
            stdpopsim.MutationType(
                dominance_coeff_list=[0.2],
                dominance_coeff_breaks=[],
            )
        # list must be one longer than breaks
        for x in ([], [0.0], dcl):
            with pytest.raises(ValueError, match="dominance.*equal to"):
                stdpopsim.MutationType(
                    dominance_coeff_list=dcl,
                    dominance_coeff_breaks=x,
                )
        # bad coefficients
        for x in (np.inf, np.nan, "abc", [], {}):
            with pytest.raises(ValueError, match="dominance.coeff"):
                stdpopsim.MutationType(
                    dominance_coeff_list=[x] + dcl[1:],
                    dominance_coeff_breaks=dcb,
                )
        # bad breaks
        for x in (np.inf, np.nan, "abc", [], {}):
            with pytest.raises(ValueError, match="dominance.*break"):
                stdpopsim.MutationType(
                    dominance_coeff_list=dcl,
                    dominance_coeff_breaks=[x] + dcb[1:],
                )
        with pytest.raises(ValueError, match="nondecreasing"):
            stdpopsim.MutationType(
                dominance_coeff_list=dcl,
                dominance_coeff_breaks=list(reversed(dcb)),
            )
        # must only affect fitness
        with pytest.raises(ValueError, match="non-fitness traits."):
            stdpopsim.MutationType(
                trait_ids=["height"],
                distribution_type="f",
                distribution_args=[1.0],
                dominance_coeff_list=dcl,
                dominance_coeff_breaks=dcb,
            )
        with pytest.raises(ValueError, match="non-fitness traits."):
            stdpopsim.MutationType(
                trait_ids=["fitness", "height"],
                distribution_type="f",
                distribution_args=[1.0, 2.0],
                dominance_coeff_list=dcl,
                dominance_coeff_breaks=dcb,
            )

    def test_copies(self):
        # making an mutation type should take a copy of its arguments;
        # test we aren't still referring to the original list
        args = [1, 2, 3]
        mt = stdpopsim.MutationType(
            trait_ids=["a", "b", "c"],
            distribution_type="f",
            distribution_args=args,
        )
        assert mt.distribution_args == args
        orig_args = args.copy()
        args[0] = 99
        assert mt.distribution_args == orig_args
        assert mt.distribution_args != args


class TestDFEOutput:
    """
    Tests for DFE ouputs.
    """

    def test_str(self):
        d = stdpopsim.DFE(
            id="xyz",
            description="abc",
            long_description="ABC",
            mutation_types=[stdpopsim.MutationType(convert_to_substitution=True)],
            proportions=[1.0],
        )
        ds = str(d)
        assert "xyz" in ds
        assert "abc" in ds
        assert "ABC" in ds

    def test_wrap_long_lines(self):
        d = stdpopsim.DFE(
            id="xyz",
            description="abc",
            long_description="ABC " * 10,
            mutation_types=[stdpopsim.MutationType(convert_to_substitution=True)],
            proportions=[1.0],
        )
        ds = str(d)
        expected = """\
        DFE:
        ║  id               = xyz
        ║  description      = abc
        ║  long_description = ABC ABC ABC ABC ABC ABC ABC ABC ABC ABC
        ║  citations        = []
        """
        assert textwrap.dedent(expected) in ds


class TestCreateDFE:
    """
    Tests for creating a DFE instance.
    """

    def test_default_dfe(self):
        d = stdpopsim.DFE(
            id="test",
            description="test",
            long_description="test",
        )
        assert d.mutation_types == []
        assert d.proportions == []
        assert d.citations == []
        assert d.qc_dfe is None
        assert d.is_neutral is True

    def test_basic_dfe(self):
        desc = "test test"
        long_desc = "test test 🐢"
        for props in ([0.4, 0.6], [1.0, 0.0], [1.0], [1 / 3, 1 / 3, 1 / 3]):
            mt = [stdpopsim.MutationType() for _ in props]
            d = stdpopsim.DFE(
                id="test",
                description=desc,
                long_description=long_desc,
                citations=["555"],
                mutation_types=mt,
                proportions=props,
            )
            assert d.id == "test"
            assert d.description == desc
            assert d.long_description == long_desc
            assert len(d.citations) == 1
            assert d.citations[0] == "555"
            for a, b in zip(mt, d.mutation_types):
                assert a == b
            for a, b in zip(props, d.proportions):
                assert a == b
            assert d.qc_dfe is None
            assert d.is_neutral is True

    def test_create_dfe_without_citations(self):
        d = stdpopsim.DFE(
            id="test",
            description="test",
            long_description="test",
            mutation_types=[stdpopsim.MutationType(), stdpopsim.MutationType()],
            proportions=[0.5, 0.5],
        )
        assert d.citations == []

    def test_create_dfe_without_mutation_types(self):
        d = stdpopsim.DFE(
            id="test",
            description="test",
            long_description="test",
        )

        assert d.mutation_types == []

    def test_create_dfe_without_proportions(self):
        d = stdpopsim.DFE(
            id="test",
            description="test",
            long_description="test",
            mutation_types=[stdpopsim.MutationType()],
        )
        assert d.proportions == [1]

    def test_create_bad_dfes_message(self):
        with pytest.raises(TypeError, match="required .* arguments"):
            # id, description, long_description are required
            stdpopsim.DFE()

        with pytest.raises(
            ValueError, match="proportions must be a list or numpy array."
        ):
            # proportions must be a list
            stdpopsim.DFE(
                id="test", description="test", long_description="test", proportions=1
            )

    def test_dfe_is_neutral(self):
        d = stdpopsim.DFE(
            id="test",
            description="test",
            long_description="test",
            mutation_types=[stdpopsim.MutationType(), stdpopsim.MutationType()],
            proportions=[0.5, 0.5],
        )
        assert d.is_neutral is True

        d = stdpopsim.DFE(
            id="test",
            description="test",
            long_description="test",
            mutation_types=[],
            proportions=[],
        )
        assert d.is_neutral is True

        for neutral in (True, False):
            for dist in ("f", "e"):
                props = [0.3, 0.7]
                if neutral:
                    svals = [0.0, 0.0]
                else:
                    svals = [0.0, 0.1]
                mt = [
                    stdpopsim.MutationType(
                        distribution_type=dist, distribution_args=[s]
                    )
                    for s in svals
                ]
                d = stdpopsim.DFE(
                    id=0,
                    description="test",
                    long_description="test test",
                    proportions=props,
                    mutation_types=mt,
                )
                assert d.is_neutral is (neutral and dist == "f")

    def test_no_msprime_dfe(self):
        # test we cannot simulate a non-neutral DFE with msprime
        m1 = stdpopsim.MutationType(
            dominance_coeff=0.2,
            distribution_type="e",
            distribution_args=[0.1],
        )
        desc = "test test"
        long_desc = "test test 🐢"
        d = stdpopsim.DFE(
            id="abc",
            description=desc,
            long_description=long_desc,
            mutation_types=[m1],
        )
        contig = stdpopsim.Contig.basic_contig(
            length=10000,
            mutation_rate=1e-6,
            ploidy=2,
        )
        contig.add_dme(
            intervals=np.array([[0, contig.length / 2]], dtype="int"),
            DME=d,
        )
        model = stdpopsim.PiecewiseConstantSize(1000)
        samples = {"pop_0": 1}
        engine = stdpopsim.get_engine("msprime")
        with pytest.raises(ValueError, match="but you are using .* msprime"):
            _ = engine.simulate(
                model,
                contig,
                samples,
            )


class TestCreateDME:
    """
    Tests for DistributionOfMutationEffects.
    """

    def test_default_dme(self):
        dme = stdpopsim.DistributionOfMutationEffects()
        assert dme.mutation_types == []
        assert dme.proportions == []

    def test_create_dme(self):
        mt1 = stdpopsim.MutationType()
        mt2 = stdpopsim.MutationType(
            trait_ids=["abc", "xyz"],
        )
        mt3 = stdpopsim.MutationType(
            trait_ids=["abc", "xyz"],
            distribution_type="f",
            distribution_args=[1, 2],
        )
        for mts, props in [
            ([mt1], [1]),
            ([mt1, mt2, mt3], [0.5, 0.25, 0.25]),
            ([mt1, mt2, mt3], [1 / 3, 2 / 3, 0]),
        ]:
            dme = stdpopsim.DistributionOfMutationEffects(
                mutation_types=mts, proportions=props
            )
            assert dme.mutation_types == mts
            assert dme.proportions == props

    def test_mutation_types_proportions_errors(self):
        with pytest.raises(ValueError, match="mutation_types must be a list."):
            # mutation_types must be a list
            stdpopsim.DistributionOfMutationEffects(
                mutation_types=stdpopsim.MutationType(),
                proportions=[0.5, 0.5],
            )

        with pytest.raises(
            ValueError,
            match="proportions and mutation_types must be lists of the same length.",
        ):
            # proportions and mutation_types must be the same length
            stdpopsim.DistributionOfMutationEffects(
                mutation_types=[stdpopsim.MutationType(), stdpopsim.MutationType()],
                proportions=[1.0],
            )

        with pytest.raises(
            ValueError, match="proportions must be nonnegative numbers."
        ):
            # proportions must be numbers
            stdpopsim.DistributionOfMutationEffects(
                mutation_types=[stdpopsim.MutationType(), stdpopsim.MutationType()],
                proportions=[True, "666"],
            )

        with pytest.raises(
            ValueError, match="proportions must be nonnegative numbers."
        ):
            # proportions must be positive
            stdpopsim.DistributionOfMutationEffects(
                mutation_types=[stdpopsim.MutationType(), stdpopsim.MutationType()],
                proportions=[-1, 0],
            )

        with pytest.raises(ValueError, match="proportions must sum to 1.0."):
            # proportions must sum 1.0
            stdpopsim.DistributionOfMutationEffects(
                mutation_types=[stdpopsim.MutationType(), stdpopsim.MutationType()],
                proportions=[1, 1],
            )

        with pytest.raises(
            ValueError, match="mutation_types must be a list of MutationType objects."
        ):
            # mutation_types must be a list of MutationType objects
            stdpopsim.DistributionOfMutationEffects(
                mutation_types=["neutral"],
                proportions=[1],
            )

        m1 = stdpopsim.MutationType()
        m2 = stdpopsim.MutationType()
        with pytest.raises(ValueError, match="must be lists of the same length"):
            stdpopsim.DistributionOfMutationEffects(
                proportions=[],
                mutation_types=[m1],
            )
        for bad_props in [
            ["abc"],
            [1.25, -0.25],
            1.0,
            [1.0],
            [0.2, 0.4, 0.4],
            [-0.1, -0.1],
            [0.8, 0.8],
        ]:
            with pytest.raises(ValueError):
                stdpopsim.DistributionOfMutationEffects(
                    proportions=bad_props,
                    mutation_types=[m1, m2],
                )
        for bad_mut_types in ["abc", {}, [1.0, 2.0], [m1], m1, ["a", "b"]]:
            with pytest.raises(ValueError):
                stdpopsim.DistributionOfMutationEffects(
                    proportions=[0.6, 0.4],
                    mutation_types=bad_mut_types,
                )
        for bad_sums in [[-0.4, 0.5], [0.6, 0.8], [139487135987, 0.0], [0.2, 0.3]]:
            with pytest.raises(ValueError):
                stdpopsim.DistributionOfMutationEffects(
                    proportions=bad_sums,
                    mutation_types=[m1, m2],
                )


class TestCreateNeutralDFE:
    """
    Tests for creating a neutral DFE instance.
    """

    def test_create_neutral_dfe(self):
        nd = stdpopsim.neutral_dfe()
        assert isinstance(nd, stdpopsim.DFE)
        assert nd.id == "neutral"
        assert nd.description == "neutral DFE"
        assert nd.long_description == "strictly neutral mutations"
        assert len(nd.mutation_types) == 1
        assert nd.mutation_types[0].is_neutral is True
        assert len(nd.proportions) == 1
        assert nd.proportions[0] == 1.0
        assert nd.is_neutral is True
