"""
Tests for infrastructure around traits.
"""

import pytest
import stdpopsim
import numpy as np
import copy


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
        assert tm.traits == []
        assert tm.environments == []
        assert tm.fitness_functions == []
        traits = [stdpopsim.Trait(id=u, type="additive") for u in "abc"]
        tm = stdpopsim.TraitsModel(traits=traits)
        assert tm.traits == traits
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
            }
            for u in "abc"
        ]
        for f in fl:
            tm.add_fitness_function(**f)
        for u, ff in zip("abc", tm.fitness_functions):
            assert ff.id == u
            assert ff.trait_ids == [u]
            assert ff.function_type == "gaussian"
            assert ff.function_args == [0, 1]

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
            }
            for u in "abc"
        ]
        for e in el:
            tm.add_environment(**e)
        for u, env in zip("abc", tm.environments):
            assert env.id == u
            assert env.trait_ids == [u]
            assert env.distribution_type == "n"
            assert env.distribution_args == [0, 1]

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

    def test_make_environment_copies(self):
        # making an environment should take a copy of its arguments;
        # test we aren't still referring to the original list
        args = {
            "id": "abc",
            "trait_ids": ["height", "boop", "num_nostrils"],
            "distribution_type": "f",
            "distribution_args": [1, 2, 3],
        }
        check_arg_copies(stdpopsim.Environment, args, "trait_ids")
        check_arg_copies(stdpopsim.Environment, args, "distribution_args")

    def test_make_envionment_errors(self):
        with pytest.raises(TypeError, match="required keyword-only"):
            stdpopsim.Environment()
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

    def test_trait_ids_errors(self):
        args = {
            "id": "abc",
            "distribution_type": "f",
            "distribution_args": [0],
        }
        check_trait_ids_errors(stdpopsim.Environment, args)


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
        for bad_type in (123, [], None):
            with pytest.raises(ValueError, match="Unknown trait type"):
                stdpopsim.Trait(id="foo", type=bad_type)
        for bad_transform in (123, [], None):
            with pytest.raises(ValueError, match="transform must be"):
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

    def test_make_fitness_function_errors(self):
        with pytest.raises(TypeError, match="required keyword-only arg"):
            stdpopsim.FitnessFunction()

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
        }
        check_arg_copies(stdpopsim.FitnessFunction, args, "trait_ids")
        check_arg_copies(stdpopsim.FitnessFunction, args, "function_args")
