"""
Tests for infrastructure around traits.
"""

import pytest
import stdpopsim
import numpy as np


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
        tids = ["height", "boop", "num_nostrils"]
        env = stdpopsim.Environment(
            id="abc",
            trait_ids=tids,
            distribution_type="f",
            distribution_args=[1, 2, 3],
        )
        assert env.trait_ids == tids
        orig_tids = tids.copy()
        tids[0] = "sakdljfsdljlj"
        assert env.trait_ids == orig_tids
        assert env.trait_ids != tids

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
        with pytest.raises(ValueError, match="at least one trait"):
            stdpopsim.Environment(
                id="abc",
                trait_ids=[],
                distribution_type="f",
                distribution_args=[0],
            )


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
        # TODO
        pass

    def test_make_trait_errors(self):
        with pytest.raises(TypeError, match="required .* argument"):
            stdpopsim.Trait()
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
            trait_ids=["fitness"],
            function_type="gaussian",
            function_args=[np.array([0]), np.array([[1]])],
        )
        assert f.trait_ids == ["fitness"]
        assert f.function_type == "gaussian"
        assert f.function_args == [0, 1]

    def test_make_fitness_function_errors(self):
        with pytest.raises(TypeError, match="required keyword-only arg"):
            stdpopsim.FitnessFunction()
        with pytest.raises(ValueError, match="At least one"):
            stdpopsim.FitnessFunction(
                trait_ids=[],
                function_type="gaussian",
                function_args=[np.array([0]), np.array([[1]])],
            )

    def test_make_fitness_function_copies(self):
        # TODO
        pass
