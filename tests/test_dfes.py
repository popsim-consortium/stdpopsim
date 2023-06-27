"""
Tests for simulation dfe infrastructure
"""
import sys
import pytest
import stdpopsim
import textwrap
import numpy as np
from stdpopsim import dfe
from stdpopsim import utils

IS_WINDOWS = sys.platform.startswith("win")


class TestCreateMutationType:
    """
    Tests for creating a MutationType instance.
    """

    def test_default_mutation_type(self):
        mt = dfe.MutationType()
        assert mt.dominance_coeff == 0.5
        assert mt.distribution_type == "f"
        assert mt.distribution_args == [0]
        assert mt.convert_to_substitution is True

    def test_Q_scaled_index(self):
        mut_params = {
            "f": ([], [0]),
            "e": ([1], [0]),
            "g": ([0.014, 0.19], [0]),
            "n": ([0.5, 1], [0, 1]),
            "ln": ([0.5, 1], []),
            "lp": ([0.5, 1], []),
        }
        for t in mut_params:
            print(f"{t}\t{mut_params[t]}")
            if t == "f":
                mt = dfe.MutationType(
                    distribution_type=t,
                )
            else:
                mt = dfe.MutationType(
                    distribution_type=t, distribution_args=mut_params[t][0]
                )
            assert mt.Q_scaled_index == mut_params[t][1]

    def test_create_bad_mutation_type_message(self):
        # dominance_coeff must be a number
        with pytest.raises(ValueError, match="dominance_coeff must be a number."):
            dfe.MutationType(
                dominance_coeff="abc",
            )

        # distribution_type must be str
        with pytest.raises(ValueError, match="distribution_type must be str."):
            dfe.MutationType(
                distribution_type=1,
            )

        # distribution_args must be list
        with pytest.raises(ValueError, match="distribution_args must be list."):
            dfe.MutationType(distribution_args=dict())

        # elements in distribution_args must be numbers
        with pytest.raises(ValueError, match="is not a number."):
            dfe.MutationType(
                distribution_type="g",
                distribution_args=[0.5, "1"],
            )

        # elements in distribution_args must be valid.
        with pytest.raises(ValueError, match="is an invalid parameter."):
            dfe.MutationType(
                distribution_type="n",
                distribution_args=[1, np.inf],
            )

        # convert_to_substitution must be bool
        with pytest.raises(ValueError, match="convert_to_substitution must be bool."):
            dfe.MutationType(
                convert_to_substitution=1,
            )

        for dc in [np.inf, np.nan, np.NINF]:
            with pytest.raises(
                ValueError, match=f"Invalid dominance coefficient {dc}."
            ):
                dfe.MutationType(
                    dominance_coeff=dc,
                )

        # unsupported distribution type
        with pytest.raises(
            ValueError, match="abc is not a supported distribution type."
        ):
            dfe.MutationType(
                distribution_type="abc",
            )

        # fixed-value selection coefficient
        with pytest.raises(ValueError, match="take a single"):
            dfe.MutationType(
                distribution_type="f",
                distribution_args=[1, 2],
            )

        # gamma-distributed selection coefficient
        with pytest.raises(ValueError, match="use a .mean, shape."):
            dfe.MutationType(
                distribution_type="g",
                distribution_args=[1],
            )

        with pytest.raises(ValueError, match="The shape parameter must be positive."):
            dfe.MutationType(
                distribution_type="g",
                distribution_args=[1, -1],
            )

        # exponentially-distributed selection coefficients
        with pytest.raises(ValueError, match="use a .mean"):
            dfe.MutationType(
                distribution_type="e",
                distribution_args=[1, 2],
            )

        # normally-distributed selection coefficients
        with pytest.raises(ValueError, match="use a .mean, sd."):
            dfe.MutationType(
                distribution_type="n",
                distribution_args=[1, 2, 3],
            )

        with pytest.raises(ValueError, match="The sd parameter must be nonnegative."):
            dfe.MutationType(
                distribution_type="n",
                distribution_args=[1, -1],
            )

        # Weibull-distributed selection coefficients
        with pytest.raises(ValueError, match="use a .scale, shape. parameterisation."):
            dfe.MutationType(
                distribution_type="w",
                distribution_args=[1, 2, 3, 4],
            )

        with pytest.raises(ValueError, match="The scale parameter must be positive."):
            dfe.MutationType(
                distribution_type="w",
                distribution_args=[-1, 2],
            )

        with pytest.raises(ValueError, match="The shape parameter must be positive."):
            dfe.MutationType(
                distribution_type="w",
                distribution_args=[1, -2],
            )

        # Lognormally-distributed selection coefficients
        for dt in ["lp", "ln"]:
            with pytest.raises(
                ValueError, match="use a .meanlog, sdlog. parameterisation"
            ):
                dfe.MutationType(
                    distribution_type=dt,
                    distribution_args=[1, 2, 3, 4],
                )

            with pytest.raises(
                ValueError, match="The sdlog parameter must be nonnegative."
            ):
                dfe.MutationType(
                    distribution_type=dt,
                    distribution_args=[1, -2],
                )

    def test_mutation_type_is_neutral(self):
        mt = dfe.MutationType()
        assert mt.is_neutral is True

        mt = dfe.MutationType(
            distribution_type="g",
            distribution_args=[0.014, 0.19],
        )
        assert mt.is_neutral is False

        mt = dfe.MutationType(
            distribution_type="f",
            distribution_args=[1],
        )
        assert not mt.is_neutral

        mt = dfe.MutationType(
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
            "lp": ([-0.1, 0.2], [0.1, 0.1], [50, 50]),
            "ln": ([-0.1, 0.2], [0.1, 0.1], [50, 50]),
        }
        for t in mut_params:
            for p in mut_params[t]:
                mt = dfe.MutationType(distribution_type=t, distribution_args=p)
                if t in ("lp", "ln"):
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
            "lp": ([], [0.1, -1], [0.1, 0.4, 0.5], [0.1], [0.1, np.inf]),
            "ln": ([], [0.1, -1], [0.1, 0.4, 0.5], [0.1], [0.1, np.inf]),
        }
        for t in bad_mut_params:
            for p in bad_mut_params[t]:
                print(t, p)
                with pytest.raises(ValueError):
                    dfe.MutationType(distribution_type=t, distribution_args=p)

    def test_convert_to_substitution(self):
        mt = dfe.MutationType()
        assert mt.convert_to_substitution is True
        for c in (True, False):
            mt = dfe.MutationType(convert_to_substitution=c)
            assert mt.convert_to_substitution == c

    def test_dominance_coeff(self):
        mt = dfe.MutationType()
        assert mt.dominance_coeff == 0.5
        for dominance_coeff in (-10, 0, 0.5, 1, 50):
            mt = dfe.MutationType(dominance_coeff=dominance_coeff)
            assert mt.dominance_coeff == dominance_coeff

    def test_dominance_coeff_list(self):
        for dcl, dcb in (
            ([-0.1, 0.7, 1.2], [-2.1, 1.0]),
            ([-0.1, -0.7], [-2.1]),
        ):
            mt = dfe.MutationType(dominance_coeff_list=dcl, dominance_coeff_breaks=dcb)
            assert np.allclose(dcl, mt.dominance_coeff_list)
            assert np.allclose(dcb, mt.dominance_coeff_breaks)

    def test_pass_by_value(self):
        # make sure that for the arguments that are lists
        # we can't post-hoc modify them (and thus bypass validation)
        val = 0.5
        x = [val]
        mt = dfe.MutationType(distribution_args=x)
        x[0] = 2 * val + 1
        assert mt.distribution_args[0] == val
        x = [val, val]
        mt = dfe.MutationType(dominance_coeff_list=x, dominance_coeff_breaks=[0.0])
        x[0] = 2 * val + 1
        assert mt.dominance_coeff_list[0] == val
        x = [val]
        mt = dfe.MutationType(dominance_coeff_list=[0.0, 0.0], dominance_coeff_breaks=x)
        x[0] = 2 * val + 1
        assert mt.dominance_coeff_breaks[0] == val

    def test_bad_dominance_coeff(self):
        for dominance_coeff in (np.inf, np.nan, "abc", [], {}):
            with pytest.raises(ValueError, match="dominance.coeff"):
                dfe.MutationType(dominance_coeff=dominance_coeff)

    def test_bad_distribution_type(self):
        for distribution_type in (1, {}, None, "~", "!", "F"):
            with pytest.raises(ValueError):
                dfe.MutationType(distribution_type=distribution_type)

    def test_bad_dominance_coeff_list(self):
        dcl = [-0.1, 0.7, 1.2]
        dcb = [-2.1, 1.0]
        # can't specify both dominance_coeff and list
        with pytest.raises(ValueError, match="both dominance_coeff and"):
            dfe.MutationType(
                dominance_coeff=0.5,
                dominance_coeff_list=dcl,
                dominance_coeff_breaks=dcb,
            )
        with pytest.raises(ValueError, match="both dominance_coeff and"):
            dfe.MutationType(
                dominance_coeff=0.5,
                dominance_coeff_list=dcl,
            )
        with pytest.raises(ValueError, match="both dominance_coeff and"):
            dfe.MutationType(
                dominance_coeff=0.5,
                dominance_coeff_breaks=dcb,
            )
        # must have both coeffs and breaks
        with pytest.raises(ValueError, match="dominance.*no breaks"):
            dfe.MutationType(dominance_coeff_list=dcl)
        # must have at least 2 bins
        with pytest.raises(ValueError, match="dominance.*at least 2"):
            dfe.MutationType(
                dominance_coeff_list=[0.2],
                dominance_coeff_breaks=[],
            )
        # list must be one longer than breaks
        for x in ([], [0.0], dcl):
            with pytest.raises(ValueError, match="dominance.*equal to"):
                dfe.MutationType(
                    dominance_coeff_list=dcl,
                    dominance_coeff_breaks=x,
                )
        # bad coefficients
        for x in (np.inf, np.nan, "abc", [], {}):
            with pytest.raises(ValueError, match="dominance.coeff"):
                dfe.MutationType(
                    dominance_coeff_list=[x] + dcl[1:],
                    dominance_coeff_breaks=dcb,
                )
        # bad breaks
        for x in (np.inf, np.nan, "abc", [], {}):
            with pytest.raises(ValueError, match="dominance.*break"):
                dfe.MutationType(
                    dominance_coeff_list=dcl,
                    dominance_coeff_breaks=[x] + dcb[1:],
                )
        with pytest.raises(ValueError, match="nondecreasing"):
            dfe.MutationType(
                dominance_coeff_list=dcl,
                dominance_coeff_breaks=list(reversed(dcb)),
            )


class TestAllDFEs:
    """
    Tests for registered DFEs.
    """

    def test_non_empty(self):
        assert len(list(stdpopsim.all_dfes())) > 0

    @pytest.mark.parametrize("d", stdpopsim.all_dfes())
    def test_all_instances(self, d):
        assert isinstance(d, dfe.DFE)
        assert len(d.id) > 0
        assert len(d.description) > 0
        assert len(d.long_description) > 0
        assert len(d.citations) > 0
        assert len(d.mutation_types) > 0
        assert len(d.proportions) > 0


class TestDFEOutput:
    """
    Tests for DFE ouputs.
    """

    def test_str(self):
        d = dfe.DFE(
            id="xyz",
            description="abc",
            long_description="ABC",
            mutation_types=[dfe.MutationType(convert_to_substitution=True)],
            proportions=[1.0],
        )
        ds = str(d)
        assert "xyz" in ds
        assert "abc" in ds
        assert "ABC" in ds

    def test_wrap_long_lines(self):
        d = dfe.DFE(
            id="xyz",
            description="abc",
            long_description="ABC " * 10,
            mutation_types=[dfe.MutationType(convert_to_substitution=True)],
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
        d = dfe.DFE(
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
            mt = [dfe.MutationType() for _ in props]
            d = dfe.DFE(
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
        d = dfe.DFE(
            id="test",
            description="test",
            long_description="test",
            mutation_types=[dfe.MutationType(), dfe.MutationType()],
            proportions=[0.5, 0.5],
        )
        assert d.citations == []

    def test_create_dfe_without_mutation_types(self):
        d = dfe.DFE(
            id="test",
            description="test",
            long_description="test",
        )

        assert d.mutation_types == []

    def test_create_dfe_without_proportions(self):
        d = dfe.DFE(
            id="test",
            description="test",
            long_description="test",
            mutation_types=[dfe.MutationType()],
        )
        assert d.proportions == [1]

    def test_create_bad_dfes_message(self):
        with pytest.raises(TypeError, match="required .* arguments"):
            # id, description, long_description are required
            dfe.DFE()

        with pytest.raises(
            ValueError, match="proportions must be a list or numpy array."
        ):
            # proportions must be a list
            dfe.DFE(
                id="test", description="test", long_description="test", proportions=1
            )

        with pytest.raises(ValueError, match="mutation_types must be a list."):
            # mutation_types must be a list
            dfe.DFE(
                id="test",
                description="test",
                long_description="test",
                mutation_types=dfe.MutationType(),
                proportions=[0.5, 0.5],
            )

        with pytest.raises(
            ValueError,
            match="proportions and mutation_types must be lists of the same length.",
        ):
            # proportions and mutation_types must be the same length
            dfe.DFE(
                id="test",
                description="test",
                long_description="test",
                mutation_types=[dfe.MutationType(), dfe.MutationType()],
                proportions=[1.0],
            )

        with pytest.raises(
            ValueError, match="proportions must be nonnegative numbers."
        ):
            # proportions must be numbers
            dfe.DFE(
                id="test",
                description="test",
                long_description="test",
                mutation_types=[dfe.MutationType(), dfe.MutationType()],
                proportions=[True, "666"],
            )

        with pytest.raises(
            ValueError, match="proportions must be nonnegative numbers."
        ):
            # proportions must be positive
            dfe.DFE(
                id="test",
                description="test",
                long_description="test",
                mutation_types=[dfe.MutationType(), dfe.MutationType()],
                proportions=[-1, 0],
            )

        with pytest.raises(ValueError, match="proportions must sum to 1.0."):
            # proportions must sum 1.0
            dfe.DFE(
                id="test",
                description="test",
                long_description="test",
                mutation_types=[dfe.MutationType(), dfe.MutationType()],
                proportions=[1, 1],
            )

        with pytest.raises(
            ValueError, match="mutation_types must be a list of MutationType objects."
        ):
            # mutation_types must be a list of MutationType objects
            dfe.DFE(
                id="test",
                description="test",
                long_description="test",
                mutation_types=["neutral"],
            )

    def test_dfe_errors(self):
        m1 = stdpopsim.MutationType()
        m2 = stdpopsim.MutationType()
        with pytest.raises(ValueError, match="must be lists of the same length"):
            _ = stdpopsim.DFE(
                id="abc",
                description="test test",
                long_description="test test test test",
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
                _ = stdpopsim.DFE(
                    id="abc",
                    description="test test",
                    long_description="test test test test",
                    proportions=bad_props,
                    mutation_types=[m1, m2],
                )
        for bad_mut_types in ["abc", {}, [1.0, 2.0], [m1], m1, ["a", "b"]]:
            with pytest.raises(ValueError):
                _ = stdpopsim.DFE(
                    id="abc",
                    description="test test",
                    long_description="test test test test",
                    proportions=[0.6, 0.4],
                    mutation_types=bad_mut_types,
                )
        for bad_sums in [[-0.4, 0.5], [0.6, 0.8], [139487135987, 0.0], [0.2, 0.3]]:
            with pytest.raises(ValueError):
                _ = stdpopsim.DFE(
                    id="abc",
                    description="test test",
                    long_description="test test test test",
                    proportions=bad_sums,
                    mutation_types=[m1, m2],
                )

    def test_dfe_is_neutral(self):
        d = dfe.DFE(
            id="test",
            description="test",
            long_description="test",
            mutation_types=[dfe.MutationType(), dfe.MutationType()],
            proportions=[0.5, 0.5],
        )
        assert d.is_neutral is True

        d = dfe.DFE(
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

    @pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
    def test_no_msprime_dfe(self):
        # test we cannot simulate a non-neutral DFE with msprime
        m1 = dfe.MutationType(
            dominance_coeff=0.2,
            distribution_type="e",
            distribution_args=[0.1],
        )
        desc = "test test"
        long_desc = "test test 🐢"
        d = dfe.DFE(
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
        contig.add_dfe(
            intervals=np.array([[0, contig.length / 2]], dtype="int"),
            DFE=d,
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


class TestCreateNeutralDFE:
    """
    Tests for creating a neutral DFE instance.
    """

    def test_create_neutral_dfe(self):
        nd = dfe.neutral_dfe()
        assert isinstance(nd, dfe.DFE)
        assert nd.id == "neutral"
        assert nd.description == "neutral DFE"
        assert nd.long_description == "strictly neutral mutations"
        assert len(nd.mutation_types) == 1
        assert nd.mutation_types[0].is_neutral is True
        assert len(nd.proportions) == 1
        assert nd.proportions[0] == 1.0
        assert nd.is_neutral is True


class TestRegisterQCDFE:
    """
    Tests for registering a QC DFE.
    """

    def make_dfe(self, name):
        return dfe.DFE(
            id=name,
            description=name,
            long_description=name,
            mutation_types=[dfe.MutationType(convert_to_substitution=True)],
            proportions=[1.0],
        )

    def test_register_qc(self):
        dfe = self.make_dfe("test")
        dfe.register_qc(dfe)
        assert dfe.qc_dfe == dfe

    def test_already_registered(self):
        dfe = self.make_dfe("test")
        dfe.register_qc(dfe)
        with pytest.raises(ValueError) as e_info:
            dfe.register_qc(dfe)
        assert str(e_info.value) == "QC DFE already registered for test."

    def test_bad_qc_dfe(self):
        dfe = self.make_dfe("test")
        for not_a_dfe in [None, 15, "Zigzag_1S14"]:
            with pytest.raises(ValueError) as e_info:
                dfe.register_qc(not_a_dfe)
            assert (
                str(e_info.value) == f"Cannot register non-DFE '{not_a_dfe}' as QC DFE."
            )


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class DFETestMixin:
    """
    Mixin for testing specific DFEs. Subclass should extend
    this class and define the self.dfe (as the dfe instance).
    """

    dfe = None

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_simulation_runs(self):
        contig = stdpopsim.Contig.basic_contig(
            length=1_000_000,
            mutation_rate=1e-8,  # Ne=1e3 and length=1e6 so theta=40
            ploidy=2,
        )
        contig.clear_dfes()
        contig.add_dfe(
            intervals=np.array([[0, contig.length / 2]], dtype="int"),
            DFE=self.dfe,
        )

        model = stdpopsim.PiecewiseConstantSize(1000)
        samples = {"pop_0": 1}
        engine = stdpopsim.get_engine("slim")
        ts = engine.simulate(
            model, contig, samples, slim_scaling_factor=10, slim_burn_in=10, seed=42
        )

        mut_info = {}
        nonneutral = np.repeat(False, ts.num_mutations)
        for mut in ts.mutations():
            for j, md in zip(
                mut.derived_state.split(","), mut.metadata["mutation_list"]
            ):
                uid = f"{mut.id}_{j}"
                if md["selection_coeff"] != 0.0:
                    nonneutral[mut.id] = True
                if uid not in mut_info:
                    mut_info[uid] = md

        num_nonneutral = sum(nonneutral)
        nonneutral_positions = ts.tables.sites.position[
            ts.tables.mutations.site[nonneutral]
        ]
        assert np.all(nonneutral_positions <= ts.sequence_length / 2)
        assert len(mut_info.keys()) > 0  # number of mutations
        assert num_nonneutral > 0  # nonneutral mutations


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class CatalogDFETestMixin(DFETestMixin):
    """
    Mixin for DFEs in the catalog.
    """

    def test_id_valid(self):
        assert utils.is_valid_dfe_id(self.dfe.id)


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class QcdCatalogDFETestMixin(CatalogDFETestMixin):
    """
    Extends the tests to also check that the qc DFE is equal to
    the production DFE.
    """

    def test_mutation_types_match(self):
        mt1 = self.dfe.mutation_types
        mt2 = self.dfe.qc_dfe.mutation_types
        assert len(mt1) == len(mt2)

        for i in range(len(mt1)):
            assert mt1[i].dominance_coeff == mt2[i].dominance_coeff
            assert mt1[i].distribution_type == mt2[i].distribution_type
            assert np.allclose(mt1[i].distribution_args, mt2[i].distribution_args)
            assert mt1[i].convert_to_substitution == mt2[i].convert_to_substitution

    def test_proporitions_match(self):
        p1 = self.dfe.proportions
        p2 = self.dfe.qc_dfe.proportions
        assert np.allclose(p1, p2)


qc_test_classes = []
for species in stdpopsim.all_species():
    for d in species.dfes:
        superclasses = []
        if d.qc_dfe is not None:
            superclasses.append(QcdCatalogDFETestMixin)
        else:
            superclasses.append(CatalogDFETestMixin)
        classname = f"Test{species.id}{d.id}"
        cls = type(classname, tuple(superclasses), dict(dfe=d))
        qc_test_classes.append(cls)
# Basic sanity checks to double check that no errors get introduced
# that lead to these qc tests being skipped silently.
assert len(qc_test_classes) > 0
for cls in qc_test_classes:
    assert issubclass(cls, DFETestMixin)
    # Insert the class into the current test module's namespace.
    setattr(sys.modules[__name__], cls.__name__, cls)
