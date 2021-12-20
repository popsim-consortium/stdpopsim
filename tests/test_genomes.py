"""
Tests for classes that hold information about genomic region to be simulated.
"""
import numpy as np
import pytest
import sys
import stdpopsim
from stdpopsim import utils

IS_WINDOWS = sys.platform.startswith("win")


class TestContig(object):

    example_dfe = stdpopsim.DFE(
        id="abc",
        description="example DFE",
        long_description="test test beep boop beep",
        proportions=[0.3, 0.7],
        mutation_types=[stdpopsim.MutationType() for _ in range(2)],
    )

    def verify_dfe_breakpoints(self, contig):
        breaks, dfe_labels = contig.dfe_breakpoints()
        for j, intervals in enumerate(contig.interval_list):
            for left, right in intervals:
                assert left in breaks
                assert right in breaks
                k = np.searchsorted(breaks, left)
                assert dfe_labels[k] == j

    def test_default_dfe(self):
        contig = stdpopsim.Contig.basic_contig(length=100)
        assert len(contig.dfe_list) == 1
        assert len(contig.interval_list) == 1
        assert contig.dfe_list[0].id == "neutral"
        assert np.all(contig.interval_list[0] == np.array([[0, contig.length]]))

    def test_add_dfe_errors(self):
        contig = stdpopsim.Contig.basic_contig(length=100)
        # bad intervals
        with pytest.raises(ValueError):
            contig.add_dfe(np.array([10, 20]), self.example_dfe)
        with pytest.raises(ValueError):
            contig.add_dfe("abc", self.example_dfe)

    def test_dfe_breakpoints(self):
        contig = stdpopsim.Contig.basic_contig(length=100)
        contig.clear_dfes()
        mt = stdpopsim.MutationType()
        for j in range(3):
            dfe = stdpopsim.DFE(
                id=str(j),
                description="test",
                long_description="test test",
                mutation_types=[mt],
            )
            contig.add_dfe(
                np.array(
                    [[(j + 1) * 5, (j + 1) * 10], [(j + 1) * 20, (j + 1) * 20 + 1]],
                    dtype="int",
                ),
                dfe,
            )
        self.verify_dfe_breakpoints(contig)

    def test_add_dfe(self):
        for clear in (True, False):
            contig = stdpopsim.Contig.basic_contig(length=100)
            if clear:
                contig.clear_dfes()
            props = [0.3, 0.7]
            mt = [stdpopsim.MutationType() for _ in props]
            dfes = [
                stdpopsim.DFE(
                    id=str(j),
                    description="test",
                    long_description="test test",
                    proportions=props,
                    mutation_types=mt,
                )
                for j in range(3)
            ]
            contig.add_dfe(
                intervals=np.array([[10, 30], [50, 100]]),
                DFE=dfes[0],
            )
            contig.add_dfe(intervals=np.array([[30, 40]]), DFE=dfes[1])
            contig.add_dfe(intervals=np.array([[20, 60]]), DFE=dfes[2])
            assert len(contig.dfe_list) == 4 - clear
            assert len(contig.mutation_types()) == 7 - clear
            if clear:
                dfe_ids = []
                true_ints = []
            else:
                true_ints = [np.array([[0, 10]])]
                dfe_ids = ["neutral"]
            dfe_ids += [dfe.id for dfe in dfes]
            true_ints += [
                np.array([[10, 20], [60, 100]]),
                np.empty((0, 2)),
                np.array([[20, 60]]),
            ]
            for d, i in zip(contig.dfe_list, dfe_ids):
                assert d.id == i
            for a1, a2 in zip(contig.interval_list, true_ints):
                assert np.all(a1.shape == a2.shape)
                assert np.all(a1 == a2)
            self.verify_dfe_breakpoints(contig)

    @pytest.mark.skip(reason="TODO allow more flexible inputs")
    def test_add_dfe_interval_formats(self):
        L = 50818
        for intervals in (
            [[0, L]],
            [[0, int(0.2 * L)]],
            ([0, int(0.2 * L)], [int(0.6 * L), L]),
        ):
            contig = stdpopsim.Contig.basic_contig(length=L)
            contig.add_dfe(intervals=intervals, DFE=self.example_dfe)
            np.testing.assert_array_equal(intervals, contig.interval_list[1])

    def test_is_neutral(self):
        for neutral in (True, False):
            for dist in ("f", "e"):
                contig = stdpopsim.Contig.basic_contig(length=100)
                contig.clear_dfes()
                props = [0.3, 0.7]
                if neutral:
                    s = 0
                else:
                    s = 0.1
                mt = [
                    stdpopsim.MutationType(
                        distribution_type=dist, distribution_args=[s]
                    )
                    for _ in props
                ]
                dfes = [
                    stdpopsim.DFE(
                        id=str(j),
                        description="test",
                        long_description="test test",
                        proportions=props,
                        mutation_types=mt,
                    )
                    for j in range(2)
                ]
                contig.add_dfe(
                    intervals=np.array([[10, 30], [50, 100]]),
                    DFE=dfes[0],
                )
                contig.add_dfe(intervals=np.array([[30, 40]]), DFE=dfes[1])
                # exponential with mean zero doesn't count as neutral!
                assert contig.is_neutral is (neutral and dist == "f")


class TestMutationTypes(object):
    def test_netural_mutation_type(self):
        mt = stdpopsim.MutationType()
        assert mt.is_neutral

    def test_mutation_types(self):
        mut_params = {
            "f": ([-0.1], [0], [0.1], [50]),
            "g": ([-0.1, 0.1], [0.1, 0.1], [50, 50]),
            "e": ([0.1], [10], [5000], [0]),
            "n": ([-0.1, 0.2], [0.1, 0.1], [50, 50]),
            "w": ([0.1, 0.2], [0.1, 0.1], [50, 50]),
            "l": ([-0.1, 0.2], [0.1, 0.1], [50, 50]),
        }
        for t in mut_params:
            for p in mut_params[t]:
                mt = stdpopsim.MutationType(distribution_type=t, distribution_args=p)
                if t == "l":
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
            "n": ([], [0.1, -1], [0.1, 0.4, 0.5], [0.1], [0.1, 0.0], [0.3, np.inf]),
            "w": ([], [-0.1, 1], [0.1, -1], [0.1, 0.4, 0.5], [0.1], [np.inf, 2.3]),
            "l": ([], [0.1, -1], [0.1, 0.4, 0.5], [0.1], [0.1, np.inf]),
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

    def test_bad_dominance_coeff(self):
        for dominance_coeff in (np.inf, np.nan):
            with pytest.raises(ValueError):
                stdpopsim.MutationType(dominance_coeff=dominance_coeff)

    def test_bad_distribution_type(self):
        for distribution_type in (1, {}, None, "~", "!", "F"):
            with pytest.raises(ValueError):
                stdpopsim.MutationType(distribution_type=distribution_type)

    def test_not_neutral_mutation_type(self):
        mt = stdpopsim.MutationType(
            distribution_type="f",
            distribution_args=[1],
        )
        assert not mt.is_neutral


class TestAll_dfe_Models:
    """
    Tests for registered simulation dfe models.
    """

    def test_non_empty(self):
        assert len(list(stdpopsim.all_dfes())) > 0

    @pytest.mark.parametrize("DFE", stdpopsim.all_dfes())
    def test_all_instances(self, DFE):
        assert len(DFE.id) > 0
        assert len(DFE.description) > 0
        assert len(DFE.long_description) > 0
        assert len(DFE.citations) > 0


class TestDFE:
    def test_basic_dfe(self):
        desc = "test test"
        long_desc = "test test 🐢"
        for props in ([0.4, 0.6], [1.0, 0.0], [1.0], [1 / 3, 1 / 3, 1 / 3]):
            mt = [stdpopsim.MutationType() for _ in props]
            dfe = stdpopsim.DFE(
                id="abc",
                description=desc,
                long_description=long_desc,
                proportions=props,
                mutation_types=mt,
            )
            assert dfe.id == "abc"
            assert dfe.description == desc
            assert dfe.long_description == long_desc
            for a, b in zip(props, dfe.proportions):
                assert a == b
            for a, b in zip(mt, dfe.mutation_types):
                assert a == b
        assert dfe.is_neutral

    def test_dfe_defaults(self):
        m1 = stdpopsim.MutationType()
        desc = "test test"
        long_desc = "test test 🐢"
        dfe = stdpopsim.DFE(
            id="abc",
            description=desc,
            long_description=long_desc,
            mutation_types=[m1],
            proportions=[],
        )
        assert isinstance(dfe.proportions, list)
        assert len(dfe.proportions) == 1
        assert dfe.proportions[0] == 1
        assert dfe.is_neutral

    def test_dfe_is_neutral(self):
        for neutral in (True, False):
            for dist in ("f", "e"):
                props = [0.3, 0.7]
                if neutral:
                    s = 0
                else:
                    s = 0.1
                mt = [
                    stdpopsim.MutationType(
                        distribution_type=dist, distribution_args=[s]
                    )
                    for _ in props
                ]
                dfe = stdpopsim.DFE(
                    id=0,
                    description="test",
                    long_description="test test",
                    proportions=props,
                    mutation_types=mt,
                )
                assert dfe.is_neutral is (neutral and dist == "f")

    @pytest.mark.usefixtures("capsys")
    def test_printing_dfe(self, capsys):
        m1 = stdpopsim.MutationType()
        desc = "test test"
        long_desc = "test test 🐢"
        dfe = stdpopsim.DFE(
            id="abc",
            description=desc,
            long_description=long_desc,
            mutation_types=[m1],
        )
        print(dfe)
        captured = capsys.readouterr()
        assert "DFE:\n" in captured.out
        assert "║" in captured.out
        assert "║  id               = abc\n" in captured.out

    def test_dfe_errors(self):
        m1 = stdpopsim.MutationType()
        m2 = stdpopsim.MutationType()
        for bad_props in [["abc"], 1.0, [1.0], [0.2, 0.4, 0.4], [-0.1, -0.1]]:
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
            print(bad_sums)
            with pytest.raises(ValueError):
                _ = stdpopsim.DFE(
                    id="abc",
                    description="test test",
                    long_description="test test test test",
                    proportions=bad_sums,
                    mutation_types=[m1, m2],
                )

    @pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
    def test_no_msprime_dfe(self):
        # test we cannot simulate a non-neutral DFE with msprime
        m1 = stdpopsim.MutationType(
            dominance_coeff=0.2,
            distribution_type="e",
            distribution_args=[0.1],
        )
        desc = "test test"
        long_desc = "test test 🐢"
        dfe = stdpopsim.DFE(
            id="abc",
            description=desc,
            long_description=long_desc,
            mutation_types=[m1],
        )
        contig = stdpopsim.Contig.basic_contig(
            length=10000,
            mutation_rate=1e-6,
        )
        contig.add_dfe(
            intervals=np.array([[0, contig.length / 2]], dtype="int"),
            DFE=dfe,
        )
        model = stdpopsim.PiecewiseConstantSize(1000)
        samples = model.get_samples(2)
        engine = stdpopsim.get_engine("msprime")
        with pytest.raises(ValueError):
            _ = engine.simulate(
                model,
                contig,
                samples,
            )


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class DFEModelTestMixin:
    """
    Mixin for testing specific DFEs. Subclasses should extend
    this class and define the self.model (as the model instance).
    """

    # To be defined in subclasses.
    model = None

    # @pytest.mark.filterwarnings("ignore:.*IncompletePopulationMetadataWarning.*")
    @pytest.mark.skip("need to remove genomic_elements from slim-engine")
    def test_simulation_runs(self):
        contig = stdpopsim.Contig.basic_contig(
            length=10000,
            mutation_rate=1e-6,  # Ne=1e3 and length=1e4 so theta=40
        )
        contig.clear_dfes()
        contig.add_dfe(
            intervals=np.array([[0, contig.length / 2]], dtype="int"),
            DFE=self.dfe,
        )

        # Generate vector with 2 samples for each pop with sampling enabled
        sample_count = []
        for p in self.model.populations:
            if p.allow_samples:
                sample_count.append(2)
            else:
                sample_count.append(0)
        samples = self.model.get_samples(*sample_count)
        engine = stdpopsim.get_engine("slim")
        ts = engine.simulate(
            self.model, contig, samples, slim_scaling_factor=10, slim_burn_in=10
        )

        mut_info = {}
        nonneutral = np.repeat(False, ts.num_mutations)
        for mut in ts.mutations():
            for j, md in zip(
                mut.derived_state.split(","), mut.metadata["mutation_list"]
            ):
                if md["selection_coeff"] != 0.0:
                    nonneutral[mut.id] = True
                if j not in mut_info:
                    mut_info[int(j)] = md

        num_nonneutral = sum(nonneutral)
        nonneutral_positions = ts.tables.sites.position[
            ts.tables.mutations.site[nonneutral]
        ]
        assert np.all(nonneutral_positions <= ts.sequence_length / 2)
        assert ts.num_populations == self.model.num_populations
        assert len(mut_info.keys()) > 0  # number of mutations
        assert num_nonneutral > 0  # nonneutral mutations


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class CatalogDFEModelTestMixin(DFEModelTestMixin):
    """
    Mixin for demographic models in the catalog.
    """

    def test_id_valid(self):
        assert utils.is_valid_dfe_id(self.dfe.id)


qc_test_classes = []
for species in stdpopsim.all_species():
    for dfe in species.dfes:
        model = stdpopsim.PiecewiseConstantSize(1000)
        superclasses = []
        superclasses.append(CatalogDFEModelTestMixin)
        classname = f"Test{species.id}{model.id}{dfe.id}"
        cls = type(classname, tuple(superclasses), dict(model=model, dfe=dfe))
        qc_test_classes.append(cls)

# Basic sanity checks to double check that no errors get introduced
# that lead to these qc tests being skipped silently.
assert len(qc_test_classes) > 0
for cls in qc_test_classes:
    assert issubclass(cls, DFEModelTestMixin)
    # Insert the class into the current test module's namespace.
    setattr(sys.modules[__name__], cls.__name__, cls)
