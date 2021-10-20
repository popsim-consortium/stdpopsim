"""
Tests for classes that hold information about genomic region to be simulated.
"""
import numpy as np
import pytest
import sys
import stdpopsim
from stdpopsim import utils

IS_WINDOWS = sys.platform.startswith("win")


class TestGenomicElementType:
    def test_bad_proportions(self):
        # no proportion but mut type, prop > 1, sum(prop) > 1
        proportions = ([], [2], [0.5, 0.5, 0.5])
        mut_type_ids = ([0], [0], [0])
        for props in proportions:
            with pytest.raises(ValueError):
                stdpopsim.GenomicElementType(
                    intervals=np.array([[0, 1]]),
                    mutation_type_ids=mut_type_ids,
                    proportions=props,
                )


class TestContig:
    def test_all_intervals_array(self):
        contig = stdpopsim.Contig.basic_contig(length=100)
        assert len(contig.all_intervals_array) == 1

        truth = np.array([[20, 30, 1], [30, 50, 0]])
        contig.clear_genomic_mutation_types()
        contig.add_genomic_element_type(
            intervals=np.array([[30, 50, 0]]),
            mutation_types=[stdpopsim.ext.MutationType()],
            proportions=np.array([0]),
        )
        contig.add_genomic_element_type(
            intervals=np.array([[20, 30, 1]]),
            mutation_types=[stdpopsim.ext.MutationType()],
            proportions=np.array([0]),
        )
        assert (contig.all_intervals_array == truth).all()

    def test_add_mutation_type_errors(self):
        contig = stdpopsim.Contig.basic_contig(length=100)
        # diff length mutation types and proportions
        with pytest.raises(ValueError):
            contig.add_mutation_types([stdpopsim.ext.MutationType()], [], 0)
        # invalid ge type id
        with pytest.raises(ValueError):
            contig.add_mutation_types([stdpopsim.ext.MutationType()], [0], 1)
        contig.fully_neutral()
        # sum props > 1
        with pytest.raises(ValueError):
            contig.add_mutation_types([stdpopsim.ext.MutationType()], [1.1], 0)

    def test_add_genomic_element_type_errors(self):
        contig = stdpopsim.Contig.basic_contig(length=100)
        # bad intervals
        with pytest.raises(ValueError):
            contig.add_genomic_element_type(np.array([10, 20]), [], [])

    def test_add_genomic_element_type(self):
        contig = stdpopsim.Contig.basic_contig(length=100)
        contig.clear_genomic_mutation_types()
        assert contig.mutation_types == []
        assert contig.genomic_element_types == []
        intervals1 = np.array([[10, 20], [20, 40], [100, 200]])
        intervals2 = np.array([[5, 10], [45, 50], [200, 250]])
        m1 = stdpopsim.ext.MutationType()
        m2 = stdpopsim.ext.MutationType()
        contig.add_genomic_element_type(intervals1, [m1, m2], [0.5, 0.3])
        contig.add_genomic_element_type(intervals2, [m1], [0.2])
        assert len(contig.genomic_element_types) == 2
        assert len(contig.mutation_types) == 2

    def test_add_DFE_errors(self):
        contig = stdpopsim.Contig.basic_contig(length=100)
        props = [0.3, 0.7]
        mt = [stdpopsim.ext.MutationType() for _ in props]
        # bad intervals
        dfe = stdpopsim.DFE(
            id="abc",
            description="test",
            long_description="test test",
            proportions=props,
            mutation_types=mt,
        )
        with pytest.raises(ValueError):
            contig.add_DFE(np.array([10, 20]), dfe)

    def test_add_DFE(self):
        contig = stdpopsim.Contig.basic_contig(length=100)
        contig.clear_genomic_mutation_types()
        props = [0.3, 0.7]
        mt = [stdpopsim.ext.MutationType() for _ in props]
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
        contig.add_DFE(intervals=np.array([[10, 30], [50, 100]]), DFE=dfes[0])
        contig.add_DFE(intervals=np.array([[30, 40]]), DFE=dfes[1])
        assert len(contig.genomic_element_types) == 2
        assert len(contig.mutation_types) == len(mt)

    def test_add_DFE_neutral_background(self):
        contig = stdpopsim.Contig.basic_contig(length=100)
        contig.clear_genomic_mutation_types()
        props = [0.3, 0.7]
        mt = [stdpopsim.ext.MutationType() for _ in props]
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
        contig.add_DFE_neutral_background(
            intervals=np.array([[10, 30], [50, 70], [98, 99]]), DFE=dfes[0]
        )
        assert len(contig.genomic_element_types) == 2
        assert len(contig.mutation_types) == len(mt) + 1
        # also test case where last element ends at contig end
        contig.clear_genomic_mutation_types()
        contig.add_DFE_neutral_background(
            intervals=np.array([[10, 30], [50, 70], [98, 100]]), DFE=dfes[0]
        )
        assert len(contig.genomic_element_types) == 2
        assert len(contig.mutation_types) == len(mt) + 1

    def test_add_DFE_neutral_background_errors(self):
        contig = stdpopsim.Contig.basic_contig(length=100)
        contig.clear_genomic_mutation_types()
        props = [0.3, 0.7]
        mt = [stdpopsim.ext.MutationType() for _ in props]
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
        # add bad interval
        with pytest.raises(ValueError):
            contig.add_DFE_neutral_background(
                intervals=np.array([[10, 10]]), DFE=dfes[0]
            )


class TestAll_DFE_Models:
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
    def test_basic_DFE(self):
        desc = "test test"
        long_desc = "test test 🐢"
        for props in ([0.4, 0.6], [1.0, 0.0], [1.0], [1 / 3, 1 / 3, 1 / 3]):
            mt = [stdpopsim.ext.MutationType() for _ in props]
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

    def test_DFE_defaults(self):
        m1 = stdpopsim.ext.MutationType()
        desc = "test test"
        long_desc = "test test 🐢"
        dfe = stdpopsim.DFE(
            id="abc",
            description=desc,
            long_description=long_desc,
            mutation_types=[m1],
        )
        assert isinstance(dfe.proportions, list)
        assert len(dfe.proportions) == 1
        assert dfe.proportions[0] == 1.0

    def test_printing_DFE(self, capsys):
        m1 = stdpopsim.ext.MutationType()
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

    def test_DFE_errors(self):
        m1 = stdpopsim.ext.MutationType()
        m2 = stdpopsim.ext.MutationType()
        for bad_props in [["abc"], 1.0, [1.0], [0.2, 0.4, 0.4]]:
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

    @pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
    @pytest.mark.skip("this isn't tested for yet (issue #1016)")
    def test_no_msprime_DFE(self):
        # test we cannot simulate a non-neutral DFE with msprime
        m1 = stdpopsim.ext.MutationType(
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
        contig.clear_genomic_mutation_types()
        contig.add_DFE(
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

    @pytest.mark.filterwarnings("ignore:.*IncompletePopulationMetadataWarning.*")
    def test_simulation_runs(self):
        contig = stdpopsim.Contig.basic_contig(
            length=10000,
            mutation_rate=1e-6,  # Ne=1e3 and length=1e4 so theta=40
        )
        contig.clear_genomic_mutation_types()
        contig.add_DFE(
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
