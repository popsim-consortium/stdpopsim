"""
Tests for DFEs in the catalog. Tests for dfe machinery are in
tests/test_traits.py
"""

import sys
import pytest
import pyslim
import stdpopsim
import numpy as np
from stdpopsim import dfe
from stdpopsim import utils


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
        dfe_two = self.make_dfe("test")
        dfe.register_qc(dfe_two)
        assert dfe.qc_dfe == dfe_two

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
        contig.clear_dmes()
        contig.add_dme(
            intervals=np.array([[0, contig.length / 2]], dtype="int"),
            DME=self.dfe,
        )

        model = stdpopsim.PiecewiseConstantSize(1000)
        samples = {"pop_0": 1}
        engine = stdpopsim.get_engine("slim")
        ts = engine.simulate(
            model,
            contig,
            samples,
            slim_scaling_factor=10,
            slim_burn_in=10,
            seed=42,
        )

        mut_info = {}
        nonneutral = np.repeat(False, ts.num_mutations)

        mut_metadata = pyslim.mutation_metadata(ts)
        assert ts.metadata["SLiM"]["traits"][0]["name"] == "fitnessT"
        for mut in ts.mutations():
            for j in mut.metadata["derived_states"]:
                md = mut_metadata[int(j)]
                uid = f"{mut.id}_{j}"
                sel_coeff = md["per_trait"][0]["effect_size"]
                if sel_coeff != 0.0:
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


class CatalogDFETestMixin(DFETestMixin):
    """
    Mixin for DFEs in the catalog.
    """

    def test_id_valid(self):
        assert utils.is_valid_dfe_id(self.dfe.id)


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
            if all(isinstance(x, str) for x in mt1[i].distribution_args):
                assert mt1[i].distribution_args == mt2[i].distribution_args
            else:
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
