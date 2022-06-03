"""
Tests for classes that hold information about genomic region to be simulated.
"""
import numpy as np
import pytest
import stdpopsim


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
