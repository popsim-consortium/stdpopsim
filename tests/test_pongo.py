import unittest
import io

import msprime

from stdpopsim import pongo


class TestPongoIM(unittest.TestCase):
    """
    Basic tests for the LockeEtAlPongoIM model.
    """

    def test_simulation_runs(self):
        model = pongo.LockeEtAlPongoIM()
        ts = msprime.simulate(
            samples=[msprime.Sample(pop, 0) for pop in range(2)],
            **model.asdict())
        self.assertEqual(ts.num_populations, 2)

    def test_debug_runs(self):
        model = pongo.LockeEtAlPongoIM()
        output = io.StringIO()
        model.debug(output)
        s = output.getvalue()
        self.assertGreater(len(s), 0)


class TestGenericConstantSize(unittest.TestCase):
    """
    Basic tests for the GenericConstantSize model.
    """

    def test_simulation_runs(self):
        model = pongo.GenericConstantSize()
        ts = msprime.simulate(
            samples=[msprime.Sample(0, 0), msprime.Sample(0, 0)],
            **model.asdict())
        self.assertEqual(ts.num_populations, 1)

    def test_debug_runs(self):
        model = pongo.GenericConstantSize()
        output = io.StringIO()
        model.debug(output)
        s = output.getvalue()
        self.assertGreater(len(s), 0)


class TestGenericTwoEpoch(unittest.TestCase):
    """
    Basic tests for the GenericTwoEpoch model.
    """

    def test_simulation_runs(self):
        model = pongo.GenericTwoEpoch()
        ts = msprime.simulate(
            samples=[msprime.Sample(0, 0), msprime.Sample(0, 0)],
            **model.asdict())
        self.assertEqual(ts.num_populations, 1)

    def test_debug_runs(self):
        model = pongo.GenericTwoEpoch()
        output = io.StringIO()
        model.debug(output)
        s = output.getvalue()
        self.assertGreater(len(s), 0)

    def test_debug_runs_v2(self):
        model = pongo.GenericTwoEpoch(100, 4)
        output = io.StringIO()
        model.debug(output)
        s = output.getvalue()
        self.assertGreater(len(s), 0)
