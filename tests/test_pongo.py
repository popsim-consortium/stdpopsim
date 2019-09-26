import unittest
import io

import msprime

from stdpopsim import pongo


# TODO this needs to be brought up to date with the rest of the API when we
# do the QC process.


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
