import unittest
import io


import stdpopsim
from stdpopsim import pongo


class TestPongoIM(unittest.TestCase):
    """
    Basic tests for the LockeEtAlPongoIM model.
    """

    def test_simulation_runs(self):
        model = pongo.LockeEtAlPongoIM()
        contig = stdpopsim.Contig()
        samples = model.get_samples(2)
        engine = stdpopsim.get_default_engine()
        ts = engine.simulate(model, contig, samples)
        self.assertEqual(ts.num_populations, 2)

    def test_debug_runs(self):
        model = pongo.LockeEtAlPongoIM()
        output = io.StringIO()
        model.debug(output)
        s = output.getvalue()
        self.assertGreater(len(s), 0)
