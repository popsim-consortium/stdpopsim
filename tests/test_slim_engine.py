"""
Tests for SLiM simulation engine.
"""
import os
import sys
import unittest
import tempfile
from unittest import mock

import tskit

import stdpopsim
import stdpopsim.cli
from . test_cli import capture_output

IS_WINDOWS = sys.platform.startswith("win")


@unittest.skipIf(IS_WINDOWS, "SLiM not available on windows")
class TestAPI(unittest.TestCase):

    def test_script_generation(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr1")

        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        samples = model.get_samples(10)
        model.generation_time = species.generation_time
        out, _ = capture_output(
                engine.simulate,
                demographic_model=model, contig=contig, samples=samples,
                slim_script=True)
        self.assertTrue("sim.registerLateEvent" in out)

        model = species.get_demographic_model("AncientEurasia_9K19")
        samples = model.get_samples(1, 2, 3, 4, 5, 6, 7)
        out, _ = capture_output(
                engine.simulate,
                demographic_model=model, contig=contig, samples=samples,
                slim_script=True)
        self.assertTrue("sim.registerLateEvent" in out)

        model = species.get_demographic_model("AmericanAdmixture_4B11")
        samples = model.get_samples(10, 10, 10)
        out, _ = capture_output(
                engine.simulate,
                demographic_model=model, contig=contig, samples=samples,
                slim_script=True)
        self.assertTrue("sim.registerLateEvent" in out)

    def test_recombination_map(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr1", genetic_map="HapMapII_GRCh37")
        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        samples = model.get_samples(10)
        model.generation_time = species.generation_time
        out, _ = capture_output(
                engine.simulate,
                demographic_model=model, contig=contig, samples=samples,
                slim_script=True)
        self.assertTrue("sim.registerLateEvent" in out)

    def test_simulate(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("AraTha")
        contig = species.get_contig("chr5", length_multiplier=0.001)
        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        model.generation_time = species.generation_time
        samples = model.get_samples(10)
        ts = engine.simulate(
                demographic_model=model, contig=contig, samples=samples)
        self.assertEqual(ts.num_samples, 10)
        self.assertTrue(all(tree.num_roots == 1 for tree in ts.trees()))


@unittest.skipIf(IS_WINDOWS, "SLiM not available on windows")
class TestCLI(unittest.TestCase):

    def docmd(self, _cmd):
        cmd = f"-e slim {_cmd} -l 0.00001 -c chr1 -s 1234 -q 10".split()
        return capture_output(stdpopsim.cli.stdpopsim_main, cmd)

    def test_script_generation(self):
        out, _ = self.docmd("--slim-script HomSap")
        self.assertTrue("sim.registerLateEvent" in out)

        # msprime.MassMigration demographic events, with proportion<1.0
        # low level migration
        out, _ = self.docmd("--slim-script HomSap -d AncientEurasia_9K19")
        self.assertTrue("sim.registerLateEvent" in out)
        # simultaneous mass migrations, with proportions summing to 1.0
        out, _ = self.docmd("--slim-script HomSap -d AmericanAdmixture_4B11")
        self.assertTrue("sim.registerLateEvent" in out)

    def test_simulate(self):
        saved_slim_env = os.environ.get("SLIM")
        with tempfile.NamedTemporaryFile(mode="w") as f:
            self.docmd(f"--slim-path slim HomSap -o {f.name}")
            ts = tskit.load(f.name)
        self.assertEqual(ts.num_samples, 10)
        self.assertTrue(all(tree.num_roots == 1 for tree in ts.trees()))

        if saved_slim_env is None:
            del os.environ["SLIM"]
        else:
            os.environ["SLIM"] = saved_slim_env

        with tempfile.NamedTemporaryFile(mode="w") as f:
            self.docmd(f"--slim-no-recapitation HomSap -o {f.name}")
            ts = tskit.load(f.name)
        self.assertEqual(ts.num_samples, 10)

        with tempfile.NamedTemporaryFile(mode="w") as f:
            self.docmd(f"--slim-no-recapitation --slim-no-burnin HomSap -o {f.name}")
            ts = tskit.load(f.name)
        self.assertEqual(ts.num_samples, 10)

        # verify sample counts for a multipopulation demographic model
        with tempfile.NamedTemporaryFile(mode="w") as f:
            cmd = (f"-e slim HomSap -o {f.name} -l 0.00001 -c chr1 -s 1234 -q "
                   "-d OutOfAfrica_3G09 0 0 8").split()
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)
            ts = tskit.load(f.name)
        self.assertEqual(ts.num_populations, 3)
        observed_counts = [0, 0, 0]
        for sample in ts.samples():
            observed_counts[ts.get_population(sample)] += 1
        self.assertEqual(observed_counts[0], 0)
        self.assertEqual(observed_counts[1], 0)
        self.assertEqual(observed_counts[2], 8)
        self.assertTrue(all(tree.num_roots == 1 for tree in ts.trees()))

    def test_bad_slim_environ_var(self):
        saved_slim_env = os.environ.get("SLIM")

        os.environ["SLIM"] = "nonexistent"
        with self.assertRaises(FileNotFoundError):
            self.docmd("HomSap")

        if saved_slim_env is None:
            del os.environ["SLIM"]
        else:
            os.environ["SLIM"] = saved_slim_env

    def test_bad_slim_path(self):
        saved_slim_env = os.environ.get("SLIM")

        with self.assertRaises(FileNotFoundError):
            self.docmd("--slim-path nonexistent HomSap")

        if saved_slim_env is None:
            del os.environ["SLIM"]
        else:
            os.environ["SLIM"] = saved_slim_env


class TestSlimAvailable(unittest.TestCase):
    """
    Checks whether SLiM is available or not on platforms that support it.
    """
    def test_parser_has_options(self):
        parser = stdpopsim.cli.stdpopsim_cli_parser()
        with mock.patch("sys.exit"):
            _, stderr = capture_output(parser.parse_args, ["--help"])
            # On windows we should have no "slim" options
            self.assertEqual(IS_WINDOWS, "slim" not in stderr)

    def test_engine_available(self):
        all_engines = [engine.id for engine in stdpopsim.all_engines()]
        self.assertEqual(IS_WINDOWS, "slim" not in all_engines)
