"""
Tests for SLiM simulation engine.
"""
import os
import sys
import unittest
import tempfile

import tskit

import stdpopsim
import stdpopsim.cli
from . test_cli import capture_output


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

    # TODO: remove decorator once recombination maps are implemented.
    @unittest.expectedFailure
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

    @unittest.skipIf(sys.platform.startswith("win"), "no conda slim package for windows")
    def test_simulate(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("AraTha")
        contig = species.get_contig("chr5", length_multiplier=0.0001)
        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        model.generation_time = species.generation_time
        samples = model.get_samples(10)
        ts = engine.simulate(
                demographic_model=model, contig=contig, samples=samples)
        self.assertEqual(ts.num_samples, 10)


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

    @unittest.skipIf(sys.platform.startswith("win"), "no conda slim package for windows")
    def test_simulate(self):
        with tempfile.NamedTemporaryFile(mode="w") as f:
            self.docmd(f"HomSap -o {f.name}")
            ts = tskit.load(f.name)
        self.assertEqual(ts.num_samples, 10)

        saved_slim_env = os.environ.get("SLIM")

        with tempfile.NamedTemporaryFile(mode="w") as f:
            self.docmd(f"--slim-no-burnin --slim-path slim HomSap -o {f.name}")
            ts = tskit.load(f.name)
        self.assertEqual(ts.num_samples, 10)

        if saved_slim_env is None:
            del os.environ["SLIM"]
        else:
            os.environ["SLIM"] = saved_slim_env

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
