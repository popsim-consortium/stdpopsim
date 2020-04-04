"""
Tests for SLiM simulation engine.
"""
import os
import re
import sys
import unittest
import tempfile
import subprocess
from unittest import mock

import tskit
import pyslim

import stdpopsim
import stdpopsim.cli
from . test_cli import capture_output

IS_WINDOWS = sys.platform.startswith("win")


@unittest.skipIf(IS_WINDOWS, "SLiM not available on windows")
class TestAPI(unittest.TestCase):

    def test_bad_params(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr1")
        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        samples = model.get_samples(10)
        model.generation_time = species.generation_time

        for scaling_factor in (0, -1, -1e-6):
            with self.assertRaises(ValueError):
                engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_scaling_factor=scaling_factor,
                    slim_script=True)

        for burn_in in (-1, -1e-6):
            with self.assertRaises(ValueError):
                engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_burn_in=burn_in,
                    slim_script=True)

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
                demographic_model=model, contig=contig, samples=samples,
                slim_scaling_factor=10, slim_burn_in=0)
        self.assertEqual(ts.num_samples, 10)
        self.assertTrue(all(tree.num_roots == 1 for tree in ts.trees()))

    def test_recap_and_rescale(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr22", length_multiplier=0.001)
        model = species.get_demographic_model("OutOfAfrica_3G09")
        samples = model.get_samples(10, 10, 10)

        seed = 12
        ts1 = engine.simulate(
                demographic_model=model, contig=contig, samples=samples,
                slim_scaling_factor=10, slim_burn_in=0, seed=seed)

        out, _ = capture_output(
                engine.simulate,
                demographic_model=model, contig=contig, samples=samples,
                slim_script=True, slim_scaling_factor=10, slim_burn_in=0,
                seed=seed)

        match = re.search(r'"trees_file",\s*"([^"]*)"', out)
        self.assertIsNotNone(match)
        tmp_trees_file = match.group(1)

        slim_cmd = [engine.slim_path(), "-s", str(seed)]

        with tempfile.NamedTemporaryFile(mode="w") as slim_script, \
                tempfile.NamedTemporaryFile(mode="w") as trees_file:
            out = out.replace(tmp_trees_file, trees_file.name)
            slim_script.write(out)
            slim_script.flush()
            slim_cmd.append(slim_script.name)
            subprocess.check_call(slim_cmd, stdout=subprocess.DEVNULL)
            ts2_headless = pyslim.load(trees_file.name)

        ts2 = engine.recap_and_rescale(
                ts2_headless,
                demographic_model=model, contig=contig, samples=samples,
                slim_scaling_factor=10, seed=seed)

        tables1 = ts1.dump_tables()
        tables2 = ts2.dump_tables()

        self.assertEqual(tables1.nodes, tables2.nodes)
        self.assertEqual(tables1.edges, tables2.edges)
        self.assertEqual(tables1.mutations, tables2.mutations)


@unittest.skipIf(IS_WINDOWS, "SLiM not available on windows")
class TestCLI(unittest.TestCase):

    def docmd(self, _cmd):
        cmd = ("-e slim --slim-scaling-factor 20 --slim-burn-in 0 "
               f"{_cmd} -l 0.001 -c chr1 -s 1234 -q 10").split()
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
            self.docmd(f"HomSap -o {f.name}")
            ts = tskit.load(f.name)
        self.assertEqual(ts.num_samples, 10)

        # verify sample counts for a multipopulation demographic model
        with tempfile.NamedTemporaryFile(mode="w") as f:
            cmd = ("-e slim --slim-scaling-factor 20 --slim-burn-in 0 "
                   f"HomSap -o {f.name} -l 0.001 -c chr1 -s 1234 -q "
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
