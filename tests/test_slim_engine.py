"""
Tests for SLiM simulation engine.
"""
import os
import re
import io
import sys
import unittest
import tempfile
import math
from unittest import mock

import tskit
import pyslim
import msprime

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

        for scaling_factor in (0, -1, -1e-6):
            with self.assertRaises(ValueError):
                engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_scaling_factor=scaling_factor,
                    dry_run=True)

        for burn_in in (-1, -1e-6):
            with self.assertRaises(ValueError):
                engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_burn_in=burn_in,
                    dry_run=True)

    def test_script_generation(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr1")

        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        samples = model.get_samples(10)
        out, _ = capture_output(
                engine.simulate,
                demographic_model=model, contig=contig, samples=samples,
                slim_script=True)
        self.assertTrue("sim.registerLateEvent" in out)

        model = species.get_demographic_model("AncientEurasia_9K19")
        samples = model.get_samples(10, 20, 30, 40, 50, 60, 70)
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
        out, _ = capture_output(
                engine.simulate,
                demographic_model=model, contig=contig, samples=samples,
                dry_run=True)

    def test_simulate(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("AraTha")
        contig = species.get_contig("chr5", length_multiplier=0.001)
        model = stdpopsim.PiecewiseConstantSize(species.population_size)
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

        with tempfile.NamedTemporaryFile(mode="w") as slim_script, \
                tempfile.NamedTemporaryFile(mode="w") as trees_file:
            out = out.replace(tmp_trees_file, trees_file.name)
            slim_script.write(out)
            slim_script.flush()
            engine._run_slim(slim_script.name, seed=seed)
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
        cmd = ("-q -e slim --slim-scaling-factor 20 --slim-burn-in 0 "
               f"{_cmd} -l 0.001 -c chr1 -s 1234 10").split()
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
            cmd = ("-q -e slim --slim-scaling-factor 20 --slim-burn-in 0 "
                   f"HomSap -o {f.name} -l 0.001 -c chr1 -s 1234 "
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

    @mock.patch("stdpopsim.slim_engine._SLiMEngine.get_version", return_value="64")
    def test_dry_run(self, _mocked_get_version):
        # --dry-run should run slim, but not create an output file.
        with mock.patch("subprocess.Popen", autospec=True) as mocked_popen:
            # Popen is used as a context manager, so we frob the return value
            # of the returned context manager, rather than the Popen mock itself.
            proc = mocked_popen.return_value.__enter__.return_value
            proc.returncode = 0
            proc.stdout = io.StringIO()
            proc.stderr = io.StringIO()
            with tempfile.NamedTemporaryFile(mode="w") as f:
                self.docmd(f"HomSap --dry-run -o {f.name}")
        mocked_popen.assert_called_once()
        self.assertTrue("slim" in mocked_popen.call_args[0][0])
        with tempfile.NamedTemporaryFile(mode="w") as f:
            self.docmd(f"HomSap --dry-run -o {f.name}")
            self.assertEqual(os.stat(f.name).st_size, 0)

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


@unittest.skipIf(IS_WINDOWS, "SLiM not available on windows")
class TestWarningsAndErrors(unittest.TestCase):
    """
    Checks that warning messages are printed when appropriate.
    """
    def test_odd_sample_warning(self):
        cmd = "-q -e slim --slim-script HomSap -d OutOfAfrica_2T12 4 6".split()
        with mock.patch("warnings.warn", autospec=True) as mock_warning:
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)
        self.assertEqual(mock_warning.call_count, 0)

        cmd = "-q -e slim --slim-script HomSap -d OutOfAfrica_2T12 4 5".split()
        with mock.patch("warnings.warn", autospec=True) as mock_warning:
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)
        self.assertEqual(mock_warning.call_count, 1)

        cmd = "-q -e slim --slim-script HomSap -d OutOfAfrica_2T12 3 5".split()
        with mock.patch("warnings.warn", autospec=True) as mock_warning:
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)
        self.assertEqual(mock_warning.call_count, 2)

    def triplet(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr22", length_multiplier=0.001)
        return engine, species, contig

    def test_bad_population_size_addSubPop(self):
        engine, species, contig = self.triplet()
        model = stdpopsim.PiecewiseConstantSize(100)
        samples = model.get_samples(2)

        with mock.patch("warnings.warn", autospec=True) as mock_warning:
            engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_scaling_factor=10, dry_run=True)
        mock_warning.assert_called_once()

    def test_no_populations_in_generation_1(self):
        engine, species, contig = self.triplet()
        model = stdpopsim.PiecewiseConstantSize(100)
        samples = model.get_samples(2)

        with self.assertRaises(stdpopsim.SLiMException):
            engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_scaling_factor=200, dry_run=True)

    def test_bad_population_size_addSubpopSplit(self):
        engine, species, contig = self.triplet()
        model = stdpopsim.IsolationWithMigration(
                NA=1000, N1=100, N2=1000, T=1000, M12=0, M21=0)
        samples = model.get_samples(2)

        with mock.patch("warnings.warn", autospec=True) as mock_warning:
            engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_scaling_factor=10, dry_run=True)
        mock_warning.assert_called_once()

        with self.assertRaises(stdpopsim.SLiMException):
            engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_scaling_factor=200, dry_run=True)

    def test_bad_population_size_setSubpopulationSize(self):
        engine, species, contig = self.triplet()
        model = stdpopsim.PiecewiseConstantSize(100, (1000, 1000))
        samples = model.get_samples(2)

        with mock.patch("warnings.warn", autospec=True) as mock_warning:
            engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_scaling_factor=10, dry_run=True)
        mock_warning.assert_called_once()

        with self.assertRaises(stdpopsim.SLiMException):
            engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_scaling_factor=200, dry_run=True)

    def test_sample_size_too_big(self):
        engine, species, contig = self.triplet()
        model = stdpopsim.PiecewiseConstantSize(1000)
        samples = model.get_samples(300)

        with self.assertRaises(stdpopsim.SLiMException):
            engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_scaling_factor=10, dry_run=True)

    def exp_decline(self, N0=100, N1=1000, T=1000):
        """
        One population model with exponential decline in population size.
        Used for testing that growth rates are handled appropriately.
        """
        r = math.log(N0 / N1) / T
        return stdpopsim.DemographicModel(
                id="exp_decline",
                description="exp_decline",
                long_description="exp_decline",
                populations=[stdpopsim.models._pop0],
                generation_time=1,
                population_configurations=[
                    msprime.PopulationConfiguration(
                        initial_size=N0, growth_rate=r,
                        metadata=stdpopsim.models._pop0.asdict())
                    ],
                demographic_events=[
                    msprime.PopulationParametersChange(
                        time=T, initial_size=N1, growth_rate=0, population_id=0),
                    ],
                )

    def test_bad_population_size_exp_decline(self):
        engine, species, contig = self.triplet()
        model = self.exp_decline()
        samples = model.get_samples(2)

        with mock.patch("warnings.warn", autospec=True) as mock_warning:
            engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_scaling_factor=10, dry_run=True)
        mock_warning.assert_called_once()

        with self.assertRaises(stdpopsim.SLiMException):
            engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_scaling_factor=200, dry_run=False)

    def test_sample_size_too_big_exp_decline(self):
        engine, species, contig = self.triplet()
        model = self.exp_decline()
        samples = model.get_samples(30)

        with self.assertRaises(stdpopsim.SLiMException):
            engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_scaling_factor=10, dry_run=True)


class TestSlimAvailable(unittest.TestCase):
    """
    Checks whether SLiM is available or not on platforms that support it.
    """
    def test_parser_has_options(self):
        parser = stdpopsim.cli.stdpopsim_cli_parser()
        with mock.patch("sys.exit", autospec=True):
            _, stderr = capture_output(parser.parse_args, ["--help"])
            # On windows we should have no "slim" options
            self.assertEqual(IS_WINDOWS, "slim" not in stderr)

    def test_engine_available(self):
        all_engines = [engine.id for engine in stdpopsim.all_engines()]
        self.assertEqual(IS_WINDOWS, "slim" not in all_engines)
