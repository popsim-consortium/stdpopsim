"""
Test cases for the command line interfaces to stdpopsim
"""
import unittest
import tempfile
import pathlib
import subprocess
import json
import sys
import io
import argparse  # NOQA
from unittest import mock

import tskit
import msprime

import stdpopsim
import stdpopsim.cli as cli


class TestException(Exception):
    """
    Custom exception we can throw for testing.
    """


def capture_output(func, *args, **kwargs):
    """
    Runs the specified function and arguments, and returns the
    tuple (stdout, stderr) as strings.
    """
    buffer_class = io.BytesIO
    if sys.version_info[0] == 3:
        buffer_class = io.StringIO
    stdout = sys.stdout
    sys.stdout = buffer_class()
    stderr = sys.stderr
    sys.stderr = buffer_class()

    try:
        func(*args, **kwargs)
        stdout_output = sys.stdout.getvalue()
        stderr_output = sys.stderr.getvalue()
    finally:
        sys.stdout.close()
        sys.stdout = stdout
        sys.stderr.close()
        sys.stderr = stderr
    return stdout_output, stderr_output


class TestProvenance(unittest.TestCase):
    """
    Test basic provenance properties.
    """
    def test_schema_validates(self):
        d = cli.get_provenance_dict()
        tskit.validate_provenance(d)

    def test_environment(self):
        # Basic environment should be the same as tskit.
        d_stdpopsim = cli.get_provenance_dict()
        d_tskit = tskit.provenance.get_provenance_dict()
        self.assertEqual(d_stdpopsim["environment"]["os"], d_tskit["environment"]["os"])
        self.assertEqual(
            d_stdpopsim["environment"]["python"], d_tskit["environment"]["python"])

    def test_libraries(self):
        libs = cli.get_provenance_dict()["environment"]["libraries"]
        self.assertEqual(libs["tskit"]["version"], tskit.__version__)
        self.assertEqual(libs["msprime"]["version"], msprime.__version__)

    def test_software(self):
        software = cli.get_provenance_dict()["software"]
        self.assertEqual(
            software, {"name": "stdpopsim", "version": stdpopsim.__version__})

    def test_parameters(self):
        d = cli.get_provenance_dict()["parameters"]
        self.assertEqual(d["command"], sys.argv[0])
        self.assertEqual(d["args"], sys.argv[1:])


class TestHomoSapiensArgumentParser(unittest.TestCase):
    """
    Tests for the argument parsers.
    """

    def test_defaults(self):
        parser = cli.stdpopsim_cli_parser()
        cmd = "homsap"
        output = "test.trees"
        args = parser.parse_args([cmd, "2", output])
        self.assertEqual(args.output, output)
        self.assertEqual(args.seed, None)
        self.assertEqual(args.samples, [2])

    def test_seed(self):
        parser = cli.stdpopsim_cli_parser()
        cmd = "homsap"
        output = "test.trees"
        args = parser.parse_args([cmd, "2", "-s", "1234", output])
        self.assertEqual(args.output, output)
        self.assertEqual(args.samples, [2])
        self.assertEqual(args.seed, 1234)

        args = parser.parse_args([cmd, "2", "--seed", "14", output])
        self.assertEqual(args.output, output)
        self.assertEqual(args.samples, [2])
        self.assertEqual(args.seed, 14)


class TestEndToEnd(unittest.TestCase):
    """
    Checks that simulations we run from the CLI have plausible looking output.
    """
    def verify(self, cmd, num_samples, seed=1):
        with tempfile.TemporaryDirectory() as tmpdir:
            filename = pathlib.Path(tmpdir) / "output.trees"
            full_cmd = cmd + f" {filename} --seed={seed}"
            with mock.patch("stdpopsim.cli.setup_logging"):
                stdout, stderr = capture_output(cli.stdpopsim_main, full_cmd.split())
            self.assertEqual(len(stderr), 0)
            self.assertGreater(len(stdout), 0)
            ts = tskit.load(str(filename))
        self.assertEqual(ts.num_samples, num_samples)
        provenance = json.loads(ts.provenance(0).record)
        prov_seed = provenance["parameters"]["random_seed"]
        self.assertEqual(prov_seed, seed)

    def test_homsap_seed(self):
        cmd = "homsap -c chr22 -l0.1 20"
        self.verify(cmd, num_samples=20, seed=1234)

    def test_homsap_constant(self):
        cmd = "homsap -c chr22 -l0.1 20"
        self.verify(cmd, num_samples=20)

    def test_tennessen_two_pop_ooa(self):
        cmd = "homsap -c chr22 -l0.1 -m ooa_2 2 3"
        self.verify(cmd, num_samples=5)

    def test_gutenkunst_three_pop_ooa(self):
        cmd = "homsap -c chr1 -l0.01 -m ooa_3 10"
        self.verify(cmd, num_samples=10)

    def test_browning_america(self):
        cmd = "homsap -c chr1 -l0.01 -m america 10"
        self.verify(cmd, num_samples=10)

    def test_ragsdale_archaic(self):
        cmd = "homsap -c chr1 -l0.01 -m ooa_archaic 10"
        self.verify(cmd, num_samples=10)

    def test_schiffels_zigzag(self):
        cmd = "homsap -c chr1 -l0.01 -m zigzag 2"
        self.verify(cmd, num_samples=2)

    def test_dromel_constant(self):
        cmd = "dromel -c chr2L -l0.001 4"
        self.verify(cmd, num_samples=4)

    def test_li_stephan_two_population(self):
        cmd = "dromel -c chr2L -l0.001 -m ooa_2 3"
        self.verify(cmd, num_samples=3)

    def test_aratha_constant(self):
        cmd = "aratha -l 0.001 8"
        self.verify(cmd, num_samples=8)

    def test_durvusula_2017_msmc(self):
        cmd = "aratha -l 0.001 -m fixme 7"
        self.verify(cmd, num_samples=7)

    def test_lapierre_constant(self):
        cmd = "esccol -l 1e-7 2"
        self.verify(cmd, num_samples=2)


class TestEndToEndSubprocess(TestEndToEnd):
    """
    Run the commands in a subprocess so that we can verify the provenance is
    stored correctly.
    """
    def verify(self, cmd, num_samples, seed=1):
        with tempfile.TemporaryDirectory() as tmpdir:
            filename = pathlib.Path(tmpdir) / "output.trees"
            full_cmd = f"{sys.executable} -m stdpopsim {cmd} {filename} -s {seed} -q"
            subprocess.run(full_cmd, shell=True, check=True)
            ts = tskit.load(str(filename))
        self.assertEqual(ts.num_samples, num_samples)
        provenance = json.loads(ts.provenance(ts.num_provenances - 1).record)
        tskit.validate_provenance(provenance)
        stored_cmd = provenance["parameters"]["args"]
        self.assertEqual(stored_cmd[-1], "-q")
        self.assertEqual(stored_cmd[-2], str(seed))
        self.assertEqual(stored_cmd[-3], "-s")
        self.assertEqual(stored_cmd[:-4], cmd.split())
        provenance = json.loads(ts.provenance(0).record)
        prov_seed = provenance["parameters"]["random_seed"]
        self.assertEqual(prov_seed, seed)


class TestSetupLogging(unittest.TestCase):
    """
    Tests that setup logging has the desired effect.
    """
    basic_cmd = ["homsap", "10", "tmp.trees"]

    def test_default(self):
        parser = cli.stdpopsim_cli_parser()
        args = parser.parse_args(self.basic_cmd)
        with mock.patch("daiquiri.setup") as mocked_setup:
            cli.setup_logging(args)
            mocked_setup.assert_called_once_with(level="WARN")

    def test_verbose(self):
        parser = cli.stdpopsim_cli_parser()
        args = parser.parse_args(["-v"] + self.basic_cmd)
        with mock.patch("daiquiri.setup") as mocked_setup:
            cli.setup_logging(args)
            mocked_setup.assert_called_once_with(level="INFO")

    def test_very_verbose(self):
        parser = cli.stdpopsim_cli_parser()
        args = parser.parse_args(["-vv"] + self.basic_cmd)
        with mock.patch("daiquiri.setup") as mocked_setup:
            cli.setup_logging(args)
            mocked_setup.assert_called_once_with(level="DEBUG")


class TestErrors(unittest.TestCase):

    # Need to mock out setup_logging here or we spew logging to the console
    # in later tests.
    @mock.patch("stdpopsim.cli.setup_logging")
    def run_stdpopsim(self, command, mock_setup_logging):
        stdout, stderr = capture_output(cli.stdpopsim_main, command)
        self.assertEqual(stderr, "")
        self.assertEqual(stdout, "")
        self.assertTrue(mock_setup_logging.called)

    def test_exit(self):
        with mock.patch("sys.exit", side_effect=TestException) as mocked_exit:
            with self.assertRaises(TestException):
                cli.exit("XXX")
            mocked_exit.assert_called_once()
            args = mocked_exit.call_args[0]
            self.assertEqual(len(args), 1)
            self.assertTrue(args[0].endswith("XXX"))

    # Need to mock out setup_logging here or we spew logging to the console
    # in later tests.
    @mock.patch("stdpopsim.cli.setup_logging")
    def verify_bad_samples(self, cmd, mock_setup_logging):
        with mock.patch("stdpopsim.cli.exit", side_effect=TestException) as mocked_exit:
            with self.assertRaises(TestException):
                cli.stdpopsim_main(cmd.split())
            mocked_exit.assert_called_once()

    def test_default(self):
        self.verify_bad_samples("homsap -q 2 3 tmp.trees")

    def test_tennessen_model(self):
        self.verify_bad_samples("homsap  -q -m ooa_2 2 3 4 tmp.trees")

    def test_gutenkunst_three_pop_ooa(self):
        self.verify_bad_samples("homsap -q -m ooa_3 2 3 4 5 tmp.trees")

    def test_browning_america(self):
        self.verify_bad_samples("homsap -q -m america 2 3 4 5 6 tmp.trees")


class TestHelp(unittest.TestCase):

    def run_stdpopsim(self, command):
        with mock.patch(
                "argparse.ArgumentParser.exit",
                side_effect=TestException) as mocked_exit:
            with self.assertRaises(TestException):
                capture_output(cli.stdpopsim_main, command.split())
            mocked_exit.assert_called_once()

    def test_basic_help(self):
        self.run_stdpopsim("--help")

    def test_homsap_help(self):
        self.run_stdpopsim("homsap --help")

    def test_homsap_models_help(self):
        self.run_stdpopsim("homsap --help-models")
        self.run_stdpopsim("homsap --help-models ooa_3")

    def test_all_species_model_help(self):
        for species in stdpopsim.all_species():
            self.run_stdpopsim(f"{species} --help-models")


class TestWriteCitations(unittest.TestCase):
    """
    Make sure all models can write citation information.
    """
    def test_model(self):
        contig = stdpopsim.Contig()
        species = stdpopsim.get_species("homsap")
        model = species.get_model("ooa_3")
        stdout, stderr = capture_output(cli.write_citations, contig, model)
        self.assertEqual(len(stderr), 0)
        # TODO Parse out the output for the model and check that the text is
        # in there.
        self.assertGreater(len(stdout), 0)

    def test_genetic_map(self):
        species = stdpopsim.get_species("homsap")
        contig = species.get_contig("chr22", genetic_map="HapmapII_GRCh37")
        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        stdout, stderr = capture_output(cli.write_citations, contig, model)
        self.assertEqual(len(stderr), 0)
        # TODO Parse out the output for the model and check that the text is
        # in there.
        self.assertGreater(len(stdout), 0)
