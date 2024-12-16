"""
Test cases for the command line interfaces to stdpopsim
"""

import tempfile
import pathlib
import subprocess
import json
import sys
import io
import argparse  # NOQA
import os
import logging
from unittest import mock

import tskit
import msprime
import kastore
import pytest

import stdpopsim
import stdpopsim.cli as cli


class ExceptionForTesting(Exception):
    """
    Custom exception we can throw for testing.
    """


def capture_output(func, *args, **kwargs):
    """
    Runs the specified function and arguments, and returns the
    tuple (stdout, stderr) as strings.
    """
    stdout = sys.stdout
    sys.stdout = io.StringIO()
    stderr = sys.stderr
    sys.stderr = io.StringIO()

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


class TestProvenance:
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
        assert d_stdpopsim["environment"]["os"] == d_tskit["environment"]["os"]
        assert d_stdpopsim["environment"]["python"] == d_tskit["environment"]["python"]

    def test_libraries(self):
        libs = cli.get_provenance_dict()["environment"]["libraries"]
        assert libs["tskit"]["version"] == tskit.__version__
        assert libs["msprime"]["version"] == msprime.__version__

    def test_software(self):
        software = cli.get_provenance_dict()["software"]
        assert software == {"name": "stdpopsim", "version": stdpopsim.__version__}

    def test_parameters(self):
        d = cli.get_provenance_dict()["parameters"]
        assert d["command"] == sys.argv[0]
        assert d["args"] == sys.argv[1:]


class TestDownloadGeneticMapsArgumentParser:
    """
    Tests for the download-genetic-maps parser
    """

    def test_defaults(self):
        parser = cli.stdpopsim_cli_parser()
        cmd = "download-genetic-maps"
        args = parser.parse_args([cmd])
        assert args.species is None
        assert len(args.genetic_maps) == 0

    def test_species_no_maps(self):
        parser = cli.stdpopsim_cli_parser()
        cmd = "download-genetic-maps some_species"
        args = parser.parse_args(cmd.split())
        assert args.species == "some_species"
        assert len(args.genetic_maps) == 0

    def test_species_one_map(self):
        parser = cli.stdpopsim_cli_parser()
        cmd = "download-genetic-maps some_species map1"
        args = parser.parse_args(cmd.split())
        assert args.species == "some_species"
        assert args.genetic_maps == ["map1"]

    def test_species_two_maps(self):
        parser = cli.stdpopsim_cli_parser()
        cmd = "download-genetic-maps some_species map1 map2"
        args = parser.parse_args(cmd.split())
        assert args.species == "some_species"
        assert args.genetic_maps == ["map1", "map2"]

    def test_verbosity(self):
        parser = cli.stdpopsim_cli_parser()
        cmd = "download-genetic-maps"
        args = parser.parse_args(cmd.split())
        assert args.verbose == 1

        cmd = "-v download-genetic-maps"
        args = parser.parse_args(cmd.split())
        assert args.verbose == 2

        cmd = "-vv download-genetic-maps"
        args = parser.parse_args(cmd.split())
        assert args.verbose == 3

        cmd = "-q download-genetic-maps"
        args = parser.parse_args(cmd.split())
        assert args.verbose == 0


class TestHomoSapiensArgumentParser:
    """
    Tests for the argument parsers.
    """

    def test_defaults(self):
        parser = cli.stdpopsim_cli_parser()
        cmd = "HomSap"
        args = parser.parse_args([cmd, "pop_0:2"])
        assert args.output is None
        assert args.seed is None
        assert args.samples == ["pop_0:2"]

    def test_output(self):
        parser = cli.stdpopsim_cli_parser()
        cmd = "HomSap"
        output = "/stuff/tmp.trees"

        args = parser.parse_args([cmd, "pop_0:2", "-o", output])
        assert args.output == output
        assert args.samples == ["pop_0:2"]

        args = parser.parse_args([cmd, "-o", output, "pop_0:2"])
        assert args.output == output
        assert args.samples == ["pop_0:2"]

        args = parser.parse_args([cmd, "--output", output, "pop_0:2"])
        assert args.output == output
        assert args.samples == ["pop_0:2"]

    def test_seed(self):
        parser = cli.stdpopsim_cli_parser()
        cmd = "HomSap"
        args = parser.parse_args([cmd, "pop_0:2", "-s", "1234"])
        assert args.samples == ["pop_0:2"]
        assert args.seed == 1234

        args = parser.parse_args([cmd, "pop_0:2", "--seed", "14"])
        assert args.samples == ["pop_0:2"]
        assert args.seed == 14

    def test_cache_dir(self):
        parser = cli.stdpopsim_cli_parser()
        cmd = "HomSap"
        args = parser.parse_args(["-c", "cache_dir", cmd, "pop_0:2"])
        assert args.samples == ["pop_0:2"]
        assert args.cache_dir == "cache_dir"

        args = parser.parse_args(["--cache-dir", "/some/cache_dir", cmd, "pop_0:2"])
        assert args.samples == ["pop_0:2"]
        assert args.cache_dir == "/some/cache_dir"

    def test_bibtex(self):
        parser = cli.stdpopsim_cli_parser()
        cmd = "HomSap"
        output = "/stuff/tmp.trees"
        bib = "tmp.bib"

        with mock.patch.object(argparse.FileType, "__call__", autospec=True) as call:
            args = parser.parse_args([cmd, "-b", bib, "-o", output, "pop_0:2"])
            assert args.output == output
            assert args.samples == ["pop_0:2"]
            call.assert_called_with(mock.ANY, bib)


class TestEndToEnd:
    """
    Checks that simulations we run from the CLI have plausible looking output.
    """

    def verify(self, cmd, num_samples, seed=1):
        with tempfile.TemporaryDirectory() as tmpdir:
            filename = pathlib.Path(tmpdir) / "output.trees"
            full_cmd = f" -q {cmd} -o {filename} --seed={seed}"
            with mock.patch("stdpopsim.cli.setup_logging", autospec=True):
                stdout, stderr = capture_output(cli.stdpopsim_main, full_cmd.split())
            assert len(stderr) == 0
            assert len(stdout) == 0
            ts = tskit.load(str(filename))
        assert ts.num_samples == num_samples

    def test_homsap_seed(self):
        cmd = "HomSap -c chr20 -l0.1 -s 1234 pop_0:10"
        self.verify(cmd, num_samples=20, seed=1234)

    def test_homsap_constant(self):
        cmd = "HomSap -c chr20 -l0.1 pop_0:5"
        self.verify(cmd, num_samples=10)

    def test_tennessen_two_pop_ooa(self):
        cmd = "HomSap -c chr20 -l0.1 -d OutOfAfrica_2T12 AFR:2 EUR:3"
        self.verify(cmd, num_samples=10)

    def test_gutenkunst_three_pop_ooa(self):
        cmd = "HomSap -c chr1 -l0.01 -d OutOfAfrica_3G09 YRI:5"
        self.verify(cmd, num_samples=10)

    def test_browning_america(self):
        cmd = "HomSap -c chr1 -l0.01 -d AmericanAdmixture_4B18 AFR:5"
        self.verify(cmd, num_samples=10)

    def test_ragsdale_archaic(self):
        cmd = "HomSap -c chr1 -l0.01 -d OutOfAfricaArchaicAdmixture_5R19 YRI:5"
        self.verify(cmd, num_samples=10)

    def test_dromel_constant(self):
        cmd = "DroMel -c 2L -l0.001 pop_0:2"
        self.verify(cmd, num_samples=4)

    def test_li_stephan_two_population(self):
        cmd = "DroMel -c 2L -l0.001 -d OutOfAfrica_2L06 AFR:3"
        self.verify(cmd, num_samples=6)

    def test_aratha_constant(self):
        cmd = "AraTha -L 100 pop_0:4"
        self.verify(cmd, num_samples=8)

    def test_durvusula_2017_msmc(self):
        cmd = "AraTha -L 100 -d SouthMiddleAtlas_1D17 SouthMiddleAtlas:7"
        self.verify(cmd, num_samples=14)

    def test_lapierre_constant(self):
        cmd = "EscCol -L 1000 pop_0:2"
        self.verify(cmd, num_samples=2)


class TestEndToEndSubprocess(TestEndToEnd):
    """
    Run the commands in a subprocess so that we can verify the provenance is
    stored correctly.
    """

    def verify(self, cmd, num_samples, seed=1):
        with tempfile.TemporaryDirectory() as tmpdir:
            filename = pathlib.Path(tmpdir) / "output.trees"
            full_cmd = f"{sys.executable} -m stdpopsim -q {cmd} -o {filename} -s {seed}"
            subprocess.run(full_cmd, shell=True, check=True)
            ts = tskit.load(str(filename))
        assert ts.num_samples == num_samples
        provenance = json.loads(ts.provenance(ts.num_provenances - 1).record)
        tskit.validate_provenance(provenance)
        stored_cmd = provenance["parameters"]["args"]
        assert stored_cmd[0] == "-q"
        assert stored_cmd[-1] == str(seed)
        assert stored_cmd[-2] == "-s"
        assert stored_cmd[1:-4] == cmd.split()


class TestWriteOutput:
    """
    Tests the paths through the write_output function.
    """

    def test_stdout(self):
        ts = msprime.simulate(10, random_seed=2)
        parser = cli.stdpopsim_cli_parser()
        args = parser.parse_args(["AraTha", "pop_0:2"])
        with mock.patch("shutil.copyfileobj", autospec=True) as mocked_copy:
            with mock.patch("sys.stdout", autospec=True) as stdout:
                stdout.buffer = open(os.devnull, "wb")
                cli.write_output(ts, args)
                mocked_copy.assert_called_once()

    def test_to_file(self):
        ts = msprime.simulate(10, random_seed=2)
        parser = cli.stdpopsim_cli_parser()
        output_file = "mocked.trees"
        args = parser.parse_args(["HomSap", "pop_0:2", "-o", output_file])
        with mock.patch("tskit.TreeSequence.dump", autospec=True) as mocked_dump:
            cli.write_output(ts, args)
            mocked_dump.assert_called_once_with(mock.ANY, output_file)


class TestRedirection:
    """
    Tests that the tree sequence file we get from redirecting is identical to the
    the one we get from using the --output option.
    """

    def verify_files(self, filename1, filename2):

        tables1 = tskit.load(filename1).dump_tables()
        tables2 = tskit.load(filename2).dump_tables()
        tables1.provenances.clear()
        tables2.provenances.clear()
        assert tables1 == tables2

        # Load the files into kastore to do some extra checks.
        with kastore.load(filename1) as store1, kastore.load(filename2) as store2:
            assert set(store1.keys()) == set(store2.keys())

    def verify(self, cmd):

        with tempfile.TemporaryDirectory() as tmpdir:
            filename1 = pathlib.Path(tmpdir) / "output1.trees"
            full_cmd = f"{sys.executable} -m stdpopsim {cmd} -o {filename1}"
            result = subprocess.run(
                full_cmd,
                shell=True,
                check=True,
                stderr=subprocess.PIPE,
                stdout=subprocess.PIPE,
            )
            assert len(result.stdout) == 0

            filename2 = pathlib.Path(tmpdir) / "output2.trees"
            full_cmd = f"{sys.executable} -m stdpopsim {cmd}"
            with open(filename2, "wb") as output:
                subprocess.run(
                    full_cmd,
                    shell=True,
                    check=True,
                    stdout=output,
                    stderr=subprocess.PIPE,
                )

            self.verify_files(filename1, filename2)

    def test_quiet(self):
        cmd = "-q HomSap pop_0:5 -s 2 -c chr20 -l 0.001"
        self.verify(cmd)

    def test_no_quiet(self):
        cmd = "HomSap pop_0:5 -s 3 -c chr20 -l 0.001"
        self.verify(cmd)


class TestArgumentParsing:
    """
    Tests that basic argument parsing works as expected.
    """

    basic_cmd = ["HomSap", "pop_0:10"]

    def test_quiet_verbose(self):
        parser = cli.stdpopsim_cli_parser()
        args = parser.parse_args(["-q"] + self.basic_cmd)
        assert args.verbose == 0
        args = parser.parse_args(["--quiet"] + self.basic_cmd)
        assert args.verbose == 0
        args = parser.parse_args(self.basic_cmd)
        assert args.verbose == 1
        args = parser.parse_args(["-v"] + self.basic_cmd)
        assert args.verbose == 2
        args = parser.parse_args(["--verbose"] + self.basic_cmd)
        assert args.verbose == 2
        args = parser.parse_args(["-v", "-v"] + self.basic_cmd)
        assert args.verbose == 3
        args = parser.parse_args(["-v", "--verbose"] + self.basic_cmd)
        assert args.verbose == 3
        args = parser.parse_args(["--verbose", "--verbose"] + self.basic_cmd)
        assert args.verbose == 3


class TestLogging:
    """
    Tests that logging has the desired effect.
    """

    basic_cmd = ["HomSap", "pop_0:10"]

    def test_quiet(self):
        parser = cli.stdpopsim_cli_parser()
        args = parser.parse_args(["-q"] + self.basic_cmd)
        with mock.patch("logging.basicConfig", autospec=True) as mocked_config:
            cli.setup_logging(args)
        _, kwargs = mocked_config.call_args
        assert kwargs["level"] == "ERROR"

    def test_default(self):
        parser = cli.stdpopsim_cli_parser()
        args = parser.parse_args(self.basic_cmd)
        with mock.patch("logging.basicConfig", autospec=True) as mocked_config:
            cli.setup_logging(args)
        _, kwargs = mocked_config.call_args
        assert kwargs["level"] == "WARN"

    def test_verbose(self):
        parser = cli.stdpopsim_cli_parser()
        args = parser.parse_args(["-v"] + self.basic_cmd)
        with mock.patch("logging.basicConfig", autospec=True) as mocked_config:
            cli.setup_logging(args)
        _, kwargs = mocked_config.call_args
        assert kwargs["level"] == "INFO"

    def test_very_verbose(self):
        parser = cli.stdpopsim_cli_parser()
        args = parser.parse_args(["-vv"] + self.basic_cmd)
        with mock.patch("logging.basicConfig", autospec=True) as mocked_config:
            cli.setup_logging(args)
        _, kwargs = mocked_config.call_args
        assert kwargs["level"] == "DEBUG"

    # XXX: we're now using pytest instead of unittest, so we may be able to
    # clean up this test.
    def test_logging_formatter(self):
        # unittest has its grubby mits all over logging, making it difficult
        # to test log output without using unittest.TestCase.assertLogs().
        # But assertLogs uses a hard-coded log format. So here we directly
        # construct logging.LogRecord()s for testing stdpopsim.cli.CLIFormatter.
        debug_str = "baz"
        debug_dict = dict(
            name="test_logger",
            pathname=__file__,
            levelno=logging.DEBUG,
            levelname=logging.getLevelName(logging.DEBUG),
            func=sys._getframe().f_code.co_name,
            lineno=sys._getframe().f_lineno,
            exc_info=None,
            sinfo=None,
            msg="%s",
            args=(debug_str,),
        )
        info_str = "foo"
        info_dict = dict(
            name="test_logger",
            pathname=__file__,
            levelno=logging.INFO,
            levelname=logging.getLevelName(logging.INFO),
            func=sys._getframe().f_code.co_name,
            lineno=sys._getframe().f_lineno,
            exc_info=None,
            sinfo=None,
            msg="%s",
            args=(info_str,),
        )
        warn_str = "bar"
        warn_dict = dict(
            name="py.warnings",
            pathname=__file__,
            levelno=logging.WARN,
            levelname=logging.getLevelName(logging.WARN),
            func=sys._getframe().f_code.co_name,
            lineno=sys._getframe().f_lineno,
            exc_info=None,
            sinfo=None,
            msg="%s",
            args=(
                f"{__file__}:564: UserWarning: {warn_str}\n  "
                f'warnings.warn("{warn_str}")\n',
            ),
        )
        warn_dict2 = dict(
            name="test_logger",
            pathname=__file__,
            levelno=logging.WARN,
            levelname=logging.getLevelName(logging.WARN),
            func=sys._getframe().f_code.co_name,
            lineno=sys._getframe().f_lineno,
            exc_info=None,
            sinfo=None,
            msg="%s",
            args=(warn_str,),
        )

        fmt = cli.CLIFormatter()
        formatted_debug_str = fmt.format(logging.makeLogRecord(debug_dict))
        formatted_info_str = fmt.format(logging.makeLogRecord(info_dict))
        formatted_warn_str = fmt.format(logging.makeLogRecord(warn_dict))
        formatted_warn_str2 = fmt.format(logging.makeLogRecord(warn_dict2))
        assert formatted_debug_str == "DEBUG: " + debug_str
        assert formatted_info_str == "INFO: " + info_str
        assert formatted_warn_str == "WARNING: " + warn_str
        assert formatted_warn_str2 == warn_str


class TestErrors:

    # Need to mock out setup_logging here or we spew logging to the console
    # in later tests.
    @mock.patch("stdpopsim.cli.setup_logging", autospec=True)
    def run_stdpopsim(self, command, mock_setup_logging):
        stdout, stderr = capture_output(cli.stdpopsim_main, command)
        assert stderr == ""
        assert stdout == ""
        assert mock_setup_logging.called

    def test_exit(self):
        with mock.patch(
            "sys.exit", side_effect=ExceptionForTesting, autospec=True
        ) as mocked_exit:
            with pytest.raises(ExceptionForTesting):
                cli.exit("XXX")
            mocked_exit.assert_called_once()
            args = mocked_exit.call_args[0]
            assert len(args) == 1
            assert args[0].endswith("XXX")

    # Need to mock out setup_logging here or we spew logging to the console
    # in later tests.
    @mock.patch("stdpopsim.cli.setup_logging", autospec=True)
    def verify_bad_samples(self, cmd, mock_setup_logging):
        with mock.patch(
            "stdpopsim.cli.exit", side_effect=ExceptionForTesting, autospec=True
        ) as mocked_exit:
            with pytest.raises(ExceptionForTesting):
                cli.stdpopsim_main(cmd.split())
            mocked_exit.assert_called_once()

    def test_default(self):
        self.verify_bad_samples("HomSap 2 3 ")

    def test_tennessen_model(self):
        self.verify_bad_samples("HomSap -d OutOfAfrica_2T12 2 3 4")

    def test_gutenkunst_three_pop_ooa(self):
        self.verify_bad_samples("HomSap -d OutOfAfrica_3G09 2 3 4 5")

    def test_browning_america(self):
        self.verify_bad_samples("HomSap -d AmericanAdmixture_4B18 2 3 4 5 6")

    IS_WINDOWS = sys.platform.startswith("win")

    @pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
    def test_browning_america_dfe(self):
        self.verify_bad_samples(
            "HomSap -d AmericanAdmixture_4B18 --dfe Gamma_K17 2 3 4 5 6"
        )


class TestHelp:
    def run_stdpopsim(self, command):
        with mock.patch(
            "argparse.ArgumentParser.exit",
            side_effect=ExceptionForTesting,
            autospec=True,
        ) as mocked_exit:
            with pytest.raises(ExceptionForTesting):
                capture_output(cli.stdpopsim_main, command.split())
            mocked_exit.assert_called_once()

    def test_basic_help(self):
        self.run_stdpopsim("--help")

    def test_homsap_help(self):
        self.run_stdpopsim("HomSap --help")

    def test_homsap_models_help(self):
        self.run_stdpopsim("HomSap --help-models")
        self.run_stdpopsim("HomSap --help-models OutOfAfrica_3G09")

    def test_all_species_model_help(self):
        for species in stdpopsim.all_species():
            self.run_stdpopsim(f"{species} --help-models")

    def test_homsap_genetic_maps_help(self):
        self.run_stdpopsim("HomSap --help-genetic-maps")
        self.run_stdpopsim("HomSap --help-genetic-maps HapMapII_GRCh37")

    def test_all_species_genetic_maps_help(self):
        for species in stdpopsim.all_species():
            self.run_stdpopsim(f"{species} --help-genetic-maps")

    def test_all_species_annots_help(self):
        for species in stdpopsim.all_species():
            self.run_stdpopsim(f"{species} --help-annotations")

    def test_DroMel_HomSap_dfe_help(self):
        for species in ["DroMel", "HomSap"]:
            self.run_stdpopsim(f"{species} --help-dfes")
        self.run_stdpopsim("HomSap --help-dfes Gamma_K17")

    def test_homsap_annotations_help(self):
        self.run_stdpopsim("HomSap --help-annotations")
        self.run_stdpopsim("HomSap --help-annotations ensembl_havana_104_exons")

    def test_dromel_annotations_help(self):
        self.run_stdpopsim("DroMel --help-annotations")
        self.run_stdpopsim("DroMel --help-annotations FlyBase_BDGP6.32.51_exons")


class TestWriteBibtex:
    """
    Test that citations are able to be converted to bibtex
    and written to file."""

    def test_whole_bibtex(self):
        # Test end to end
        seed = 1
        with tempfile.TemporaryDirectory() as tmpdir:
            filename = pathlib.Path(tmpdir) / "output.trees"
            bibfile = pathlib.Path(tmpdir) / "bib.bib"
            full_cmd = (
                f"HomSap -c chr20 -l0.1 YRI:10 "
                f"-o {filename} -d OutOfAfrica_3G09 --seed={seed} "
                f"--bibtex={bibfile}"
            )
            with mock.patch("stdpopsim.cli.setup_logging", autospec=True):
                with mock.patch.object(
                    stdpopsim.citations.Citation, "fetch_bibtex", autospec=True
                ) as mocked_bib:
                    with mock.patch("argparse.FileType", autospec=True):
                        stdout, stderr = capture_output(
                            cli.stdpopsim_main, full_cmd.split()
                        )
                        mocked_bib.assert_called()

    def test_number_of_calls(self):
        # Test that genetic map citations are converted.
        species = stdpopsim.get_species("HomSap")
        genetic_map = species.get_genetic_map("HapMapII_GRCh38")
        contig = species.get_contig("chr20", genetic_map=genetic_map.id)
        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        engine = stdpopsim.get_default_engine()
        dfe = species.get_dfe("Gamma_K17")
        local_cites = stdpopsim.Citation.merge(
            [
                stdpopsim.citations._stdpopsim_citation,
                stdpopsim.citations._catalog_citation,
            ]
            + genetic_map.citations
            + model.citations
            + engine.citations
            + species.genome.citations
            + species.citations
            + dfe.citations
        )
        dois = set([ref.doi for ref in local_cites])
        ncite = len(dois)
        assert ncite == len(local_cites)
        cli_cites = cli.get_citations(engine, model, contig, species, dfe)
        assert len(cli_cites) == len(local_cites)

        # Patch out writing to a file, then
        # ensure that the method is called
        # the correct number of times.
        with mock.patch("builtins.open", mock.mock_open()):
            with open("tmp.bib", "w") as bib:
                with mock.patch.object(
                    stdpopsim.citations.Citation, "fetch_bibtex", autospec=True
                ) as mock_bib:
                    cli.write_bibtex(engine, model, contig, species, bib, dfe)
                    assert mock_bib.call_count == ncite


class TestWriteCitations:
    """
    Make sure citation information is written.
    """

    @pytest.mark.usefixtures("caplog")
    def test_model_citations(self, caplog):
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("22")
        model = species.get_demographic_model("OutOfAfrica_3G09")
        engine = stdpopsim.get_default_engine()
        dfe = None
        stdout, stderr = capture_output(
            cli.write_citations, engine, model, contig, species, dfe
        )
        assert len(stdout) == 0
        genetic_map = None
        self.check_citations(engine, species, genetic_map, model, caplog.text)

        species = stdpopsim.get_species("BosTau")
        contig = species.get_contig("25")
        model = species.get_demographic_model("HolsteinFriesian_1M13")
        # this model has recombination rate, which should knock out the default
        # recombination rate and associated citation
        engine = stdpopsim.get_default_engine()
        dfe = None
        stdout, stderr = capture_output(
            cli.write_citations, engine, model, contig, species, dfe
        )
        assert len(stdout) == 0
        genetic_map = None
        self.check_citations(engine, species, genetic_map, model, caplog.text)

    @pytest.mark.usefixtures("caplog")
    def test_genetic_map_citations(self, caplog):
        species = stdpopsim.get_species("HomSap")
        genetic_map = species.get_genetic_map("HapMapII_GRCh37")
        contig = species.get_contig("chr20", genetic_map=genetic_map.id)
        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        engine = stdpopsim.get_default_engine()
        dfe = None
        stdout, stderr = capture_output(
            cli.write_citations, engine, model, contig, species, dfe
        )
        assert len(stdout) == 0
        self.check_citations(engine, species, genetic_map, model, caplog.text)

    @pytest.mark.usefixtures("caplog")
    def test_dfe_citations(self, caplog):
        species = stdpopsim.get_species("HomSap")
        genetic_map = species.get_genetic_map("HapMapII_GRCh37")
        dfe = species.get_genetic_map("HapMapII_GRCh37")
        contig = species.get_contig("chr20", genetic_map=genetic_map.id)
        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        engine = stdpopsim.get_default_engine()
        dfe = species.get_dfe("Gamma_K17")
        stdout, stderr = capture_output(
            cli.write_citations, engine, model, contig, species, dfe
        )
        assert len(stdout) == 0
        assert "[distribution of fitness effects]" in caplog.text
        assert "Kim et al., 2017" in caplog.text

    def check_citations(self, engine, species, genetic_map, model, output):
        if genetic_map is None:
            genetic_map = stdpopsim.GeneticMap(
                species,
                id="test",
                url="http://example.com/test.tgz",
                sha256="1234",
                citations=[],
            )
        for citations, assert_msg in zip(
            (engine.citations, model.citations, genetic_map.citations),
            (
                f"engine citation not written for {engine.id}",
                f"model citation not written for {model.id}",
                f"genetic map citation not written for {genetic_map.id}",
            ),
        ):
            for citation in citations:
                assert citation.author in output, assert_msg
                assert str(citation.year) in output, assert_msg
                assert citation.doi in output, assert_msg


class TestCacheDir:
    """
    Tests for setting the cache directory.
    """

    @mock.patch("stdpopsim.cli.setup_logging", autospec=True)
    @mock.patch("stdpopsim.cli.run", autospec=True)
    def run_stdpopsim(self, command, mock_setup_logging, mock_run):
        stdout, stderr = capture_output(cli.stdpopsim_main, command)
        assert stderr == ""
        assert stdout == ""
        mock_setup_logging.assert_called_once()
        mock_run.assert_called_once()

    def check_cache_dir_set(self, cmd, cache_dir):
        with mock.patch(
            "stdpopsim.set_cache_dir", autospec=True
        ) as mocked_set_cache_dir:
            self.run_stdpopsim(cmd.split())
            mocked_set_cache_dir.assert_called_once_with(cache_dir)

    def test_homsap_simulation(self):
        cache_dir = "/some/cache/dir"
        cmd = f"-c {cache_dir} HomSap pop_0:2 -o tmp.trees"
        self.check_cache_dir_set(cmd, cache_dir)

    def test_dromel_simulation(self):
        cache_dir = "cache_dir"
        cmd = f"--cache-dir {cache_dir} DroMel pop_0:2 -o tmp.trees"
        self.check_cache_dir_set(cmd, cache_dir)

    def test_download_genetic_maps(self):
        cache_dir = "/some/other/cache/dir"
        cmd = f"-c {cache_dir} download-genetic-maps"
        self.check_cache_dir_set(cmd, cache_dir)


class TestDownloadGeneticMaps:
    """
    Tests for the download genetic maps function.
    """

    def run_download(self, cmd_args, expected_num_downloads):
        parser = cli.stdpopsim_cli_parser()
        args = parser.parse_args(["download-genetic-maps"] + cmd_args.split())
        with mock.patch(
            "stdpopsim.GeneticMap.download", autospec=True
        ) as mocked_download:
            cli.run_download_genetic_maps(args)
            assert mocked_download.call_count == expected_num_downloads

    def test_defaults(self):
        num_maps = sum(len(species.genetic_maps) for species in stdpopsim.all_species())
        assert num_maps > 0
        self.run_download("", num_maps)

    def test_homsap_defaults(self):
        species = stdpopsim.get_species("HomSap")
        num_maps = len(species.genetic_maps)
        assert num_maps > 0
        self.run_download("HomSap", num_maps)

    def test_homsap_specify_maps(self):
        species = stdpopsim.get_species("HomSap")
        maps = [gmap.id for gmap in species.genetic_maps]
        for j in range(len(maps)):
            args = " ".join(maps[: j + 1])
            self.run_download("HomSap " + args, j + 1)


class TestSearchWrappers:
    """
    Tests that the search wrappers for species etc work correctly.
    """

    def test_bad_species(self):
        with mock.patch("stdpopsim.cli.exit", autospec=True) as mocked_exit:
            cli.get_species_wrapper("XXX")
            available_species = ", ".join(stdpopsim.species.registered_species)
            mocked_exit.assert_called_once_with(
                f"Species 'XXX' not in catalog ({available_species})"
            )

    def test_bad_model(self):
        species = stdpopsim.get_species("HomSap")
        with mock.patch("stdpopsim.cli.exit", autospec=True) as mocked_exit:
            cli.get_model_wrapper(species, "XXX")
            available_models = ", ".join([dm.id for dm in species.demographic_models])
            mocked_exit.assert_called_once_with(
                f"DemographicModel 'HomSap/XXX' not in catalog ({available_models})"
            )

    def test_bad_genetic_map(self):
        species = stdpopsim.get_species("HomSap")
        with mock.patch("stdpopsim.cli.exit", autospec=True) as mocked_exit:
            cli.get_genetic_map_wrapper(species, "XXX")
            available_maps = ", ".join([gm.id for gm in species.genetic_maps])
            mocked_exit.assert_called_once_with(
                f"GeneticMap 'HomSap/XXX' not in catalog ({available_maps})"
            )

    def test_bad_dfe(self):
        species = stdpopsim.get_species("HomSap")
        with mock.patch("stdpopsim.cli.exit", autospec=True) as mocked_exit:
            cli.get_dfe_wrapper(species, "XXX")
            available_dfes = [dm.id for dm in species.dfes]
            avail_dfes_str = ", ".join(available_dfes)
            mocked_exit.assert_called_once_with(
                f"DFE 'HomSap/XXX' not in catalog ({avail_dfes_str})"
            )

    def test_bad_annotation(self):
        species = stdpopsim.get_species("HomSap")
        with mock.patch("stdpopsim.cli.exit", autospec=True) as mocked_exit:
            cli.get_annotation_wrapper(species, "666_foo_666")
            available_annots = [dm.id for dm in species.annotations]
            avail_annots_str = ", ".join(available_annots)
            mocked_exit.assert_called_once_with(
                f"Annotations 'HomSap/666_foo_666' not in catalog ({avail_annots_str})"
            )


class TestDryRun:
    """
    Checks that simulations we run from the CLI with the --dry-run option have no output
    """

    def test_dry_run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = pathlib.Path(tmpdir)
            with open(path / "stderr", "w+") as stderr:
                filename = path / "output.trees"
                cmd = (
                    f"{sys.executable} -m stdpopsim HomSap -D -L 1000 -o "
                    f"{filename} pop_0:1"
                )
                subprocess.run(cmd, stderr=stderr, shell=True, check=True)
                assert stderr.tell() > 0
            assert not os.path.isfile(filename)

    def test_dry_run_quiet(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = pathlib.Path(tmpdir)
            with open(path / "stderr", "w+") as stderr:
                filename = path / "output.trees"
                cmd = (
                    f"{sys.executable} -m stdpopsim -q HomSap -D -L 1000 "
                    f"-o {filename} pop_0:1"
                )
                subprocess.run(cmd, stderr=stderr, shell=True, check=True)
                assert stderr.tell() == 0
            assert not os.path.isfile(filename)


class TestMsprimeEngine:
    def docmd(self, _cmd):
        with tempfile.TemporaryDirectory() as tmpdir:
            filename = pathlib.Path(tmpdir) / "output.trees"
            cmd = f"-q -e msprime {_cmd} AraTha -L 100 --seed 1 -o {filename} pop_0:5"
            return capture_output(stdpopsim.cli.stdpopsim_main, cmd.split())

    def test_simulate(self):
        self.docmd("")
        self.docmd("--msprime-model hudson")
        self.docmd("--msprime-model smc")
        self.docmd("--msprime-model smc_prime")
        self.docmd("--msprime-model dtwf --msprime-change-model 49.6 hudson")

        self.docmd(
            "--msprime-model hudson "
            "--msprime-change-model 10 dtwf "
            "--msprime-change-model 20 hudson "
            "--msprime-change-model 30 dtwf "
            "--msprime-change-model 40 hudson"
        )

    def test_invalid_CLI_parameters(self):
        with pytest.raises(SystemExit):
            self.docmd("--msprime-model notamodel")
        with pytest.raises(SystemExit):
            self.docmd("--msprime-model dtwf --msprime-change-model 50 notamodel")
        with pytest.raises(SystemExit):
            self.docmd("--msprime-model dtwf --msprime-change-model notanumber hudson")
        with pytest.raises(SystemExit):
            self.docmd("--msprime-model hudson --msprime-change-model dtwf")
        with pytest.raises(SystemExit):
            self.docmd("--msprime-model hudson --msprime-change-model 10")

    def test_invalid_API_parameters(self):
        engine = stdpopsim.get_engine("msprime")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr20")
        model = species.get_demographic_model("OutOfAfrica_2T12")
        samples = {"AFR": 5}
        with pytest.raises(ValueError):
            engine.simulate(model, contig, samples, msprime_model="notamodel")
        with pytest.raises(ValueError):
            engine.simulate(
                model,
                contig,
                samples,
                msprime_change_model=[
                    (10, "notamodel"),
                ],
            )


class TestNonAutosomal:
    # TODO: This test should be removed when #383 is fixed.
    # https://github.com/popsim-consortium/stdpopsim/issues/383
    def test_chrX_gives_a_warning(self):
        cmd = "HomSap -D -c chrX -o /dev/null pop_0:10".split()
        # setup_logging() interferes with pytest.warns().
        with mock.patch("stdpopsim.cli.setup_logging", autospec=True):
            with pytest.warns(stdpopsim.NonAutosomalWarning):
                capture_output(stdpopsim.cli.stdpopsim_main, cmd)


class TestRecombinationRate:
    @pytest.mark.usefixtures("caplog")
    def test_recomb_rate(self, caplog):
        sp = stdpopsim.get_species("BosTau")
        mod = sp.get_demographic_model("HolsteinFriesian_1M13")
        cmd = (
            "BosTau -D -c 1 -o /dev/null -d HolsteinFriesian_1M13 Holstein_Friesian:10"
        )
        with mock.patch("stdpopsim.cli.setup_logging", autospec=True):
            out, err = capture_output(stdpopsim.cli.stdpopsim_main, cmd.split())
            assert len(out) == 0
            assert len(err) == 0
            log_output = caplog.text
            for ll in log_output.split("\n"):
                if "recombination" in ll:
                    assert (
                        ll.split()
                        == f"Mean recombination rate: {mod.recombination_rate}".split()
                    )


class TestNoQCWarning:
    species = stdpopsim.get_species("EscCol")
    model = stdpopsim.DemographicModel(
        id="FakeModel",
        description="FakeModel",
        long_description="FakeModel",
        citations=[
            stdpopsim.Citation(
                author="Farnsworth",
                year=3000,
                doi="https://doi.org/10.1000/xyz123",
                reasons={stdpopsim.CiteReason.DEM_MODEL},
            )
        ],
        generation_time=10,
        populations=[stdpopsim.Population("Omicronians", "Popplers, etc.")],
        population_configurations=[msprime.PopulationConfiguration(initial_size=1000)],
    )

    @classmethod
    def setup_class(cls):
        cls.species.add_demographic_model(cls.model)

    @classmethod
    def teardown_class(cls):
        cls.species.demographic_models.remove(cls.model)

    def verify_noQC_warning(self, cmd):
        # setup_logging() interferes with pytest.warns().
        with mock.patch("stdpopsim.cli.setup_logging", autospec=True):
            with pytest.warns(stdpopsim.QCMissingWarning):
                capture_output(stdpopsim.cli.stdpopsim_main, cmd.split())

    def test_noQC_warning(self):
        self.verify_noQC_warning("EscCol -d FakeModel -D -L 1000 Omicronians:10")

    def test_noQC_warning_quiet(self):
        self.verify_noQC_warning("-q EscCol -d FakeModel -D -L 1000 Omicronians:10")

    def verify_noQC_citations_not_written(self, cmd, caplog):
        # Non-QCed models shouldn't be used in publications, so citations
        # shouldn't be offered to the user.
        out, err = capture_output(stdpopsim.cli.stdpopsim_main, cmd.split())
        log_output = caplog.text
        # citations are written out via the logging mechanism
        assert len(out) == 0
        assert len(err) == 0
        for citation in self.model.citations:
            assert not (citation.author in log_output)
            assert not (citation.doi in log_output)
            assert not (citation.author in log_output)
            assert not (citation.doi in log_output)
            assert not (citation.author in log_output)
            assert not (citation.doi in log_output)

    # The following two tests use the "caplog" pytest fixture, which captures
    # the logging output. The caplog param is automatically passed to test_*()
    # methods by pytest, which we pass through to verify_noQC_citations_not_written().

    @pytest.mark.filterwarnings("ignore::stdpopsim.QCMissingWarning")
    @pytest.mark.usefixtures("caplog")
    def test_noQC_citations_not_written(self, caplog):
        self.verify_noQC_citations_not_written(
            "EscCol -d FakeModel -D -L 1000 Omicronians:10", caplog
        )

    @pytest.mark.filterwarnings("ignore::stdpopsim.QCMissingWarning")
    @pytest.mark.usefixtures("caplog")
    def test_noQC_citations_not_written_verbose(self, caplog):
        self.verify_noQC_citations_not_written(
            "-vv EscCol -d FakeModel -D -L 1000 Omicronians:10", caplog
        )


@pytest.mark.parametrize("species_id", [s.id for s in stdpopsim.all_species()])
def test_species_simulation(species_id):
    # L must be longer than any mean gene conversion tract length;
    # since gene conversion length is quite long for a few species
    # we can't run with L=1000 for all species, so (confusingly) adjust
    # for those long-GC species here
    L = 1000
    species = stdpopsim.get_species(species_id)
    for chrom in species.genome.chromosomes:
        if chrom.gene_conversion_length is not None:
            L = max(L, chrom.gene_conversion_length + 100)
    cmd = f"-q {species_id} -L {L} --seed 1234 pop_0:10"
    # Just check to see if the simulation runs
    with mock.patch("sys.stdout", autospec=True) as stdout:
        stdout.buffer = open(os.devnull, "wb")
        stdpopsim.cli.stdpopsim_main(cmd.split())


class TestSampleCountParser:
    def verify(self, cmd, num_samples, seed=1):
        with tempfile.TemporaryDirectory() as tmpdir:
            filename = pathlib.Path(tmpdir) / "output.trees"
            full_cmd = f" -q {cmd} -o {filename} --seed={seed}"
            with mock.patch("stdpopsim.cli.setup_logging", autospec=True):
                stdout, stderr = capture_output(cli.stdpopsim_main, full_cmd.split())
            assert len(stderr) == 0
            assert len(stdout) == 0
            ts = tskit.load(str(filename))
        for i in range(len(num_samples)):
            assert len(ts.samples(population=i)) == num_samples[i]

    @pytest.mark.filterwarnings("ignore::stdpopsim.DeprecatedFeatureWarning")
    def test_deprecated_positional_samples(self):
        cmd = "HomSap -c chr1 -l0.01 -d OutOfAfrica_3G09 5"
        self.verify(cmd, num_samples=[5, 0, 0])
        cmd = "HomSap -c chr1 -l0.01 -d OutOfAfrica_3G09 5 10"
        self.verify(cmd, num_samples=[5, 10, 0])
        cmd = "HomSap -c chr1 -l0.01 -d OutOfAfrica_3G09 5 0 10"
        self.verify(cmd, num_samples=[5, 0, 10])

    def test_population_sample_pairs(self):
        cmd = "HomSap -c chr1 -l0.01 -d OutOfAfrica_3G09 YRI:5"
        self.verify(cmd, num_samples=[10, 0, 0])
        cmd = "HomSap -c chr1 -l0.01 -d OutOfAfrica_3G09 CHB:5 YRI:10"
        self.verify(cmd, num_samples=[20, 0, 10])
        cmd = "HomSap -c chr1 -l0.01 -d OutOfAfrica_3G09 CHB:5 YRI:10 CEU:0"
        self.verify(cmd, num_samples=[20, 0, 10])

    def test_bad_sample_specification(self):
        for cmd in [
            "HomSap -c chr1 -l0.01 -d OutOfAfrica_3G09 5 bad",
            "HomSap -c chr1 -l0.01 -d OutOfAfrica_3G09 5 CEU:5",
            "HomSap -c chr1 -l0.01 -d OutOfAfrica_3G09 YRI:5 very:bad",
            "HomSap -c chr1 -l0.01 -d OutOfAfrica_3G09 YRI=5",
        ]:
            with pytest.raises(
                ValueError, match="Sample specification must be in the form"
            ):
                self.verify(cmd, num_samples=[5, 0, 0])

    def test_duplicate_sample_specification(self):
        cmd = "HomSap -c chr1 -l0.01 -d OutOfAfrica_3G09 YRI:5 YRI:5"
        with pytest.raises(ValueError, match="specified more than once"):
            self.verify(cmd, num_samples=[20, 0, 0])
