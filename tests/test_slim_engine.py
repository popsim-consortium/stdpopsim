"""
Tests for SLiM simulation engine.
"""
import os
import io
import sys
import tempfile
import math
from unittest import mock
import numpy as np
from numpy.testing import assert_array_equal

import pytest
import tskit
import msprime

import stdpopsim
import stdpopsim.cli
from .test_cli import capture_output

IS_WINDOWS = sys.platform.startswith("win")
slim_path = os.environ.get("SLIM", "slim")


def count_mut_types(ts):
    mut_info = {}
    for mut in ts.mutations():
        for j, md in zip(mut.derived_state.split(","), mut.metadata["mutation_list"]):
            if j not in mut_info:
                mut_info[int(j)] = md
    num_neutral = sum([mut_info[j]["selection_coeff"] == 0.0 for j in mut_info])
    return [num_neutral, abs(len(mut_info) - num_neutral)]


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class TestAPI:
    def test_bad_params(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr1")
        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        samples = model.get_samples(10)

        for scaling_factor in (0, -1, -1e-6):
            with pytest.raises(ValueError):
                engine.simulate(
                    demographic_model=model,
                    contig=contig,
                    samples=samples,
                    slim_scaling_factor=scaling_factor,
                    dry_run=True,
                )

        for burn_in in (-1, -1e-6):
            with pytest.raises(ValueError):
                engine.simulate(
                    demographic_model=model,
                    contig=contig,
                    samples=samples,
                    slim_burn_in=burn_in,
                    dry_run=True,
                )

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    @pytest.mark.filterwarnings(
        "ignore:.*model has mutation rate.*but this simulation used.*"
    )
    def test_script_generation(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr1")

        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        samples = model.get_samples(10)
        out, _ = capture_output(
            engine.simulate,
            demographic_model=model,
            contig=contig,
            samples=samples,
            slim_script=True,
        )
        assert "sim.registerLateEvent" in out

        model = species.get_demographic_model("AncientEurasia_9K19")
        samples = model.get_samples(10, 20, 30, 40, 50, 60, 70)
        out, _ = capture_output(
            engine.simulate,
            demographic_model=model,
            contig=contig,
            samples=samples,
            slim_script=True,
        )
        assert "sim.registerLateEvent" in out

        model = species.get_demographic_model("AmericanAdmixture_4B11")
        samples = model.get_samples(10, 10, 10)
        out, _ = capture_output(
            engine.simulate,
            demographic_model=model,
            contig=contig,
            samples=samples,
            slim_script=True,
        )
        assert "sim.registerLateEvent" in out

    @pytest.mark.filterwarnings("ignore:Recombination map has length:UserWarning")
    def test_recombination_map(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr1", genetic_map="HapMapII_GRCh37")
        model = stdpopsim.PiecewiseConstantSize(100)
        samples = model.get_samples(10)
        engine.simulate(
            demographic_model=model,
            contig=contig,
            samples=samples,
            slim_burn_in=0.1,
        )

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_simulate(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("AraTha")
        contig = species.get_contig("5", length_multiplier=0.001)
        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        samples = model.get_samples(10)
        ts = engine.simulate(
            demographic_model=model,
            contig=contig,
            samples=samples,
            slim_scaling_factor=10,
            slim_burn_in=0,
        )
        assert ts.num_samples == 10
        assert all(tree.num_roots == 1 for tree in ts.trees())

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_simulate_verbosity(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("AraTha")
        contig = species.get_contig("5", length_multiplier=0.001)
        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        samples = model.get_samples(10)
        for v in [0, 1, 2, 3]:
            ts = engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                slim_scaling_factor=10,
                slim_burn_in=0,
                verbosity=v,
            )
            assert ts.num_samples == 10
            assert all(tree.num_roots == 1 for tree in ts.trees())

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    @pytest.mark.filterwarnings(
        "ignore:.*model has mutation rate.*but this simulation used.*"
    )
    def test_recap_and_rescale(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        model = species.get_demographic_model("OutOfAfrica_3G09")
        samples = model.get_samples(10, 10, 10)
        for proportion, seed in zip((0, 1), (1234, 2345)):
            contig = species.get_contig("chr22", length_multiplier=0.001)
            # need selected mutations so that SLiM produces some
            contig.add_dfe(
                intervals=np.array([[0, contig.length / 2]], dtype="int"),
                DFE=stdpopsim.DFE(
                    id="test",
                    description="test",
                    long_description="test",
                    mutation_types=[
                        stdpopsim.MutationType(
                            distribution_type="n",
                            distribution_args=[0, 0.01],
                        )
                    ],
                ),
            )
            if proportion:
                extended_events = None
            else:
                extended_events = []
            ts1 = engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                extended_events=extended_events,
                slim_scaling_factor=10,
                slim_burn_in=0,
                seed=seed,
            )
            ts2_headless = engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                extended_events=extended_events,
                slim_scaling_factor=10,
                slim_burn_in=0,
                seed=seed,
                _recap_and_rescale=False,
            )
            ts2 = engine.recap_and_rescale(
                ts2_headless,
                demographic_model=model,
                contig=contig,
                samples=samples,
                extended_events=extended_events,
                slim_scaling_factor=10,
                seed=seed,
            )

            tables1 = ts1.dump_tables()
            tables2 = ts2.dump_tables()

            assert tables1.nodes == tables2.nodes
            assert tables1.edges == tables2.edges
            assert tables1.mutations == tables2.mutations

    def test_assert_min_version(self):
        engine = stdpopsim.get_engine("slim")
        with mock.patch(
            "stdpopsim.slim_engine._SLiMEngine.get_version", return_value="3.4"
        ):
            with pytest.raises(RuntimeError):
                engine._assert_min_version("3.5", engine.slim_path())
            with pytest.raises(RuntimeError):
                engine._assert_min_version("4.0", None)
        with mock.patch(
            "stdpopsim.slim_engine._SLiMEngine.get_version", return_value="4.0"
        ):
            engine._assert_min_version("3.5", engine.slim_path())
            engine._assert_min_version("3.6", None)


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class TestCLI:
    def docmd(self, _cmd):
        cmd = (
            f"-q -e slim --slim-burn-in 0 {_cmd} -l 0.001 -c chr1 -s 1234 10"
        ).split()
        return capture_output(stdpopsim.cli.stdpopsim_main, cmd)

    def test_script_generation(self):
        out, _ = self.docmd("--slim-script HomSap")
        assert "sim.registerLateEvent" in out

        # msprime.MassMigration demographic events, with proportion<1.0
        # low level migration
        out, _ = self.docmd("--slim-script HomSap -d AncientEurasia_9K19")
        assert "sim.registerLateEvent" in out
        # simultaneous mass migrations, with proportions summing to 1.0
        out, _ = self.docmd("--slim-script HomSap -d AmericanAdmixture_4B11")
        assert "sim.registerLateEvent" in out

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    @pytest.mark.filterwarnings("ignore:.*has only.*individuals alive")
    def test_simulate(self):
        saved_slim_env = os.environ.get("SLIM")
        slim_path = os.environ.get("SLIM", "slim")
        with tempfile.NamedTemporaryFile(mode="w") as f:
            self.docmd(
                f"--slim-scaling-factor 20 --slim-path {slim_path} HomSap -o {f.name}"
            )
            ts = tskit.load(f.name)
        assert ts.num_samples == 10
        assert all(tree.num_roots == 1 for tree in ts.trees())

        if saved_slim_env is None:
            del os.environ["SLIM"]
        else:
            os.environ["SLIM"] = saved_slim_env

        with tempfile.NamedTemporaryFile(mode="w") as f:
            self.docmd(f"--slim-scaling-factor 20 HomSap -o {f.name}")
            ts = tskit.load(f.name)
        assert ts.num_samples == 10

        # verify sample counts for a multipopulation demographic model
        with tempfile.NamedTemporaryFile(mode="w") as f:
            cmd = (
                "-q -e slim --slim-scaling-factor 20 --slim-burn-in 0 "
                f"HomSap -o {f.name} -l 0.001 -c chr1 -s 1234 "
                "-d OutOfAfrica_3G09 0 0 8"
            ).split()
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)
            ts = tskit.load(f.name)
        assert ts.num_populations == 3
        observed_counts = [0, 0, 0]
        for sample in ts.samples():
            observed_counts[ts.get_population(sample)] += 1
        assert observed_counts == [0, 0, 8]
        assert all(tree.num_roots == 1 for tree in ts.trees())

    def verify_slim_sim(self, ts, num_samples):
        assert ts.num_samples == num_samples
        assert all(tree.num_roots == 1 for tree in ts.trees())
        n_mut_types = count_mut_types(ts)
        assert n_mut_types[0] > 0
        assert n_mut_types[1] > 0

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_dfe_no_demography(self):
        with tempfile.NamedTemporaryFile(mode="w") as f:
            cmd = (
                f"-q -e slim --slim-scaling-factor 20 --slim-path {slim_path} "
                f"HomSap -c chr22 -l 0.02 -o {f.name} --dfe Gamma_K17 -s 24 10"
            ).split()
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)
            ts = tskit.load(f.name)
        self.verify_slim_sim(ts, num_samples=10)

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_dfe_interval(self):
        with tempfile.NamedTemporaryFile(mode="w") as f:
            cmd = (
                f"-q -e slim --slim-scaling-factor 20 --slim-path {slim_path} "
                f"HomSap -c chr22 -l 0.01 -o {f.name} --dfe Gamma_K17 -s 984 "
                f"--dfe-interval 1000,100000 10"
            ).split()
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)
            ts = tskit.load(f.name)
        self.verify_slim_sim(ts, num_samples=10)

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    @pytest.mark.filterwarnings("ignore:.*has only.*individuals alive")
    def test_dfe_demography(self):
        with tempfile.NamedTemporaryFile(mode="w") as f:
            cmd = (
                f"-q -e slim --slim-scaling-factor 20 --slim-path {slim_path} "
                f"HomSap -c chr22 -l 0.01 -o {f.name} "
                "-d OutOfAfrica_3G09 --dfe Gamma_K17 -s 148 "
                "--dfe-interval 1000,100000 10"
            ).split()
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)
            ts = tskit.load(f.name)
        self.verify_slim_sim(ts, num_samples=10)

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_dfe_annotation(self):
        with tempfile.NamedTemporaryFile(mode="w") as f:
            cmd = (
                f"-q -e slim --slim-scaling-factor 20 --slim-path {slim_path} "
                f"HomSap -c chr22 -o {f.name} --dfe Gamma_K17 -s 913 "
                "--dfe-annotation ensembl_havana_104_CDS 10"
            ).split()
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)
            ts = tskit.load(f.name)
        self.verify_slim_sim(ts, num_samples=10)

    # tmp_path is a pytest fixture
    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_dfe_bed_file(self, tmp_path):
        lines = [
            "\t".join(["chr22", "100000", "145000"]),
            "\t".join(["chr22", "150000", "302425"]),
        ]
        bedfile = open(tmp_path / "ex.bed", "w")
        for lin in lines:
            bedfile.write(lin + "\n")
        bedfile.close()
        fname = tmp_path / "sim.trees"
        cmd = (
            f"-q -e slim --slim-scaling-factor 20 --slim-path {slim_path} "
            f"HomSap -c chr22 -s 1234 -l 0.01 -o {fname} --dfe Gamma_K17 -s 183 "
            f"--dfe-bed-file {tmp_path / 'ex.bed'} 10"
        ).split()
        capture_output(stdpopsim.cli.stdpopsim_main, cmd)
        ts = tskit.load(fname)
        self.verify_slim_sim(ts, num_samples=10)

    # tmp_path is a pytest fixture
    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_errors(self, tmp_path):
        lines = [
            "\t".join(["chr22", "100000", "145000"]),
            "\t".join(["chr22", "150000", "302425"]),
        ]
        bedfile = open(tmp_path / "ex.bed", "w")
        for lin in lines:
            bedfile.write(lin + "\n")
        bedfile.close()
        fname = tmp_path / "sim.trees"
        base_cmd = (
            f"-q -e slim --slim-scaling-factor 20 --slim-path {slim_path} "
            f"HomSap -c chr22 -s 1234 -l 0.01 -o {fname} "
            "-d OutOfAfrica_3G09 "
        )

        # Intervals but no DFE
        cmd = (base_cmd + "--dfe-interval 1000,100000 10").split()
        with pytest.raises(
            SystemExit, match="interval has been assigned " "without a DFE"
        ):
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)

        # Annotation but no DFE
        cmd = (base_cmd + "--dfe-annotation ensembl_havana_104_exons 10").split()
        with pytest.raises(
            SystemExit, match="A DFE annotation has been assigned without a DFE."
        ):
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)

        # bed file but no DFE
        cmd = (base_cmd + f"--dfe-bed-file {tmp_path / 'ex.bed'} 10").split()
        with pytest.raises(
            SystemExit, match="A DFE bed file has been assigned without a DFE."
        ):
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)

        # annotation and dfe interval
        cmd = (
            base_cmd + "--dfe Gamma_K17 "
            "--dfe-interval 999,1000 "
            "--dfe-annotation ensembl_havana_104_exons 10"
        ).split()
        with pytest.raises(
            SystemExit, match="A DFE annotation and a DFE interval have been"
        ):
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)

        # bed file and interval
        cmd = (
            base_cmd + "--dfe Gamma_K17 "
            "--dfe-interval 999,1000 "
            f"--dfe-bed-file {tmp_path / 'ex.bed'} 10"
        ).split()
        with pytest.raises(
            SystemExit, match="A DFE bed file and a DFE interval have been"
        ):
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)

        # bed file and annotation
        cmd = (
            base_cmd + "--dfe Gamma_K17 "
            f"--dfe-bed-file {tmp_path / 'ex.bed'} "
            "--dfe-annotation ensembl_havana_104_exons 10"
        ).split()
        with pytest.raises(SystemExit, match="A DFE bed file and a DFE annotation"):
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)

    @mock.patch("stdpopsim.slim_engine._SLiMEngine.get_version", return_value="64.64")
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
        slim_path = os.environ.get("SLIM", "slim")
        assert slim_path in mocked_popen.call_args[0][0]
        with tempfile.NamedTemporaryFile(mode="w") as f:
            self.docmd(f"HomSap --dry-run -o {f.name}")
            assert os.stat(f.name).st_size == 0

    def test_bad_slim_environ_var(self):
        saved_slim_env = os.environ.get("SLIM")

        os.environ["SLIM"] = "nonexistent"
        with pytest.raises(FileNotFoundError):
            self.docmd("HomSap")

        if saved_slim_env is None:
            del os.environ["SLIM"]
        else:
            os.environ["SLIM"] = saved_slim_env

    def test_bad_slim_path(self):
        saved_slim_env = os.environ.get("SLIM")

        with pytest.raises(FileNotFoundError):
            self.docmd("--slim-path nonexistent HomSap")

        if saved_slim_env is None:
            del os.environ["SLIM"]
        else:
            os.environ["SLIM"] = saved_slim_env


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class TestWarningsAndErrors:
    """
    Checks that warning messages are printed when appropriate.
    """

    def triplet(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr22", length_multiplier=0.001)
        return engine, species, contig

    @pytest.mark.filterwarnings("error::stdpopsim.SLiMOddSampleWarning")
    @pytest.mark.filterwarnings(
        "ignore:.*model has mutation rate.*but this simulation used.*"
    )
    def test_no_odd_sample_warning_for_even_samples(self):
        engine, species, contig = self.triplet()
        model = species.get_demographic_model("OutOfAfrica_2T12")
        samples = model.get_samples(4, 6)
        engine.simulate(
            demographic_model=model,
            contig=contig,
            samples=samples,
            dry_run=True,
        )

    def test_odd_sample_warning(self):
        engine, species, contig = self.triplet()
        model = stdpopsim.PiecewiseConstantSize(100)
        samples = model.get_samples(5)
        with pytest.warns(stdpopsim.SLiMOddSampleWarning):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                dry_run=True,
            )

        model = species.get_demographic_model("OutOfAfrica_2T12")
        samples = model.get_samples(2, 5)
        with pytest.warns(stdpopsim.SLiMOddSampleWarning):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                dry_run=True,
            )

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_bad_population_size_addSubPop_warning(self):
        engine, species, contig = self.triplet()
        model = stdpopsim.PiecewiseConstantSize(100)
        samples = model.get_samples(2)

        with pytest.warns(
            stdpopsim.UnspecifiedSLiMWarning, match="has only.*individuals alive"
        ):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                slim_scaling_factor=10,
                dry_run=True,
            )

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_no_populations_in_generation1_error(self):
        engine, species, contig = self.triplet()
        model = stdpopsim.PiecewiseConstantSize(100)
        samples = model.get_samples(2)

        with pytest.raises(stdpopsim.SLiMException):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                slim_scaling_factor=200,
                dry_run=True,
            )

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_bad_population_size_addSubpopSplit_warning(self):
        engine, species, contig = self.triplet()
        model = stdpopsim.IsolationWithMigration(
            NA=1000, N1=100, N2=1000, T=1000, M12=0, M21=0
        )
        samples = model.get_samples(2)
        with pytest.warns(
            stdpopsim.UnspecifiedSLiMWarning, match="has only.*individuals alive"
        ):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                slim_scaling_factor=10,
                dry_run=True,
            )

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    @pytest.mark.filterwarnings("ignore:.*has only.*individuals alive")
    def test_bad_population_size_addSubpopSplit_error(self):
        engine, species, contig = self.triplet()
        model = stdpopsim.IsolationWithMigration(
            NA=1000, N1=100, N2=1000, T=1000, M12=0, M21=0
        )
        samples = model.get_samples(2)
        with pytest.raises(stdpopsim.SLiMException):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                slim_scaling_factor=200,
                dry_run=True,
            )

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_bad_population_size_setSubpopulationSize_warning(self):
        engine, species, contig = self.triplet()
        model = stdpopsim.PiecewiseConstantSize(100, (1000, 1000))
        samples = model.get_samples(2)
        with pytest.warns(
            stdpopsim.UnspecifiedSLiMWarning, match="has only.*individuals alive"
        ):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                slim_scaling_factor=10,
                dry_run=True,
            )

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    @pytest.mark.filterwarnings("ignore:.*has only.*individuals alive")
    def test_bad_population_size_setSubpopulationSize_error(self):
        engine, species, contig = self.triplet()
        model = stdpopsim.PiecewiseConstantSize(100, (1000, 1000))
        samples = model.get_samples(2)
        with pytest.raises(stdpopsim.SLiMException):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                slim_scaling_factor=200,
                dry_run=True,
            )

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_sample_size_too_big_error(self):
        engine, species, contig = self.triplet()
        model = stdpopsim.PiecewiseConstantSize(1000)
        samples = model.get_samples(300)

        with pytest.raises(stdpopsim.SLiMException):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                slim_scaling_factor=10,
                dry_run=True,
            )

    def exp_decline(self, N0=100, N1=1000, T=1000):
        """
        One population model with exponential decline in population size.
        Used for testing that growth rates are handled appropriately.
        """
        r = math.log(N0 / N1) / T
        pop0 = stdpopsim.models.Population(id="pop0", description="")
        return stdpopsim.DemographicModel(
            id="exp_decline",
            description="exp_decline",
            long_description="exp_decline",
            populations=[pop0],
            generation_time=1,
            population_configurations=[
                msprime.PopulationConfiguration(
                    initial_size=N0,
                    growth_rate=r,
                    metadata=pop0.asdict(),
                )
            ],
            demographic_events=[
                msprime.PopulationParametersChange(
                    time=T, initial_size=N1, growth_rate=0, population_id=0
                ),
            ],
        )

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_bad_population_size_exp_decline_warning(self):
        engine, species, contig = self.triplet()
        model = self.exp_decline()
        samples = model.get_samples(2)
        with pytest.warns(
            stdpopsim.UnspecifiedSLiMWarning, match="has only.*individuals alive"
        ):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                slim_scaling_factor=10,
                dry_run=True,
            )

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    @pytest.mark.filterwarnings("ignore:.*has only.*individuals alive")
    def test_bad_population_size_exp_decline_error(self):
        engine, species, contig = self.triplet()
        model = self.exp_decline()
        samples = model.get_samples(2)
        with pytest.raises(stdpopsim.SLiMException):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                slim_scaling_factor=200,
                dry_run=False,
            )

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    @pytest.mark.filterwarnings("ignore:.*has only.*individuals alive")
    def test_sample_size_too_big_exp_decline_error(self):
        engine, species, contig = self.triplet()
        model = self.exp_decline()
        samples = model.get_samples(30)

        with pytest.raises(stdpopsim.SLiMException):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                slim_scaling_factor=10,
                dry_run=True,
            )

    @pytest.mark.filterwarnings("error::stdpopsim.SLiMScalingFactorWarning")
    def test_no_warning_when_not_scaling(self):
        engine, species, contig = self.triplet()
        model = stdpopsim.PiecewiseConstantSize(10000)
        samples = model.get_samples(100)
        engine.simulate(
            demographic_model=model,
            contig=contig,
            samples=samples,
            dry_run=True,
        )

        engine.simulate(
            demographic_model=model,
            contig=contig,
            samples=samples,
            dry_run=True,
            slim_scaling_factor=1,
        )

        engine.simulate(
            demographic_model=model,
            contig=contig,
            samples=samples,
            dry_run=True,
            slim_scaling_factor=1.0,
        )

    @pytest.mark.parametrize("scaling_factor", [2, 4.3])
    def test_warning_when_scaling(self, scaling_factor):
        engine, species, contig = self.triplet()
        model = stdpopsim.PiecewiseConstantSize(10000)
        samples = model.get_samples(100)
        with pytest.warns(stdpopsim.SLiMScalingFactorWarning):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                slim_scaling_factor=scaling_factor,
                dry_run=True,
            )


class TestSlimAvailable:
    """
    Checks whether SLiM is available or not on platforms that support it.
    """

    def test_parser_has_options(self):
        parser = stdpopsim.cli.stdpopsim_cli_parser()
        with mock.patch("sys.exit", autospec=True):
            _, stderr = capture_output(parser.parse_args, ["--help"])
            # On windows we should have no "slim" options
            assert IS_WINDOWS == ("slim" not in stderr)

    def test_engine_available(self):
        all_engines = [engine.id for engine in stdpopsim.all_engines()]
        assert IS_WINDOWS == ("slim" not in all_engines)


def get_test_contig(spp="HomSap", chrom="chr22", length_multiplier=0.001):
    species = stdpopsim.get_species(spp)
    contig = species.get_contig(chrom, length_multiplier=length_multiplier)
    return contig


class PiecewiseConstantSizeMixin(object):
    """
    Mixin that sets up a simple demographic model used by multiple unit tests.
    """

    N0 = 1000  # size in the present
    N1 = 500  # ancestral size
    T = 500  # generations since size change occurred
    T_mut = 300  # introduce a mutation at this generation
    model = stdpopsim.PiecewiseConstantSize(N0, (T, N1))
    model.generation_time = 1
    samples = model.get_samples(100)
    contig = get_test_contig()
    mut_id = "mut"
    contig.add_single_site(
        id=mut_id,
        coordinate=100,
        description="🧬",
        long_description="👽",
    )

    def allele_frequency(self, ts):
        """
        Get the allele frequency of the drawn mutation.
        """
        # surely there's a simpler way!
        assert ts.num_mutations == 1
        alive = ts.samples(time=0)
        mut = next(ts.mutations())
        tree = ts.at(ts.site(mut.site).position)
        have_mut = [u for u in alive if tree.is_descendant(u, mut.node)]
        af = len(have_mut) / len(alive)
        return af


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class TestRecombinationMap(PiecewiseConstantSizeMixin):
    def verify_recombination_map(self, contig, ts):
        Q = ts.metadata["SLiM"]["user_metadata"]["Q"]
        ends = ts.metadata["SLiM"]["user_metadata"]["recombination_ends"]
        rates = ts.metadata["SLiM"]["user_metadata"]["recombination_rates"]
        rm = contig.recombination_map
        rm_rates = rm.rate.copy()
        rm_rates[np.isnan(rm_rates)] = 0.0
        rescaled_rates = (1 - (1 - 2 * rm_rates) ** Q) / 2
        assert np.allclose(ends, rm.right - 1)
        assert np.allclose(rates, rescaled_rates)
        # check we have no recombinations where they aren't allowed
        # ... and we need zero rates for this to be a good check
        assert np.min(rm_rates) == 0
        breaks = [x for x in ts.breakpoints() if x > 0 and x < ts.sequence_length]
        # this will have rm.position[i-1] <= breaks < rm.position[i]
        i = np.searchsorted(
            rm.position,
            breaks,
            side="right",
        )
        assert np.max(i) < rm.num_intervals
        assert np.min(i) > 0
        assert np.all(rm.rate[i - 1] > 0)

    @pytest.mark.parametrize("Q", [1, 12])
    @pytest.mark.filterwarnings("ignore:Recombination map has length:UserWarning")
    def test_chr1(self, Q):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr1", genetic_map="HapMapII_GRCh37")
        model = stdpopsim.PiecewiseConstantSize(100)
        samples = model.get_samples(10)
        ts = engine.simulate(
            demographic_model=model,
            contig=contig,
            samples=samples,
            slim_burn_in=0.1,
            verbosity=3,
        )
        self.verify_recombination_map(contig, ts)

    def test_off_by_one(self):
        # make an extreme example that tests whether we've got the endpoints
        # of recombination rates right
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("AraTha")
        contig = species.get_contig("5")
        midpoint = int(contig.length / 2)
        contig.recombination_map = msprime.RateMap(
            position=np.array([0.0, midpoint, midpoint + 1, contig.length]),
            rate=np.array([0.0, 0.1, 0.0]),
        )
        model = stdpopsim.PiecewiseConstantSize(100)
        samples = model.get_samples(10)
        ts = engine.simulate(
            demographic_model=model,
            contig=contig,
            samples=samples,
            slim_burn_in=0.1,
            verbosity=3,
        )
        self.verify_recombination_map(contig, ts)
        assert list(ts.breakpoints()) == [0.0, midpoint, contig.length]


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class TestGenomicElementTypes(PiecewiseConstantSizeMixin):

    mut_params = {
        "f": ([-0.1], [0]),
        "g": ([-0.1, 0.1], [1, 2]),
        "e": ([0.1],),
        "n": ([-0.1, 0.2],),
        "w": ([0.1, 0.4],),
        "l": ([-0.1, 0.2],),
    }
    example_mut_types = [stdpopsim.MutationType()] + [
        stdpopsim.MutationType(distribution_type=t, distribution_args=p)
        for t, params in mut_params.items()
        for p in params
    ]
    assert len(example_mut_types) == 9

    def get_example_dfes(self):
        # this is in a function because scoping is weird
        example_dfes = [
            stdpopsim.dfe.neutral_dfe(),
            stdpopsim.DFE(
                id="simple",
                description="just one mut type",
                long_description="🐺 🦊 🐕 🦝 🦮 🐩",
                mutation_types=[self.example_mut_types[4]],
            ),
            stdpopsim.DFE(
                id="less_simple",
                description="two mutation types",
                long_description="🦕 🦖 🐳 🐋 🐬 🦭🐠 🐡 🦈 🐙",
                proportions=[0.5, 0.5],
                mutation_types=self.example_mut_types[2:4],
            ),
            stdpopsim.DFE(
                id="everything",
                description="all of them",
                long_description="this has one of each distribution type",
                proportions=[
                    2
                    * j
                    / (len(self.example_mut_types) * (len(self.example_mut_types) - 1))
                    for j in range(len(self.example_mut_types))
                ],
                mutation_types=self.example_mut_types,
            ),
        ]
        return example_dfes

    def test_slim_simulates_dfes(self):
        contig = get_test_contig()
        contig.mutation_rate = 10 * contig.mutation_rate
        dfes = [
            stdpopsim.DFE(
                id=str(j),
                description=f"mut_{j}",
                long_description=f"mut_{j}",
                mutation_types=[mt],
            )
            for j, mt in enumerate(self.example_mut_types)
        ]
        n = len(dfes)
        breaks = np.linspace(0, contig.length, n + 1)
        for j, dfe in enumerate(dfes):
            interval = [breaks[j], breaks[j + 1]]
            contig.add_dfe(intervals=np.array([interval], dtype="int"), DFE=dfe)
        engine = stdpopsim.get_engine("slim")
        ts = engine.simulate(
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            seed=123,
            verbosity=3,
        )
        for j, mt in enumerate(self.example_mut_types):
            sites = np.where(
                np.logical_and(
                    ts.tables.sites.position >= breaks[j],
                    ts.tables.sites.position < breaks[j + 1],
                )
            )[0]
            assert len(sites) > 0
            for k in sites:
                s = ts.site(k)
                for mut in s.mutations:
                    for md in mut.metadata["mutation_list"]:
                        if not mt.is_neutral:
                            assert md["selection_coeff"] != 0

    def slim_metadata_key0(self, metadata, key):
        # Everything in SLiM is a vector, so you can't just put a singleton
        # into the JSON dictionary, you've got to have a vector of length 1;
        # so this method is in place of doing `metadata[key][0]` a lot:
        # see discussion at https://github.com/MesserLab/SLiM/issues/236
        assert key in metadata and len(metadata[key]) == 1
        return metadata[key][0]

    def verify_rates_match(self, match_rate, intervals, slim_rates, slim_ends):
        for left, right in intervals:
            for rate, end in zip(slim_rates, slim_ends):
                # SLiM's endpoints are inclusive
                if left <= end:
                    assert np.isclose(rate, match_rate)
                if right - 1 <= end:
                    break
            assert right - 1 <= end

    def verify_mutation_rates(self, contig, ts):
        # compare information about the mutation rate map as recorded
        # by SLiM itself in metadata to what we think it should be
        mut_rate_metadata = self.slim_metadata_key0(
            ts.metadata["SLiM"]["user_metadata"], "mutationRates"
        )
        slim_rates = mut_rate_metadata["rates"]
        slim_ends = mut_rate_metadata["ends"]
        Q = self.slim_metadata_key0(ts.metadata["SLiM"]["user_metadata"], "Q")

        for dfe, intervals in zip(contig.dfe_list, contig.interval_list):
            prop_nonneutral = 0.0
            for mt, p in zip(dfe.mutation_types, dfe.proportions):
                if not mt.is_neutral:
                    prop_nonneutral += p
            mut_rate = contig.mutation_rate * prop_nonneutral * Q
            self.verify_rates_match(mut_rate, intervals, slim_rates, slim_ends)

        # also verify any bits not covered by DFEs have no mutations
        empty_intervals = np.array([[0, int(contig.length)]], dtype="int")
        for intervals in contig.interval_list:
            empty_intervals = stdpopsim.utils.mask_intervals(
                empty_intervals,
                intervals,
            )
        self.verify_rates_match(0.0, empty_intervals, slim_rates, slim_ends)

    def verify_dfes_metadata(self, ts, ge_types):
        def _slim_proportions(dfe):
            is_dfe_neutral = True
            for mt in dfe["mutation_types"]:
                is_dfe_neutral &= mt["is_neutral"]
            # Recall a mutation type will have proportion 0.0 in SLiM if it is
            # neutral and not all the mutation types in the corresponding DFE are
            # neutral. SLiM does not allow a genomicelementtype to have all 0.0
            # proportions.
            return [
                prop if (not mt["is_neutral"]) or is_dfe_neutral else 0.0
                for prop, mt in zip(dfe["proportions"], dfe["mutation_types"])
            ]

        for i, dfe in enumerate(ts.metadata["stdpopsim"]["DFEs"]):
            if dfe["id"] == "recapitation":
                continue
            ge = self.slim_metadata_key0(ge_types, str(i))
            int_ends = []
            int_starts = []
            for inter in dfe["intervals"]:
                int_starts.append(inter[0])
                # in SLiM ends are inclusive
                int_ends.append(inter[1] - 1)
            ge_from_meta = {
                "intervalEnds": int_ends,
                "intervalStarts": int_starts,
                "mutationFractions": _slim_proportions(dfe),
                "mutationTypes": [
                    mt["slim_mutation_type_id"] for mt in dfe["mutation_types"]
                ],
            }
            ge = self.slim_metadata_key0(ge_types, str(i))
            assert ge == ge_from_meta

    def verify_genomic_elements(self, contig, ts):
        # compare information about genomic elements as recorded
        # by SLiM itself in metadata to what we think it should be
        mut_types = self.slim_metadata_key0(
            ts.metadata["SLiM"]["user_metadata"], "mutationTypes"
        )
        ge_types = self.slim_metadata_key0(
            ts.metadata["SLiM"]["user_metadata"], "genomicElementTypes"
        )
        self.verify_dfes_metadata(ts, ge_types)
        assert len(contig.dfe_list) == len(ge_types)
        for j, (dfe, intervals) in enumerate(
            zip(contig.dfe_list, contig.interval_list)
        ):
            assert str(j) in ge_types
            ge = self.slim_metadata_key0(ge_types, str(j))
            # checking that the neutral mutations have 0.0 proportion unless
            # all the mutations are neutral in this dfe
            for mt, dfe_prop, slim_prop in zip(
                dfe.mutation_types, dfe.proportions, ge["mutationFractions"]
            ):
                if mt.is_neutral and not dfe.is_neutral:
                    assert slim_prop == 0.0
                else:
                    assert slim_prop == dfe_prop
            # "+1" because SLiM's intervals are closed on both ends,
            # stdpopsim's are closed on the left, open on the right
            slim_intervals = np.column_stack(
                [
                    ge["intervalStarts"],
                    [x + 1 for x in ge["intervalEnds"]],
                ]
            )
            assert_array_equal(intervals, slim_intervals)
            slim_mut_types = ge["mutationTypes"]
            assert len(dfe.mutation_types) == len(slim_mut_types)
            for mt, mt_id in zip(dfe.mutation_types, slim_mut_types):
                assert str(mt_id) in mut_types
                slim_mt = self.slim_metadata_key0(mut_types, str(mt_id))
                assert mt.dominance_coeff == self.slim_metadata_key0(
                    slim_mt, "dominanceCoeff"
                )
                assert mt.distribution_type == self.slim_metadata_key0(
                    slim_mt, "distributionType"
                )
                assert len(mt.distribution_args) == len(slim_mt["distributionParams"])
                for a, b in zip(mt.distribution_args, slim_mt["distributionParams"]):
                    assert a == b

    def test_no_dfe(self):
        contig = get_test_contig()
        contig.clear_dfes()
        assert len(contig.dfe_list) == 0
        engine = stdpopsim.get_engine("slim")
        with pytest.raises(ValueError):
            _ = engine.simulate(
                demographic_model=self.model,
                contig=contig,
                samples=self.samples,
            )

    def test_default_dfe(self):
        contig = get_test_contig()
        assert len(contig.dfe_list) == 1
        engine = stdpopsim.get_engine("slim")
        ts = engine.simulate(
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            verbosity=3,  # to get metadata output
        )
        self.verify_genomic_elements(contig, ts)
        self.verify_mutation_rates(contig, ts)

    def test_multiple_dfes(self):
        contig = get_test_contig()
        L = contig.length
        example_dfes = self.get_example_dfes()
        for j, dfe in enumerate(example_dfes):
            contig.add_dfe(
                np.array([[10 * (j + 1), 100 * (10 - j)], [L / 2, L]], dtype="int"), dfe
            )
        assert len(contig.dfe_list) == 1 + len(example_dfes)
        engine = stdpopsim.get_engine("slim")
        ts = engine.simulate(
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            verbosity=3,  # to get metadata output
        )
        self.verify_genomic_elements(contig, ts)
        self.verify_mutation_rates(contig, ts)

    def test_unused_dfe(self):
        # test that even if a DFE ends up being unused then it'll still be there
        contig = get_test_contig()
        example_dfes = self.get_example_dfes()
        contig.add_dfe(np.array([[0, contig.length]], dtype="int"), example_dfes[0])
        contig.add_dfe(np.array([[0, contig.length]], dtype="int"), example_dfes[3])
        assert len(contig.dfe_list) == 3
        engine = stdpopsim.get_engine("slim")
        ts = engine.simulate(
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            verbosity=3,  # to get metadata output
        )
        self.verify_genomic_elements(contig, ts)
        self.verify_mutation_rates(contig, ts)

    def test_same_dfes(self):
        # test that we can add multiple DFEs with the same name
        # (maybe the same, maybe not)
        contig = get_test_contig()
        L = int(contig.length)
        dfe0 = stdpopsim.DFE(
            id="dfe",
            description="the first one",
            long_description="hello world",
            proportions=[1.0],
            mutation_types=[self.example_mut_types[5]],
        )
        dfe1 = stdpopsim.DFE(
            id="dfe",
            description="the second one",
            long_description="I'm different but have the same name! =( =(",
            proportions=[0.2, 0.8],
            mutation_types=self.example_mut_types[7:9],
        )
        contig.add_dfe(np.array([[0, 0.5 * L]], dtype="int"), dfe0)
        contig.add_dfe(np.array([[0.2 * L, 0.4 * L]], dtype="int"), dfe0)
        contig.add_dfe(np.array([[0.45 * L, L]], dtype="int"), dfe1)
        contig.add_dfe(np.array([[0.7 * L, 0.9 * L]], dtype="int"), dfe0)
        assert len(contig.dfe_list) == 5
        engine = stdpopsim.get_engine("slim")
        ts = engine.simulate(
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            verbosity=3,  # to get metadata output
        )
        self.verify_genomic_elements(contig, ts)
        self.verify_mutation_rates(contig, ts)

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_slim_produces_mutations(self):
        engine = stdpopsim.get_engine("slim")
        contig = get_test_contig()
        dfe = stdpopsim.DFE(
            id="test",
            description="non-neutral",
            long_description="",
            mutation_types=[
                stdpopsim.MutationType(),
                stdpopsim.MutationType(
                    distribution_type="l", distribution_args=[0.01, 0.2]
                ),
            ],
            proportions=[0.5, 0.5],
        )
        contig.add_dfe(np.array([[0, contig.length]], dtype="int"), dfe)
        ts = engine.simulate(
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            slim_scaling_factor=10,
            slim_burn_in=0.1,
            verbosity=3,  # to get metadata output
            _recap_and_rescale=False,
        )
        assert ts.num_sites > 0

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_no_neutral_mutations_are_simulated_by_slim(self):
        engine = stdpopsim.get_engine("slim")
        contig = get_test_contig()
        contig.mutation_rate = 10 * contig.mutation_rate
        ts = engine.simulate(
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            slim_scaling_factor=10,
            slim_burn_in=0.1,
            verbosity=3,
            _recap_and_rescale=False,
        )
        assert ts.num_sites == 0
        self.verify_genomic_elements(contig, ts)
        self.verify_mutation_rates(contig, ts)
        assert len(ts.metadata["stdpopsim"]["DFEs"]) == len(contig.dfe_list)

        contig.dfe_list[0].mutation_types = [
            stdpopsim.MutationType() for i in range(10)
        ]
        contig.dfe_list[0].proportions = [1 / 10 for i in range(10)]
        ts = engine.simulate(
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            slim_scaling_factor=10,
            slim_burn_in=0.1,
            verbosity=3,
            _recap_and_rescale=False,
        )
        assert ts.num_sites == 0
        self.verify_genomic_elements(contig, ts)
        self.verify_mutation_rates(contig, ts)

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_neutral_dfe_slim_proportions(self):
        engine = stdpopsim.get_engine("slim")
        contig = get_test_contig()
        Q = 10
        ts = engine.simulate(
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            slim_scaling_factor=Q,
            slim_burn_in=0.1,
            verbosity=3,
            _recap_and_rescale=False,
        )
        ge_types = self.slim_metadata_key0(
            ts.metadata["SLiM"]["user_metadata"], "genomicElementTypes"
        )
        assert np.allclose(ge_types["0"][0]["mutationFractions"], [1.0])
        assert np.allclose(
            ts.metadata["SLiM"]["user_metadata"]["mutationRates"][0]["rates"], [0.0 * Q]
        )
        contig.add_dfe(
            intervals=np.array([[0, contig.length]], dtype="int"),
            DFE=stdpopsim.DFE(
                id="neutral_sel",
                description="neutral with some selected",
                long_description="test",
                mutation_types=[
                    stdpopsim.MutationType(),
                    stdpopsim.MutationType(
                        distribution_type="e",
                        distribution_args=[-0.01],
                    ),
                ],
                proportions=[0.9, 0.1],
            ),
        )
        ts = engine.simulate(
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            slim_scaling_factor=Q,
            slim_burn_in=0.1,
            verbosity=3,
            _recap_and_rescale=False,
        )
        ge_types = self.slim_metadata_key0(
            ts.metadata["SLiM"]["user_metadata"], "genomicElementTypes"
        )
        assert np.allclose(ge_types["1"][0]["mutationFractions"], [0.0, 0.1])
        assert np.allclose(
            ts.metadata["SLiM"]["user_metadata"]["mutationRates"][0]["rates"],
            [contig.mutation_rate * 0.1 * Q],
        )


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class TestLogfile(PiecewiseConstantSizeMixin):
    # tmp_path is a pytest fixture
    def test_logfile(self, tmp_path):
        engine = stdpopsim.get_engine("slim")
        logfile = tmp_path / "slim.log"
        _ = engine.simulate(
            demographic_model=self.model,
            contig=self.contig,
            samples=self.samples,
            slim_burn_in=0.1,
            logfile=logfile,
        )
        with open(logfile, "r") as f:
            header = f.readline().strip().split(",")
            data = np.loadtxt(f, delimiter=",")
        assert header[0] == "generation"
        assert header[1][:8] == "fitness_"
        assert header[2][:8] == "fitness_"
        # neutral model, should have no fitness variation
        assert np.all(data[:, 1] == 1.0)
        assert np.all(data[:, 2] == 0.0)


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class TestDrawMutation(PiecewiseConstantSizeMixin):
    def test_draw_mutation_save(self):
        extended_events = [
            stdpopsim.ext.DrawMutation(
                time=self.T_mut,
                single_site_id=self.mut_id,
                population_id=0,
                save=True,
            ),
        ]
        engine = stdpopsim.get_engine("slim")
        engine.simulate(
            demographic_model=self.model,
            contig=self.contig,
            samples=self.samples,
            extended_events=extended_events,
            dry_run=True,
        )

    def test_draw_mutation_no_save(self):
        extended_events = [
            stdpopsim.ext.DrawMutation(
                time=self.T_mut,
                single_site_id=self.mut_id,
                population_id=0,
            ),
        ]
        engine = stdpopsim.get_engine("slim")
        engine.simulate(
            demographic_model=self.model,
            contig=self.contig,
            samples=self.samples,
            extended_events=extended_events,
            dry_run=True,
        )

    def test_invalid_mutation_id(self):
        engine = stdpopsim.get_engine("slim")
        for mut_id in ["deleterious", "sweep"]:
            extended_events = [
                stdpopsim.ext.DrawMutation(
                    time=self.T_mut,
                    single_site_id=mut_id,
                    population_id=0,
                ),
            ]
            with pytest.raises(ValueError):
                engine.simulate(
                    demographic_model=self.model,
                    contig=self.contig,
                    samples=self.samples,
                    extended_events=extended_events,
                    dry_run=True,
                )

    def test_no_mutation_types_defined(self):
        extended_events = [
            stdpopsim.ext.DrawMutation(
                time=self.T_mut,
                single_site_id=self.mut_id,
                population_id=0,
            ),
        ]
        contig = get_test_contig()
        contig.add_single_site(id=self.mut_id, coordinate=100)
        contig.dfe_list[1].mutation_types = []
        engine = stdpopsim.get_engine("slim")
        with pytest.raises(ValueError):
            engine.simulate(
                demographic_model=self.model,
                contig=contig,
                samples=self.samples,
                extended_events=extended_events,
                dry_run=True,
            )

    def test_multiple_mutation_types_defined(self):
        extended_events = [
            stdpopsim.ext.DrawMutation(
                time=self.T_mut,
                single_site_id="test",
                population_id=0,
            ),
        ]
        contig = get_test_contig()
        mt = stdpopsim.MutationType(
            distribution_type="f",
            dominance_coeff=1.0,
            distribution_args=[0.0],
            convert_to_substitution=False,
        )
        dfe = stdpopsim.DFE(
            id="test",
            mutation_types=[mt, mt],
            proportions=[0.5, 0.5],
            description="test test",
            long_description="test test test",
        )
        contig.add_dfe(
            intervals=np.array([[100, 101]], dtype="int"),
            DFE=dfe,
        )
        engine = stdpopsim.get_engine("slim")
        with pytest.raises(ValueError, match="must be unique"):
            engine.simulate(
                demographic_model=self.model,
                contig=contig,
                samples=self.samples,
                extended_events=extended_events,
                dry_run=True,
            )

    def test_fitness_distribution_not_fixed(self):
        extended_events = [
            stdpopsim.ext.DrawMutation(
                time=self.T_mut,
                single_site_id="test",
                population_id=0,
            ),
        ]
        contig = get_test_contig()
        contig.add_single_site(id="test", coordinate=100)
        mt = stdpopsim.MutationType(
            distribution_type="g",
            dominance_coeff=1.0,
            distribution_args=[1.0, 2.0],
            convert_to_substitution=False,
        )
        contig.dfe_list[1].mutation_types[0] = mt
        engine = stdpopsim.get_engine("slim")
        with pytest.raises(ValueError, match="fixed fitness coefficient"):
            engine.simulate(
                demographic_model=self.model,
                contig=contig,
                samples=self.samples,
                extended_events=extended_events,
                dry_run=True,
            )

    def test_multiple_intervals(self):
        extended_events = [
            stdpopsim.ext.DrawMutation(
                time=self.T_mut,
                single_site_id="test",
                population_id=0,
            ),
        ]
        contig = get_test_contig()
        mt = stdpopsim.MutationType(
            distribution_type="f",
            dominance_coeff=1.0,
            distribution_args=[0.0],
            convert_to_substitution=False,
        )
        dfe = stdpopsim.DFE(
            id="test",
            mutation_types=[mt],
            proportions=[1.0],
            description="test test",
            long_description="test test test",
        )
        contig.add_dfe(
            intervals=np.array([[100, 101], [103, 104]], dtype="int"),
            DFE=dfe,
        )
        engine = stdpopsim.get_engine("slim")
        with pytest.raises(ValueError, match="DFE with multiple intervals"):
            engine.simulate(
                demographic_model=self.model,
                contig=contig,
                samples=self.samples,
                extended_events=extended_events,
                dry_run=True,
            )

    def test_interval_too_large(self):
        extended_events = [
            stdpopsim.ext.DrawMutation(
                time=self.T_mut,
                single_site_id="test",
                population_id=0,
            ),
        ]
        contig = get_test_contig()
        mt = stdpopsim.MutationType(
            distribution_type="f",
            dominance_coeff=1.0,
            distribution_args=[0.0],
            convert_to_substitution=False,
        )
        dfe = stdpopsim.DFE(
            id="test",
            mutation_types=[mt],
            proportions=[1.0],
            description="test test",
            long_description="test test test",
        )
        contig.add_dfe(
            intervals=np.array([[100, 102]], dtype="int"),
            DFE=dfe,
        )
        engine = stdpopsim.get_engine("slim")
        with pytest.raises(ValueError, match="spans multiple sites"):
            engine.simulate(
                demographic_model=self.model,
                contig=contig,
                samples=self.samples,
                extended_events=extended_events,
                dry_run=True,
            )

    def test_duplicate_mutation_ids(self):
        extended_events = [
            stdpopsim.ext.DrawMutation(
                time=self.T_mut,
                single_site_id="test",
                population_id=0,
            ),
        ]
        contig = get_test_contig()
        contig.add_single_site(id="test", coordinate=100)
        contig.add_single_site(id="test", coordinate=110)
        engine = stdpopsim.get_engine("slim")
        with pytest.raises(ValueError, match="must be unique"):
            engine.simulate(
                demographic_model=self.model,
                contig=contig,
                samples=self.samples,
                extended_events=extended_events,
                dry_run=True,
            )

    def test_mutation_has_no_interval(self):
        extended_events = [
            stdpopsim.ext.DrawMutation(
                time=self.T_mut,
                single_site_id="test",
                population_id=0,
            ),
        ]
        contig = get_test_contig()
        contig.add_single_site(id="test", coordinate=100)
        contig.add_single_site(id="overlapping", coordinate=100)
        engine = stdpopsim.get_engine("slim")
        with pytest.raises(ValueError, match="has no coordinate"):
            engine.simulate(
                demographic_model=self.model,
                contig=contig,
                samples=self.samples,
                extended_events=extended_events,
                dry_run=True,
            )

    def test_bad_time(self):
        for time in (-1,):
            with pytest.raises(ValueError):
                stdpopsim.ext.DrawMutation(
                    time=time,
                    single_site_id="irrelevant",
                    population_id=0,
                )
        for time in (0, -1):
            with pytest.raises(ValueError):
                stdpopsim.ext.DrawMutation(
                    time=stdpopsim.ext.GenerationAfter(time),
                    single_site_id="irrelevant",
                    population_id=0,
                )


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class TestAlleleFrequencyConditioning(PiecewiseConstantSizeMixin):
    def test_save_point_creation(self):
        extended_events = [
            stdpopsim.ext.DrawMutation(
                time=self.T_mut,
                single_site_id=self.mut_id,
                population_id=0,
                save=True,
            ),
            stdpopsim.ext.ConditionOnAlleleFrequency(
                start_time=stdpopsim.ext.GenerationAfter(self.T_mut),
                end_time=0,
                single_site_id=self.mut_id,
                population_id=0,
                op=">",
                allele_frequency=0,
                save=True,
            ),
            stdpopsim.ext.ConditionOnAlleleFrequency(
                start_time=self.T_mut // 2,
                end_time=self.T_mut // 2,
                single_site_id=self.mut_id,
                population_id=0,
                op=">",
                allele_frequency=0,
                save=True,
            ),
        ]
        engine = stdpopsim.get_engine("slim")
        engine.simulate(
            demographic_model=self.model,
            contig=self.contig,
            samples=self.samples,
            extended_events=extended_events,
            dry_run=True,
        )

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_drawn_mutation_not_lost(self):
        engine = stdpopsim.get_engine("slim")
        ct = get_test_contig()
        ct.mutation_rate = 0.0
        ct.add_single_site(id="test", coordinate=100)
        extended_events = [
            stdpopsim.ext.DrawMutation(
                time=self.T_mut,
                single_site_id="test",
                population_id=0,
                save=True,
            ),
            stdpopsim.ext.ConditionOnAlleleFrequency(
                start_time=0,
                end_time=0,
                single_site_id="test",
                population_id=0,
                op=">",
                allele_frequency=0,
            ),
        ]
        ts = engine.simulate(
            demographic_model=self.model,
            contig=ct,
            samples=self.samples,
            extended_events=extended_events,
            slim_scaling_factor=10,
            slim_burn_in=0.1,
            _recap_and_rescale=False,
        )
        assert ts.num_mutations == 1

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_drawn_mutation_is_lost(self):
        engine = stdpopsim.get_engine("slim")
        ct = get_test_contig()
        ct.mutation_rate = 0.0
        ct.add_single_site(id="test", coordinate=100)
        extended_events = [
            stdpopsim.ext.DrawMutation(
                time=self.T_mut,
                single_site_id="test",
                population_id=0,
                save=True,
            ),
            stdpopsim.ext.ConditionOnAlleleFrequency(
                start_time=0,
                end_time=0,
                single_site_id="test",
                population_id=0,
                op="<=",
                allele_frequency=0,
            ),
        ]
        ts = engine.simulate(
            demographic_model=self.model,
            contig=ct,
            samples=self.samples,
            extended_events=extended_events,
            slim_scaling_factor=10,
            slim_burn_in=0.1,
            _recap_and_rescale=False,
        )
        assert ts.num_mutations == 0

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_drawn_mutation_meets_AF_threshold(self):
        engine = stdpopsim.get_engine("slim")
        ct = get_test_contig()
        ct.mutation_rate = 0.0
        ct.add_single_site(id="test", coordinate=100)
        for af_threshold, seed in zip((0.01, 0.1, 0.2), (1, 2, 3)):
            extended_events = [
                stdpopsim.ext.DrawMutation(
                    time=self.T_mut,
                    single_site_id="test",
                    population_id=0,
                    save=True,
                ),
                # Condition on desired AF at end of simulation.
                stdpopsim.ext.ConditionOnAlleleFrequency(
                    start_time=0,
                    end_time=0,
                    single_site_id="test",
                    population_id=0,
                    op=">=",
                    allele_frequency=af_threshold,
                ),
            ]
            ts = engine.simulate(
                demographic_model=self.model,
                contig=ct,
                samples=self.samples,
                extended_events=extended_events,
                slim_scaling_factor=10,
                slim_burn_in=0.1,
                seed=seed,
                _recap_and_rescale=False,
            )
            assert ts.num_mutations == 1
            assert self.allele_frequency(ts) >= af_threshold

    def test_bad_AF_conditioning_parameters(self):
        for op, af in [
            # bad ops
            ("=", 0.5),
            ("==", 0.5),
            ("!=", 0.5),
            ({}, 0.5),
            ("", 0.5),
            # bad allele frequencies
            ("<", -1),
            ("<=", 2.0),
            # condition is always false
            ("<", 0),
            (">", 1),
            # condition is always true
            (">=", 0),
            ("<=", 1),
        ]:
            with pytest.raises(ValueError):
                stdpopsim.ext.ConditionOnAlleleFrequency(
                    start_time=0,
                    end_time=0,
                    single_site_id="irrelevant",
                    population_id=0,
                    op=op,
                    allele_frequency=af,
                )

    def test_bad_times(self):
        for start_time, end_time in [(-1, 0), (0, -1), (1, 100)]:
            with pytest.raises(ValueError):
                stdpopsim.ext.ConditionOnAlleleFrequency(
                    start_time=start_time,
                    end_time=end_time,
                    single_site_id="irrelevant",
                    population_id=0,
                    op=">",
                    allele_frequency=0,
                )

    def test_bad_GenerationAfter_times(self):
        engine = stdpopsim.get_engine("slim")
        for start_time, end_time in [
            # Errors caught when the event is created.
            (-1, 0),
            (0, -1),
            (1, 100),
            (100, 100),
            (0, 0),
            # Errors caught when the GenerationAfter has been calculated.
            (1e-9, 0),
            (100 + 1e-9, 100),
        ]:
            with pytest.raises(ValueError):
                extended_events = [
                    stdpopsim.ext.DrawMutation(
                        time=self.T_mut,
                        single_site_id=self.mut_id,
                        population_id=0,
                        save=True,
                    ),
                    stdpopsim.ext.ConditionOnAlleleFrequency(
                        start_time=stdpopsim.ext.GenerationAfter(start_time),
                        end_time=end_time,
                        single_site_id=self.mut_id,
                        population_id=0,
                        op=">",
                        allele_frequency=0,
                    ),
                ]
                engine.simulate(
                    demographic_model=self.model,
                    contig=self.contig,
                    samples=self.samples,
                    extended_events=extended_events,
                    dry_run=True,
                )

    def test_op_id(self):
        op_types = stdpopsim.ext.ConditionOnAlleleFrequency.op_types
        for op in op_types:
            id = stdpopsim.ext.ConditionOnAlleleFrequency.op_id(op)
            assert 0 <= id < len(op_types)
        for op in ("==", "=", "!=", {}, ""):
            with pytest.raises(ValueError):
                id = stdpopsim.ext.ConditionOnAlleleFrequency.op_id(op)

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_conditioning_without_save(self):
        extended_events = [
            stdpopsim.ext.DrawMutation(
                time=self.T_mut,
                single_site_id=self.mut_id,
                population_id=0,
            ),
            stdpopsim.ext.ConditionOnAlleleFrequency(
                start_time=stdpopsim.ext.GenerationAfter(self.T_mut),
                end_time=0,
                single_site_id=self.mut_id,
                population_id=0,
                op=">=",
                allele_frequency=1,
            ),
        ]
        engine = stdpopsim.get_engine("slim")
        with pytest.raises(stdpopsim.SLiMException):
            # TODO: get this to fail using dry_run=True
            engine.simulate(
                demographic_model=self.model,
                contig=self.contig,
                samples=self.samples,
                extended_events=extended_events,
                slim_scaling_factor=10,
                slim_burn_in=0.1,
            )

    def test_no_drawn_mutation(self):
        extended_events = [
            stdpopsim.ext.ConditionOnAlleleFrequency(
                start_time=stdpopsim.ext.GenerationAfter(self.T_mut),
                end_time=0,
                single_site_id=self.mut_id,
                population_id=0,
                op=">",
                allele_frequency=0,
                save=True,
            ),
        ]
        engine = stdpopsim.get_engine("slim")
        with pytest.raises(ValueError):
            engine.simulate(
                demographic_model=self.model,
                contig=self.contig,
                samples=self.samples,
                extended_events=extended_events,
                dry_run=True,
            )


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class TestFixedSelectionCoefficient(PiecewiseConstantSizeMixin):
    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_drawn_mutation_has_correct_selection_coeff(self):
        engine = stdpopsim.get_engine("slim")
        contig = get_test_contig()
        mut_id = "sweep"
        contig.add_single_site(
            id=mut_id,
            coordinate=100,
            selection_coeff=0.1,
            dominance_coeff=0.5,
        )
        extended_events = [
            stdpopsim.ext.DrawMutation(
                time=self.T_mut,
                single_site_id=mut_id,
                population_id=0,
                save=True,
            ),
            stdpopsim.ext.ConditionOnAlleleFrequency(
                start_time=0,
                end_time=0,
                single_site_id=mut_id,
                population_id=0,
                op=">",
                allele_frequency=0,
            ),
        ]
        scaling_factor = 10
        ts = engine.simulate(
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            extended_events=extended_events,
            slim_scaling_factor=scaling_factor,
            slim_burn_in=0.1,
            _recap_and_rescale=False,
        )
        assert ts.num_mutations == 1
        mut = next(ts.mutations())
        mutation_list = mut.metadata["mutation_list"]
        assert len(mutation_list) == 1
        assert mutation_list[0]["selection_coeff"] == scaling_factor * 0.1


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class TestChangeMutationFitness(PiecewiseConstantSizeMixin):
    # Testing stdpopsim.ext.ChangeMutationFitness is challenging, because
    # the side-effects are not deterministic. But if we condition on fixation
    # of a drawn mutation, such a simulation will be very slow without strong
    # positive selection (because we effectively do rejection sampling on the
    # simulation until we get one that meets the allele frequency condition).
    # So if this test takes more than a few seconds to run, that's a good
    # indication that selection is not acting.
    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_positive_mutation_meets_AF_threshold(self):
        engine = stdpopsim.get_engine("slim")
        for af_threshold, seed in zip((0.5, 1), (1, 2)):
            extended_events = [
                stdpopsim.ext.DrawMutation(
                    time=self.T_mut,
                    single_site_id=self.mut_id,
                    population_id=0,
                    save=True,
                ),
                stdpopsim.ext.ChangeMutationFitness(
                    start_time=stdpopsim.ext.GenerationAfter(self.T_mut),
                    end_time=0,
                    single_site_id=self.mut_id,
                    population_id=0,
                    selection_coeff=0.1,
                    dominance_coeff=0.5,
                ),
                # Condition on AF > 0, to restore() immediately if the
                # allele is lost.
                stdpopsim.ext.ConditionOnAlleleFrequency(
                    start_time=0,
                    end_time=0,
                    single_site_id=self.mut_id,
                    population_id=0,
                    op=">",
                    allele_frequency=0,
                ),
                # Condition on desired AF at end of simulation.
                stdpopsim.ext.ConditionOnAlleleFrequency(
                    start_time=0,
                    end_time=0,
                    single_site_id=self.mut_id,
                    population_id=0,
                    op=">=",
                    allele_frequency=af_threshold,
                ),
            ]
            ts = engine.simulate(
                demographic_model=self.model,
                contig=self.contig,
                samples=self.samples,
                extended_events=extended_events,
                slim_scaling_factor=10,
                slim_burn_in=0.1,
                seed=seed,
                _recap_and_rescale=False,
            )
            assert ts.num_mutations == 1
            assert self.allele_frequency(ts) >= af_threshold

    def test_no_drawn_mutation(self):
        extended_events = [
            stdpopsim.ext.ChangeMutationFitness(
                start_time=stdpopsim.ext.GenerationAfter(self.T_mut),
                end_time=0,
                single_site_id=self.mut_id,
                population_id=0,
                selection_coeff=0.1,
                dominance_coeff=0.5,
            ),
        ]
        engine = stdpopsim.get_engine("slim")
        with pytest.raises(ValueError):
            engine.simulate(
                demographic_model=self.model,
                contig=self.contig,
                samples=self.samples,
                extended_events=extended_events,
                dry_run=True,
            )

    def test_bad_times(self):
        for start_time, end_time in [(-1, 0), (0, -1), (1, 100)]:
            with pytest.raises(ValueError):
                stdpopsim.ext.ChangeMutationFitness(
                    start_time=start_time,
                    end_time=end_time,
                    single_site_id="irrelevant",
                    population_id=0,
                    selection_coeff=0.1,
                    dominance_coeff=0.5,
                )

    def test_bad_GenerationAfter_times(self):
        engine = stdpopsim.get_engine("slim")
        for start_time, end_time in [
            # Errors caught when the event is created.
            (-1, 0),
            (0, -1),
            (1, 100),
            (100, 100),
            (0, 0),
            # Errors caught when the GenerationAfter has been calculated.
            (1e-9, 0),
            (100 + 1e-9, 100),
        ]:
            with pytest.raises(ValueError):
                extended_events = [
                    stdpopsim.ext.DrawMutation(
                        time=self.T_mut,
                        single_site_id=self.mut_id,
                        population_id=0,
                        save=True,
                    ),
                    stdpopsim.ext.ChangeMutationFitness(
                        start_time=stdpopsim.ext.GenerationAfter(start_time),
                        end_time=end_time,
                        single_site_id=self.mut_id,
                        population_id=0,
                        selection_coeff=0.1,
                        dominance_coeff=0.5,
                    ),
                ]
                engine.simulate(
                    demographic_model=self.model,
                    contig=self.contig,
                    samples=self.samples,
                    extended_events=extended_events,
                    dry_run=True,
                )


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class TestExtendedEvents(PiecewiseConstantSizeMixin):
    def test_bad_extended_events(self):
        engine = stdpopsim.get_engine("slim")
        for bad_ee in [
            msprime.PopulationParametersChange(time=0, initial_size=100),
            None,
            {},
            "",
        ]:
            with pytest.raises(ValueError):
                engine.simulate(
                    demographic_model=self.model,
                    contig=self.contig,
                    samples=self.samples,
                    extended_events=[bad_ee],
                    dry_run=True,
                )
