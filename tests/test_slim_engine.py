"""
Tests for SLiM simulation engine.
"""

import os
import io
import tempfile
import math
from unittest import mock
import numpy as np
from numpy.testing import assert_array_equal
import collections
import re
import logging
import subprocess

import pytest
import tskit
import msprime

import stdpopsim
import stdpopsim.cli
from .test_cli import capture_output

slim_path = os.environ.get("SLIM", "slim")


def count_mut_types(ts):
    selection_coeffs = [
        stdpopsim.selection_coeff_from_mutation(ts, mut) for mut in ts.mutations()
    ]
    num_neutral = sum([s == 0 for s in selection_coeffs])
    return [num_neutral, abs(len(selection_coeffs) - num_neutral)]

class TestAPI:

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    @pytest.mark.usefixtures("tmp_path")
    def test_recap_and_rescale_on_external_slim_run(self, tmp_path):
        # _SLiMEngine.simulate() adds metadata after SLiM is run.  If we
        # produce a script file, then run with SLiM externally, this metadata
        # would be omitted. We want to check that engine.recap_and_rescale
        # still runs in this case.
        treefile = tmp_path / "foo.trees"
        scriptfile = tmp_path / "slim.script"
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig(length=2489564)
        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        samples = {"pop_0": 5}
        seed = 1024
        out, _ = capture_output(
            engine.simulate,
            demographic_model=model,
            contig=contig,
            samples=samples,
            slim_script=True,
            slim_scaling_factor=10,
            seed=seed,
        )
        # re.sub processes backslash escapes given a string,
        # so to avoid this (for Windows filenames) we use a callable
        out = re.sub(
            r'defineConstant."trees_file.+;',
            lambda u: rf'defineConstant("trees_file", "{treefile}");',
            out,
        )
        with open(scriptfile, "w") as f:
            f.write(out)
        engine._run_slim(scriptfile, seed=seed)
        ts_external = tskit.load(treefile)
        ts_external = engine.recap_and_rescale(
            ts_external,
            demographic_model=model,
            contig=contig,
            samples=samples,
            slim_scaling_factor=10,
            seed=seed,
        )
        ts_internal = engine.simulate(
            demographic_model=model,
            contig=contig,
            samples=samples,
            slim_scaling_factor=10,
            seed=seed,
        )
        tables1 = ts_external.dump_tables()
        tables2 = ts_internal.dump_tables()
        assert tables1.nodes == tables2.nodes
        assert tables1.edges == tables2.edges
        assert tables1.mutations == tables2.mutations


class TestCLI:
    def docmd(self, _cmd):
        cmd = (
            f"-e slim --slim-burn-in 0 {_cmd} --left 10000000 --right 10010000 "
            "-c chr1 -s 1234"
        ).split()
        # return capture_output(stdpopsim.cli.stdpopsim_main, cmd)
        stdpopsim.cli.stdpopsim_main(cmd)

    def test_runs(self):
        slim_path = os.environ.get("SLIM", "slim")
        print("slim path:", slim_path)
        with subprocess.Popen(
            [slim_path, "-v"],
            bufsize=1,
            universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ) as proc:
            for line in proc.stdout:
                print("stdout", line)
            stderr = proc.stderr.read()
            for line in stderr.splitlines():
                print("stderr", line)
        print("return code:", proc.returncode)
        print("-----------")
        assert proc.returncode == 0
        with subprocess.Popen(
            [slim_path, "-d", 'trees_file="temp.trees"', "temp.slim"],
            bufsize=1,
            universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ) as proc:
            outs, errs = proc.communicate(timeout=15)
            for line in outs.splitlines():
                print("stdout", line)
            for line in errs.splitlines():
                print("stderr", line)
        print("return code:", proc.returncode)
        assert proc.returncode == 0
        print("-----------")
        print("all good?")

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    @pytest.mark.filterwarnings("ignore:.*has only.*individuals alive")
    @pytest.mark.usefixtures("tmp_path")
    def test_simulate(self, tmp_path):
        saved_slim_env = os.environ.get("SLIM")
        slim_path = os.environ.get("SLIM", "slim")
        fname = tmp_path / "sim.trees"
        print("fname", fname)
        # stdout, stderr = self.docmd(
        self.docmd(
            f"--slim-scaling-factor 20 --slim-path {slim_path} -v -v -v -v HomSap "
            f"pop_0:5 -o {fname}"
        )
        # print("::: STDOUT")
        # print(stdout)
        # print("::: STDERR")
        # print(stderr)
        ts = tskit.load(fname)
        # assert ts.num_samples == 10
        # assert all(tree.num_roots == 1 for tree in ts.trees())

        # if saved_slim_env is None:
        #     del os.environ["SLIM"]
        # else:
        #     os.environ["SLIM"] = saved_slim_env

        # fname = tmp_path / "sim2.trees"
        # self.docmd(f"--slim-scaling-factor 20 HomSap pop_0:5 -o {fname}")
        # ts = tskit.load(fname)
        # assert ts.num_samples == 10

        # # verify sample counts for a multipopulation demographic model
        # fname = tmp_path / "sim3.trees"
        # cmd = (
        #     "-q -e slim --slim-scaling-factor 20 --slim-burn-in 0 "
        #     f"HomSap -o {fname} -l 0.001 -c chr1 -s 1234 "
        #     "-d OutOfAfrica_3G09 YRI:0 CEU:0 CHB:4"
        # ).split()
        # capture_output(stdpopsim.cli.stdpopsim_main, cmd)
        # ts = tskit.load(fname)
        # assert ts.num_populations == 3
        # observed_counts = [0, 0, 0]
        # for sample in ts.samples():
        #     observed_counts[ts.get_population(sample)] += 1
        # assert observed_counts == [0, 0, 8]
        # assert all(tree.num_roots == 1 for tree in ts.trees())


class DontTestThis:

    def verify_slim_sim(self, ts, num_samples):
        assert ts.num_samples == num_samples
        # regions of missing data will have 0 roots
        assert all(tree.num_roots <= 1 for tree in ts.trees(root_threshold=2))
        n_mut_types = count_mut_types(ts)
        assert n_mut_types[0] > 0
        assert n_mut_types[1] > 0

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_dfe_no_demography(self):
        with tempfile.NamedTemporaryFile(mode="w") as f:
            cmd = (
                f"-q -e slim --slim-scaling-factor 20 --slim-path {slim_path} "
                f"HomSap -c chr22 -l 0.02 -o {f.name} --dfe Gamma_K17 -s 24 "
                f"pop_0:5"
            ).split()
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)
            ts = tskit.load(f.name)
        self.verify_slim_sim(ts, num_samples=10)

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_dfe_no_interval(self):
        with tempfile.NamedTemporaryFile(mode="w") as f:
            cmd = (
                f"-q -e slim --slim-scaling-factor 20 --slim-path {slim_path} "
                f"HomSap -c chr22 -l 0.01 -o {f.name} --dfe Gamma_K17 -s 984 "
                f"pop_0:5"
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
                f"--dfe-interval 1000,100000 pop_0:5"
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
                "--dfe-interval 1000,100000 YRI:5"
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
                "--dfe-annotation ensembl_havana_104_CDS pop_0:5"
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
            f"--dfe-bed-file {tmp_path / 'ex.bed'} pop_0:5"
        ).split()
        capture_output(stdpopsim.cli.stdpopsim_main, cmd)
        ts = tskit.load(fname)
        self.verify_slim_sim(ts, num_samples=10)

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_chromosomal_segment(self):
        left = 101024
        right = 200111
        sp = "HomSap"
        chrom = "chr22"
        with tempfile.NamedTemporaryFile(mode="w") as f:
            cmd = (
                f"-q -e slim --slim-scaling-factor 20 --slim-path {slim_path} "
                f"{sp} -c {chrom} --left {left} --right {right} -o {f.name} "
                f"--dfe Gamma_K17 -s 24 pop_0:5"
            ).split()
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)
            ts = tskit.load(f.name)
        species = stdpopsim.get_species(sp)
        contig = species.get_contig(chrom)
        assert ts.sequence_length == contig.length
        self.verify_slim_sim(ts, num_samples=10)

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_chromosomal_segment_with_dfe_interval(self):
        left = 101024
        right = 200111
        dfe_left = left - 5000
        dfe_right = left + 90000
        sp = "HomSap"
        chrom = "chr21"
        with tempfile.NamedTemporaryFile(mode="w") as f:
            cmd = (
                f"-q -e slim --slim-scaling-factor 20 --slim-path {slim_path} "
                f"{sp} -c {chrom} --left {left} --right {right} -o {f.name} "
                f"--dfe Gamma_K17 -s 984 --dfe-interval {dfe_left},{dfe_right} "
                f"pop_0:5"
            ).split()
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)
            ts = tskit.load(f.name)
        species = stdpopsim.get_species(sp)
        contig = species.get_contig(chrom)
        assert ts.sequence_length == contig.length
        self.verify_slim_sim(ts, num_samples=10)

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_chromosomal_segment_with_dfe_bed_file(self, tmp_path):
        left = 101024
        right = 300111
        sp = "HomSap"
        chrom = "chr19"
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
            f"{sp} -c {chrom} -s 1234 --left {left} --right {right} -o {fname} "
            f"--dfe Gamma_K17 -s 183 --dfe-bed-file {tmp_path / 'ex.bed'} "
            f"pop_0:5"
        ).split()
        capture_output(stdpopsim.cli.stdpopsim_main, cmd)
        ts = tskit.load(fname)
        species = stdpopsim.get_species(sp)
        contig = species.get_contig(chrom)
        assert ts.sequence_length == contig.length
        self.verify_slim_sim(ts, num_samples=10)

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_chromosomal_segment_with_dfe_annotation(self):
        left = 37500000
        right = 48000000
        sp = "HomSap"
        chrom = "chr18"
        with tempfile.NamedTemporaryFile(mode="w") as f:
            cmd = (
                f"-q -e slim --slim-scaling-factor 20 --slim-path {slim_path} "
                f"{sp} -c {chrom} --left {left} --right {right} -o {f.name} "
                f"--dfe Gamma_K17 -s 913 --dfe-annotation ensembl_havana_104_CDS "
                f"pop_0:5"
            ).split()
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)
            ts = tskit.load(f.name)
        species = stdpopsim.get_species(sp)
        contig = species.get_contig(chrom)
        assert ts.sequence_length == contig.length
        self.verify_slim_sim(ts, num_samples=10)

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    @pytest.mark.usefixtures("tmp_path")
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
            f"-d OutOfAfrica_3G09 YRI:5 "
        )

        # Intervals but no DFE
        cmd = (base_cmd + "--dfe-interval 1000,100000").split()
        with pytest.raises(
            SystemExit, match="interval has been assigned " "without a DFE"
        ):
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)

        # Annotation but no DFE
        cmd = (base_cmd + "--dfe-annotation ensembl_havana_104_exons").split()
        with pytest.raises(
            SystemExit, match="A DFE annotation has been assigned without a DFE."
        ):
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)

        # bed file but no DFE
        cmd = (base_cmd + f"--dfe-bed-file {tmp_path / 'ex.bed'}").split()
        with pytest.raises(
            SystemExit, match="A DFE bed file has been assigned without a DFE."
        ):
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)

        # annotation and dfe interval
        cmd = (
            base_cmd + "--dfe Gamma_K17 "
            "--dfe-interval 999,1000 "
            "--dfe-annotation ensembl_havana_104_exons"
        ).split()
        with pytest.raises(
            SystemExit, match="A DFE annotation and a DFE interval have been"
        ):
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)

        # bed file and interval
        cmd = (
            base_cmd + "--dfe Gamma_K17 "
            "--dfe-interval 999,1000 "
            f"--dfe-bed-file {tmp_path / 'ex.bed'}"
        ).split()
        with pytest.raises(
            SystemExit, match="A DFE bed file and a DFE interval have been"
        ):
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)

        # bed file and annotation
        cmd = (
            base_cmd + "--dfe Gamma_K17 "
            f"--dfe-bed-file {tmp_path / 'ex.bed'} "
            "--dfe-annotation ensembl_havana_104_exons"
        ).split()
        with pytest.raises(SystemExit, match="A DFE bed file and a DFE annotation"):
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)

        # keep mutation ids as alleles without slim engine
        cmd = (
            f"-q -e msprime HomSap -c chr22 -s 1234 -l 0.01 -o {fname} "
            f"-d OutOfAfrica_3G09 --keep-mutation-ids-as-alleles YRI:5"
        ).split()
        with pytest.raises(SystemExit, match="only applies to the SLiM engine"):
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
                self.docmd(f"HomSap pop_0:5 --dry-run -o {f.name}")
        mocked_popen.assert_called_once()
        slim_path = os.environ.get("SLIM", "slim")
        assert slim_path in mocked_popen.call_args[0][0]
        with tempfile.NamedTemporaryFile(mode="w") as f:
            self.docmd(f"HomSap pop_0:5 --dry-run -o {f.name}")
            assert os.stat(f.name).st_size == 0

    def test_bad_slim_environ_var(self):
        saved_slim_env = os.environ.get("SLIM")

        os.environ["SLIM"] = "nonexistent"
        with pytest.raises(FileNotFoundError):
            self.docmd("HomSap pop_0:5")

        if saved_slim_env is None:
            del os.environ["SLIM"]
        else:
            os.environ["SLIM"] = saved_slim_env

    def test_bad_slim_path(self):
        saved_slim_env = os.environ.get("SLIM")

        with pytest.raises(FileNotFoundError):
            self.docmd("--slim-path nonexistent HomSap pop_0:5")

        if saved_slim_env is None:
            del os.environ["SLIM"]
        else:
            os.environ["SLIM"] = saved_slim_env


class DontTestPloidy:
    """
    Test that population sizes used in SLiM engine are scaled correctly
    with regards to ploidy.

    - Check that (diploid) population size used in SLiM simulation equals
    model population size * ploidy / 2.

    - Check that `recap_epoch` returned by `slim_makescript`, used to parameterize
    recapitation with msprime, contains population size from the original model
    (e.g. # of individuals).

    - Check that individuals in tree sequence are haploid/diploid.
    """

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    @pytest.mark.usefixtures("caplog")
    def test_slim_population_size_haploid(self, caplog):
        N = 100
        engine = stdpopsim.get_engine("slim")
        contig = stdpopsim.Contig.basic_contig(length=1000, ploidy=1)
        model = stdpopsim.PiecewiseConstantSize(N)
        with caplog.at_level(logging.DEBUG):
            engine.simulate(model, contig, samples={"pop_0": 2}, verbosity=2)
        log_str = " ".join([rec.getMessage() for rec in caplog.records])
        # match: "1: p = sim.addSubpop(0, <SLiM population size>);"
        extract_ne = re.compile(".+1: p = sim.addSubpop\\(0, ([0-9]+)\\).+")
        match = extract_ne.match(log_str)
        assert match is not None
        (Nslim,) = match.groups()
        assert int(Nslim) == int(N / 2)

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    @pytest.mark.usefixtures("caplog")
    def test_slim_population_size_diploid(self, caplog):
        N = 100
        engine = stdpopsim.get_engine("slim")
        contig = stdpopsim.Contig.basic_contig(length=1000, ploidy=2)
        model = stdpopsim.PiecewiseConstantSize(N)
        with caplog.at_level(logging.DEBUG):
            engine.simulate(model, contig, samples={"pop_0": 2}, verbosity=2, seed=9)
        log_str = " ".join([rec.getMessage() for rec in caplog.records])
        # match: "1: p = sim.addSubpop(0, <SLiM population size>);"
        extract_ne = re.compile(".+1: p = sim.addSubpop\\(0, ([0-9]+)\\).+")
        match = extract_ne.match(log_str)
        assert match is not None
        (Nslim,) = match.groups()
        assert int(Nslim) == int(N)

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_recap_population_size(self):
        N = 100
        model = stdpopsim.PiecewiseConstantSize(N)
        for ploidy in [1, 2]:
            contig = stdpopsim.Contig.basic_contig(length=1000, ploidy=ploidy)
            sample_sets = model.get_sample_sets({"pop_0": 2}, ploidy=contig.ploidy)
            rate_map = stdpopsim.get_slim_mutation_rate_map(contig)
            with open(os.devnull, "w") as scriptfile:
                recap_epoch = stdpopsim.slim_makescript(
                    scriptfile,
                    "unused",
                    model,
                    contig,
                    sample_sets,
                    extended_events=None,
                    scaling_factor=1,
                    burn_in=0,
                    slim_rate_map=rate_map,
                )
            assert int(recap_epoch.populations[0].start_size) == N

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_individual_ploidy(self):
        N = 100
        model = stdpopsim.PiecewiseConstantSize(N)
        engine = stdpopsim.get_engine("slim")
        for ploidy in [1, 2]:
            contig = stdpopsim.Contig.basic_contig(length=1000, ploidy=ploidy)
            ts = engine.simulate(model, contig, samples={"pop_0": 2}, seed=8)
            assert ts.num_individuals == 2
            assert ts.num_samples == 2 * ploidy
            individual = ts.tables.nodes.individual
            assert len(np.unique(individual[: ts.num_samples])) == ts.num_individuals

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_haploidize_individuals(self):
        N = 100
        model = stdpopsim.PiecewiseConstantSize(N)
        engine = stdpopsim.get_engine("slim")
        contig = stdpopsim.Contig.basic_contig(length=1000, ploidy=2)
        ts = engine.simulate(model, contig, samples={"pop_0": 3}, seed=7)
        ts_hap = stdpopsim.utils.haploidize_individuals(ts)
        assert ts_hap.num_individuals == ts.num_individuals * 2
        for i, j in zip(ts.samples(), ts_hap.samples()):
            assert i == j
            ind_i = ts.nodes_individual[i]
            ind_j = ts_hap.nodes_individual[j]
            assert ts.tables.individuals[ind_i] == ts_hap.tables.individuals[ind_j]

def get_test_contig(
    spp="HomSap",
    length=50818,
):
    species = stdpopsim.get_species(spp)
    contig = species.get_contig(length=length)
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
    samples = {"pop_0": 50}
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
class TestSelectiveSweep(PiecewiseConstantSizeMixin):
    @staticmethod
    def _get_island_model(Ne=1000, migration_rate=0.01):
        model = msprime.Demography()
        model.add_population(initial_size=Ne, name="pop_0")
        model.add_population(initial_size=Ne, name="pop_1")
        model.set_migration_rate(source="pop_0", dest="pop_1", rate=migration_rate)
        model.set_migration_rate(source="pop_1", dest="pop_0", rate=migration_rate)
        return stdpopsim.DemographicModel(
            id="🏝️",
            description="island model",
            long_description="for sweep tests",
            model=model,
            generation_time=1,
        )

    @staticmethod
    def _fitness_per_generation(
        logfile, start_generation_ago, end_generation_ago, pop=0
    ):
        # NB: Be careful of rounding error in log with weak selection
        assert start_generation_ago >= end_generation_ago
        # Pull out average fitness for every generation in a specified time
        # interval from the fitness logfile. Columns are:
        #   "tick", "fitness_p0_mean", "fitness_p0_sd"
        # This test will break if the logfile format changes.
        log = np.genfromtxt(logfile, delimiter=",", names=True)
        tick = log["tick"].astype(int)
        # Because of save-restore during rejection sampling, the log may
        # contain multiple trajectories; keep the one that passed (last
        # fitness value for a given generation).
        accepted_trajectory = {
            i: x for i, x in zip(tick, log["fitness_p" + str(pop) + "_mean"])
        }
        generation = np.array(sorted(accepted_trajectory))
        mean_fitness = np.array([accepted_trajectory[i] for i in generation])
        generation = np.max(generation) - generation
        # Generations in which trajectory was rejected
        rejected = [tick[i + 1] <= tick[i] for i in range(len(tick) - 1)] + [False]
        rejections_by_generation = collections.Counter(np.max(tick) - tick[rejected])
        # Check that logging interval hits sweep boundaries exactly.  Fitness
        # is calculated in the early() block; so in the logfile, mean fitness
        # refers to the fitness of the generation **previous** to the tick.
        start_generation_ago = max(start_generation_ago - 1, 0)
        end_generation_ago = max(end_generation_ago - 1, 0)
        in_sweep = np.logical_and(
            generation <= start_generation_ago,
            generation >= end_generation_ago,
        )
        assert np.max(generation[in_sweep]) == start_generation_ago
        assert np.min(generation[in_sweep]) == end_generation_ago
        assert np.sum(in_sweep) == start_generation_ago - end_generation_ago + 1
        return (
            mean_fitness[in_sweep],
            mean_fitness[~in_sweep],
            rejections_by_generation,
        )

    @pytest.mark.usefixtures("tmp_path")
    def test_sweep_with_background_selection(self, tmp_path):
        logfile = tmp_path / "sweep_bgs.logfile"
        engine = stdpopsim.get_engine("slim")
        model = self._get_island_model()
        contig = get_test_contig()
        contig.add_dfe(
            intervals=np.array([[0, contig.length // 2]], dtype="int"),
            DFE=stdpopsim.DFE(
                id="BGS",
                description="deleterious",
                long_description="mutations",
                mutation_types=[
                    stdpopsim.MutationType(
                        distribution_type="e",
                        distribution_args=[-0.01],
                    )
                ],
            ),
        )
        locus_id = "sweep"
        mutation_generation_ago = 1000
        contig.add_single_site(
            id=locus_id,
            coordinate=100,
        )
        extended_events = stdpopsim.selective_sweep(
            single_site_id=locus_id,
            population="pop_1",
            mutation_generation_ago=mutation_generation_ago,
            selection_coeff=100.0,
            globally_adaptive=False,
        )
        engine.simulate(
            demographic_model=model,
            contig=contig,
            samples=self.samples,
            extended_events=extended_events,
            slim_burn_in=1,
            logfile=logfile,
            logfile_interval=1,
            seed=642,
        )
        p0_in_sweep, p0_outside_sweep, _ = self._fitness_per_generation(
            logfile=logfile,
            start_generation_ago=mutation_generation_ago,
            end_generation_ago=0,
            pop=0,
        )
        p1_in_sweep, p1_outside_sweep, _ = self._fitness_per_generation(
            logfile=logfile,
            start_generation_ago=mutation_generation_ago,
            end_generation_ago=0,
            pop=1,
        )
        # Population 0 = BGS only; Population 1 = BGS + sweep
        assert np.all(p0_in_sweep <= 1)
        assert np.all(p0_outside_sweep <= 1)
        assert np.all(p1_in_sweep > 1)
        assert np.all(p1_outside_sweep <= 1)
