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


class TestCLI:
    def docmd(self, _cmd):
        cmd = (
            f"-q -e slim --slim-burn-in 0 {_cmd} --left 10000000 --right 10010000 "
            "-c chr1 -s 1234"
        ).split()
        return capture_output(stdpopsim.cli.stdpopsim_main, cmd)

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    @pytest.mark.filterwarnings("ignore:.*has only.*individuals alive")
    @pytest.mark.usefixtures("tmp_path")
    def test_simulate(self, tmp_path):
        saved_slim_env = os.environ.get("SLIM")
        slim_path = os.environ.get("SLIM", "slim")
        fname = tmp_path / "sim.trees"
        self.docmd(
            f"--slim-scaling-factor 20 --slim-path {slim_path} HomSap "
            f"pop_0:5 -o {fname}"
        )
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
