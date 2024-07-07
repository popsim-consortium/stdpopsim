"""
Tests for SLiM simulation engine.
"""
import os
import tempfile
import numpy as np
import re
import contextlib

import pytest
import tskit

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
        contig = species.get_contig("chr19", length_multiplier=0.01)
        model = stdpopsim.PiecewiseConstantSize(species.population_size / 10)
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
        a, b = re.split(r'defineConstant\("trees_file.+;', out)
        out = a + 'defineConstant("trees_file", "' + str(treefile) + '");' + b
        # out = re.sub(
        #     r'defineConstant\("trees_file.+;',
        #     rf'defineConstant("trees_file", "{treefile}");',
        #     out,
        # )
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
        cmd = (f"-q -e slim --slim-burn-in 0 {_cmd} -l 0.001 -c chr1 -s 1234").split()
        return capture_output(stdpopsim.cli.stdpopsim_main, cmd)

    def verify_slim_sim(self, ts, num_samples):
        assert ts.num_samples == num_samples
        # regions of missing data will have 0 roots
        assert all(tree.num_roots <= 1 for tree in ts.trees(root_threshold=2))
        n_mut_types = count_mut_types(ts)
        assert n_mut_types[0] > 0
        assert n_mut_types[1] > 0

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    @pytest.mark.usefixtures("tmp_path")
    def test_dfe_no_demography(self, tmp_path):
        fname = tmp_path / "test_dfe_no_demography.trees"
        cmd = (
            f"-q -e slim --slim-scaling-factor 50 --slim-path {slim_path} "
            "--slim-script "
            f"HomSap -c chr22 -l 0.02 -o {fname} --dfe Gamma_K17 -s 24 "
            f"pop_0:5"
        ).split()
        capture_output(stdpopsim.cli.stdpopsim_main, cmd)
        cmd = (
            f"-q -e slim --slim-scaling-factor 50 --slim-path {slim_path} "
            f"HomSap -c chr22 -l 0.02 -o {fname} --dfe Gamma_K17 -s 24 "
            f"pop_0:5"
        ).split()
        capture_output(stdpopsim.cli.stdpopsim_main, cmd)
        ts = tskit.load(fname)
        self.verify_slim_sim(ts, num_samples=10)

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_dfe_no_interval(self):
        tempdir = tempfile.TemporaryDirectory(prefix="TestCLI_")

        @contextlib.contextmanager
        def temp_file():
            fname = os.path.join(tempdir.name, f"{os.urandom(3).hex()}.slim")
            f = open(fname, "w")
            yield f, fname
            f.close()

        with temp_file() as fn:
            f, fname = fn
            cmd = (
                f"-q -e slim --slim-scaling-factor 20 --slim-path {slim_path} "
                "--slim-script "
                f"HomSap -c chr22 -l 0.01 -o {fname} --dfe Gamma_K17 -s 984 "
                f"pop_0:5"
            ).split()
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)
            cmd = (
                f"-q -e slim --slim-scaling-factor 20 --slim-path {slim_path} "
                f"HomSap -c chr22 -l 0.01 -o {fname} --dfe Gamma_K17 -s 984 "
                f"pop_0:5"
            ).split()
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)
            ts = tskit.load(fname)
        self.verify_slim_sim(ts, num_samples=10)


def get_test_contig(
    spp="HomSap", chrom="chr22", length_multiplier=0.001, left=None, right=None
):
    species = stdpopsim.get_species(spp)
    contig = species.get_contig(
        chrom, length_multiplier=length_multiplier, left=left, right=right
    )
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
        assert header[0] == "tick"
        assert header[1][:8] == "fitness_"
        assert header[2][:8] == "fitness_"
        # neutral model, should have no fitness variation
        assert np.all(data[:, 1] == 1.0)
        assert np.all(data[:, 2] == 0.0)
