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
import collections
import re
import logging

import pytest
import tskit
import msprime

import stdpopsim
import stdpopsim.cli
from .test_cli import capture_output

IS_WINDOWS = sys.platform.startswith("win")
slim_path = os.environ.get("SLIM", "slim")


def count_mut_types(ts):
    selection_coeffs = [
        stdpopsim.selection_coeff_from_mutation(ts, mut) for mut in ts.mutations()
    ]
    num_neutral = sum([s == 0 for s in selection_coeffs])
    return [num_neutral, abs(len(selection_coeffs) - num_neutral)]


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class TestAPI:
    def test_bad_params(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr1")
        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        samples = {"pop_0": 5}

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

    def test_bad_samples(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr1")
        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        samples = [1, 2, ["foo"]]
        with pytest.raises(ValueError, match="Samples must be a dict"):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                dry_run=True,
            )
        with pytest.raises(ValueError, match="Samples must be a dict"):
            engine.recap_and_rescale(
                ts=None,
                demographic_model=model,
                contig=contig,
                samples=samples,
            )
        samples = [
            msprime.SampleSet(
                num_samples=2,
                population=0,
                ploidy=3,
            )
        ]
        with pytest.raises(ValueError, match="Sample ploidy other than 1 or 2"):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
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
        samples = {"pop_0": 5}
        out, _ = capture_output(
            engine.simulate,
            demographic_model=model,
            contig=contig,
            samples=samples,
            slim_script=True,
        )
        assert "community.registerLateEvent" in out

        model = species.get_demographic_model("AncientEurasia_9K19")
        samples = {
            "Mbuti": 5,
            "LBK": 10,
            "Sardinian": 15,
            "Loschbour": 20,
            "MA1": 25,
            "Han": 30,
            "UstIshim": 35,
        }
        out, _ = capture_output(
            engine.simulate,
            demographic_model=model,
            contig=contig,
            samples=samples,
            slim_script=True,
        )
        assert "community.registerLateEvent" in out

        model = species.get_demographic_model("AmericanAdmixture_4B11")
        samples = {"AFR": 10, "EUR": 10, "ASIA": 10}
        out, _ = capture_output(
            engine.simulate,
            demographic_model=model,
            contig=contig,
            samples=samples,
            slim_script=True,
        )
        assert "community.registerLateEvent" in out

    @pytest.mark.filterwarnings("ignore: Genetic map.*is longer than chromosome length")
    def test_recombination_map(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr1", genetic_map="HapMapII_GRCh37")
        model = stdpopsim.PiecewiseConstantSize(100)
        samples = {"pop_0": 5}
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
        samples = {"pop_0": 5}
        ts = engine.simulate(
            demographic_model=model,
            contig=contig,
            samples=samples,
            slim_scaling_factor=10,
            slim_burn_in=0,
        )
        assert ts.num_samples == 10
        assert all(tree.num_roots <= 1 for tree in ts.trees(root_threshold=2))

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_simulate_verbosity(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("AraTha")
        contig = species.get_contig("5", length_multiplier=0.001)
        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        samples = {"pop_0": 5}
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
            assert all(tree.num_roots <= 1 for tree in ts.trees(root_threshold=2))

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    @pytest.mark.filterwarnings("ignore::stdpopsim.DeprecatedFeatureWarning")
    @pytest.mark.filterwarnings(
        "ignore:.*model has mutation rate.*but this simulation used.*"
    )
    def test_recap_and_rescale(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        model = species.get_demographic_model("OutOfAfrica_3G09")
        samples_deprecated = model.get_samples(10, 10, 10)
        samples_dict = {"YRI": 5, "CEU": 5, "CHB": 5}
        for samples in [samples_deprecated, samples_dict]:
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
        contig = species.get_contig("chr1", length_multiplier=0.01)
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
        out = re.sub(
            'defineConstant\\("trees_file.+;',
            f'defineConstant("trees_file", "{treefile}");',
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

    def test_stacked_mutations(self):
        # Verify that `count_mut_types` works with stacked mutations, after
        # these have been converted to nucleotides. When it's possible to set
        # the stacking policy to "l" in SLiMMutationModel, this test will break
        # and can be removed.
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig(length=10)
        model = stdpopsim.PiecewiseConstantSize(100)
        samples = {"pop_0": 5}
        contig.mutation_rate = 1e-2
        while True:
            ts = engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                slim_burn_in=0,
                keep_mutation_ids_as_alleles=False,
            )
            is_stacked = [len(m.metadata["mutation_list"]) > 1 for m in ts.mutations()]
            if any(is_stacked):
                break
        count_mut_types(ts)

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_allele_codes(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("AraTha")
        contig = species.get_contig("5", left=0, right=100000)
        contig.add_dfe(
            intervals=[[0, contig.length // 2]],
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
        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        samples = {"pop_0": 5}
        # nucleotides
        ts = engine.simulate(
            demographic_model=model,
            contig=contig,
            samples=samples,
            slim_scaling_factor=10,
            slim_burn_in=0,
        )
        for mut in ts.mutations():
            assert mut.derived_state in "ACGT"
        # slim mutation ids (may be stacked)
        ts = engine.simulate(
            demographic_model=model,
            contig=contig,
            samples=samples,
            keep_mutation_ids_as_alleles=True,
            slim_scaling_factor=10,
            slim_burn_in=0,
        )
        for mut in ts.mutations():
            alleles = mut.derived_state.split(",")
            assert all([x.isnumeric() for x in alleles])


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class TestCLI:
    def docmd(self, _cmd):
        cmd = (f"-q -e slim --slim-burn-in 0 {_cmd} -l 0.001 -c chr1 -s 1234").split()
        return capture_output(stdpopsim.cli.stdpopsim_main, cmd)

    def test_script_generation(self):
        out, _ = self.docmd("--slim-script HomSap pop_0:5")
        assert "community.registerLateEvent" in out

        # msprime.MassMigration demographic events, with proportion<1.0
        # low level migration
        out, _ = self.docmd("--slim-script HomSap -d AncientEurasia_9K19 Mbuti:5")
        assert "community.registerLateEvent" in out
        # simultaneous mass migrations, with proportions summing to 1.0
        out, _ = self.docmd("--slim-script HomSap -d AmericanAdmixture_4B11 AFR:5")
        assert "community.registerLateEvent" in out

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    @pytest.mark.filterwarnings("ignore:.*has only.*individuals alive")
    def test_simulate(self):
        saved_slim_env = os.environ.get("SLIM")
        slim_path = os.environ.get("SLIM", "slim")
        with tempfile.NamedTemporaryFile(mode="w") as f:
            self.docmd(
                f"--slim-scaling-factor 20 --slim-path {slim_path} HomSap "
                f"pop_0:5 -o {f.name}"
            )
            ts = tskit.load(f.name)
        assert ts.num_samples == 10
        assert all(tree.num_roots == 1 for tree in ts.trees())

        if saved_slim_env is None:
            del os.environ["SLIM"]
        else:
            os.environ["SLIM"] = saved_slim_env

        with tempfile.NamedTemporaryFile(mode="w") as f:
            self.docmd(f"--slim-scaling-factor 20 HomSap pop_0:5 -o {f.name}")
            ts = tskit.load(f.name)
        assert ts.num_samples == 10

        # verify sample counts for a multipopulation demographic model
        with tempfile.NamedTemporaryFile(mode="w") as f:
            cmd = (
                "-q -e slim --slim-scaling-factor 20 --slim-burn-in 0 "
                f"HomSap -o {f.name} -l 0.001 -c chr1 -s 1234 "
                "-d OutOfAfrica_3G09 YRI:0 CEU:0 CHB:4"
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


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class TestWarningsAndErrors:
    """
    Checks that warning messages are printed when appropriate.
    """

    def triplet_diploid(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr22", length_multiplier=0.001)
        return engine, species, contig

    def triplet_haploid(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("EscCol")
        contig = species.get_contig("Chromosome", length_multiplier=0.001)
        return engine, species, contig

    @pytest.mark.filterwarnings("error::stdpopsim.SLiMOddSampleWarning")
    @pytest.mark.filterwarnings(
        "ignore:.*model has mutation rate.*but this simulation used.*"
    )
    def test_no_odd_sample_warning_for_even_samples(self):
        engine, species, contig = self.triplet_haploid()
        model = stdpopsim.PiecewiseConstantSize(100)
        samples = {"pop_0": 4}
        engine.simulate(
            demographic_model=model,
            contig=contig,
            samples=samples,
            dry_run=True,
        )

    def test_odd_sample_warning(self):
        engine, species, contig = self.triplet_haploid()
        model = stdpopsim.PiecewiseConstantSize(100)
        samples = {"pop_0": 5}
        with pytest.warns(stdpopsim.SLiMOddSampleWarning):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                dry_run=True,
            )

        model = stdpopsim.IsolationWithMigration(100, 100, 100, 100, 0.1, 0.1)
        samples = {"pop1": 2, "pop2": 5}
        with pytest.warns(stdpopsim.SLiMOddSampleWarning):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                dry_run=True,
            )

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_bad_population_size_addSubPop_warning(self):
        engine, species, contig = self.triplet_diploid()
        model = stdpopsim.PiecewiseConstantSize(100)
        samples = {"pop_0": 1}

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
        engine, species, contig = self.triplet_diploid()
        model = stdpopsim.PiecewiseConstantSize(100)
        samples = {"pop_0": 1}

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
        engine, species, contig = self.triplet_diploid()
        model = stdpopsim.IsolationWithMigration(
            NA=1000, N1=100, N2=1000, T=1000, M12=0, M21=0
        )
        samples = {"pop1": 1}
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
        engine, species, contig = self.triplet_diploid()
        model = stdpopsim.IsolationWithMigration(
            NA=1000, N1=100, N2=1000, T=1000, M12=0, M21=0
        )
        samples = {"pop1": 1}
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
        engine, species, contig = self.triplet_diploid()
        model = stdpopsim.PiecewiseConstantSize(100, (1000, 1000))
        samples = {"pop_0": 1}
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
        engine, species, contig = self.triplet_diploid()
        model = stdpopsim.PiecewiseConstantSize(100, (1000, 1000))
        samples = {"pop_0": 1}
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
        engine, species, contig = self.triplet_diploid()
        model = stdpopsim.PiecewiseConstantSize(1000)
        samples = {"pop_0": 150}

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
        pop_0 = stdpopsim.models.Population(id="pop_0", description="")
        return stdpopsim.DemographicModel(
            id="exp_decline",
            description="exp_decline",
            long_description="exp_decline",
            populations=[pop_0],
            generation_time=1,
            population_configurations=[
                msprime.PopulationConfiguration(
                    initial_size=N0,
                    growth_rate=r,
                    metadata=pop_0.asdict(),
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
        engine, species, contig = self.triplet_diploid()
        model = self.exp_decline()
        samples = {"pop_0": 1}
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
        engine, species, contig = self.triplet_diploid()
        model = self.exp_decline()
        samples = {"pop_0": 1}
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
        engine, species, contig = self.triplet_diploid()
        model = self.exp_decline()
        samples = {"pop_0": 15}

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
        engine, species, contig = self.triplet_diploid()
        model = stdpopsim.PiecewiseConstantSize(10000)
        samples = {"pop_0": 50}
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
        engine, species, contig = self.triplet_diploid()
        model = stdpopsim.PiecewiseConstantSize(10000)
        samples = {"pop_0": 50}
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
    def test_chr1(self, Q):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr1", genetic_map="HapMapII_GRCh38")
        model = stdpopsim.PiecewiseConstantSize(100)
        samples = {"pop_0": 5}
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
        samples = {"pop_0": 5}
        ts = engine.simulate(
            demographic_model=model,
            contig=contig,
            samples=samples,
            slim_burn_in=0.1,
            verbosity=3,
        )
        self.verify_recombination_map(contig, ts)
        assert list(ts.breakpoints()) == [0.0, midpoint, contig.length]

    def test_not_simulated_outside_region(self):
        # test that when left, right are specified
        # we legit don't simulate anything outside that region
        species = stdpopsim.get_species("AraTha")
        model = stdpopsim.PiecewiseConstantSize(100)
        samples = {"pop_0": 50}

        left, right = 100000, 900000
        contig = species.get_contig("1", left=left, right=right)
        dfe = species.get_dfe("Gamma_H18")
        exons = species.get_annotations("araport_11_exons")
        exon_intervals = exons.get_chromosome_annotations("1")
        contig.add_dfe(intervals=exon_intervals, DFE=dfe)

        engine = stdpopsim.get_engine("slim")
        ts = engine.simulate(
            model,
            contig,
            samples,
            seed=236,
            slim_burn_in=100,
        )

        assert ts.sequence_length > right
        assert ts.num_sites > 0
        assert left in ts.breakpoints()
        assert right in ts.breakpoints()
        assert left <= min(ts.sites_position)
        assert max(ts.sites_position) < right
        assert ts.num_trees > 2
        for t in ts.trees(root_threshold=2):
            tl, tr = t.interval
            if tl > right or tr < left:
                assert t.num_roots == 0


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class TestGenomicElementTypes(PiecewiseConstantSizeMixin):

    mut_params = {
        "f": ([-0.1], [0]),
        "g": ([-0.1, 0.1], [1, 2]),
        "e": ([0.1],),
        "n": ([-0.1, 0.2],),
        "w": ([0.1, 0.4],),
        "lp": ([-5, 0.2],),
        "ln": ([-5, 0.2],),
    }
    example_mut_types = [("f", stdpopsim.MutationType())] + [
        (t, stdpopsim.MutationType(distribution_type=t, distribution_args=p))
        for t, params in mut_params.items()
        for p in params
    ]
    for dcl, dcb in [([0.0, 1.0], [0.0]), ([-0.2, 1.4, 0.5], [-0.1, 0.1])]:
        for t in ("f", "e", "n"):
            mt = stdpopsim.MutationType(
                distribution_type=t,
                distribution_args=mut_params[t][0],
                dominance_coeff_list=dcl,
                dominance_coeff_breaks=dcb,
            )
            example_mut_types.append((t, mt))
    assert len(example_mut_types) == 16

    def get_example_dfes(self):
        # this is in a function because scoping is weird
        example_dfes = [
            stdpopsim.dfe.neutral_dfe(),
            stdpopsim.DFE(
                id="simple",
                description="just one mut type",
                long_description="🐺 🦊 🐕 🦝 🦮 🐩",
                mutation_types=[self.example_mut_types[4][1]],
            ),
            stdpopsim.DFE(
                id="less_simple",
                description="two mutation types",
                long_description="🦕 🦖 🐳 🐋 🐬 🦭🐠 🐡 🦈 🐙",
                proportions=[0.5, 0.5],
                mutation_types=[m for _, m in self.example_mut_types[2:4]],
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
                mutation_types=[m for _, m in self.example_mut_types],
            ),
        ]
        return example_dfes

    def test_slim_simulates_dfes(self):
        contig = get_test_contig()
        # make theta = 40
        contig.mutation_rate = 20 / (1000 * contig.recombination_map.sequence_length)
        dfes = [
            stdpopsim.DFE(
                id=str(j),
                description=f"mut_{j}",
                long_description=f"mut_{j}",
                mutation_types=[mt],
            )
            for j, (_, mt) in enumerate(self.example_mut_types)
        ]
        n = len(dfes)
        breaks = np.linspace(0, contig.length, n + 1)
        for j, dfe in enumerate(dfes):
            interval = [breaks[j], breaks[j + 1]]
            print(j, interval)
            contig.add_dfe(intervals=np.array([interval], dtype="int"), DFE=dfe)
        engine = stdpopsim.get_engine("slim")
        ts = engine.simulate(
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            seed=124,
            verbosity=3,
        )
        for j, (t, mt) in enumerate(self.example_mut_types):
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
                    if t == "lp":
                        assert md["selection_coeff"] > 0
                    elif t == "ln":
                        assert md["selection_coeff"] < 0

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
            # neutral. SLiM does not allow a genomicElementType to have all 0.0
            # proportions.
            prop_list = []
            for prop, mt in zip(dfe["proportions"], dfe["mutation_types"]):
                first_mt = True
                for i in mt["slim_mutation_type_id"]:
                    nonzero_prop = (
                        (not mt["is_neutral"]) or is_dfe_neutral
                    ) and first_mt
                    prop_list.append(prop if nonzero_prop else 0.0)
                    first_mt = False
            return prop_list

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
                    mid
                    for mt in dfe["mutation_types"]
                    for mid in mt["slim_mutation_type_id"]
                ],
            }
            assert ge == ge_from_meta

    def verify_discretized_dominance(self, contig, ts):
        # check that the dummy first mutation type of discretized
        # h-s relationship mutation types are absent
        mut_type_counts = collections.Counter(
            x["mutation_type"]
            for m in ts.mutations()
            for x in m.metadata["mutation_list"]
        )
        metadata_ids = [x["id"] for x in ts.metadata["stdpopsim"]["DFEs"]]
        slim_mt_info = ts.metadata["SLiM"]["user_metadata"]["mutationTypes"][0]
        has_recap = metadata_ids[-1] == "recapitation"
        slim_to_mt_map = {}
        assert len(contig.dfe_list) + has_recap == len(ts.metadata["stdpopsim"]["DFEs"])
        for dfe, ts_metadata in zip(contig.dfe_list, ts.metadata["stdpopsim"]["DFEs"]):
            assert dfe.id == ts_metadata["id"]
            assert len(dfe.mutation_types) == len(ts_metadata["mutation_types"])
            for mt, md in zip(dfe.mutation_types, ts_metadata["mutation_types"]):
                if mt.dominance_coeff_list is not None:
                    assert len(mt.dominance_coeff_list) + 1 == len(
                        md["slim_mutation_type_id"]
                    )
                    first_type = md["slim_mutation_type_id"][0]
                    assert mut_type_counts[first_type] == 0
                    for h, k in zip(
                        mt.dominance_coeff_list, md["slim_mutation_type_id"][1:]
                    ):
                        assert str(k) in slim_mt_info
                        assert np.allclose(slim_mt_info[str(k)][0]["dominanceCoeff"], h)
                        slim_to_mt_map[k] = mt

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
        self.verify_discretized_dominance(contig, ts)
        assert len(contig.dfe_list) == len(ge_types)
        for j, (dfe, intervals) in enumerate(
            zip(contig.dfe_list, contig.interval_list)
        ):
            assert str(j) in ge_types
            ge = self.slim_metadata_key0(ge_types, str(j))
            # checking that the neutral mutations have 0.0 proportion unless
            # all the mutations are neutral in this dfe
            assert len(dfe.mutation_types) == len(dfe.proportions)
            ge_index = 0
            for mt, dfe_prop in zip(dfe.mutation_types, dfe.proportions):
                slim_prop = ge["mutationFractions"][ge_index]
                ge_index += 1
                if mt.is_neutral and not dfe.is_neutral:
                    assert slim_prop == 0.0
                else:
                    assert slim_prop == dfe_prop
                if mt.dominance_coeff_list is not None:
                    for _ in mt.dominance_coeff_list:
                        slim_prop = ge["mutationFractions"][ge_index]
                        assert slim_prop == 0.0
                        ge_index += 1
            assert ge_index == len(ge["mutationFractions"])
            # "+1" because SLiM's intervals are closed on both ends,
            # stdpopsim's are closed on the left, open on the right
            slim_intervals = np.column_stack(
                [
                    ge["intervalStarts"],
                    [x + 1 for x in ge["intervalEnds"]],
                ]
            )
            assert_array_equal(intervals, slim_intervals)
            ge_index = 0
            for mt in dfe.mutation_types:
                mt_id = ge["mutationTypes"][ge_index]
                ge_index += 1
                assert str(mt_id) in mut_types
                slim_mt = self.slim_metadata_key0(mut_types, str(mt_id))
                if mt.dominance_coeff_list is None:
                    assert mt.dominance_coeff == self.slim_metadata_key0(
                        slim_mt, "dominanceCoeff"
                    )
                else:
                    assert 0.5 == self.slim_metadata_key0(slim_mt, "dominanceCoeff")
                    for h in mt.dominance_coeff_list:
                        mt_id = ge["mutationTypes"][ge_index]
                        ge_index += 1
                        assert str(mt_id) in mut_types
                        slim_mt = self.slim_metadata_key0(mut_types, str(mt_id))
                        assert np.allclose(
                            h, self.slim_metadata_key0(slim_mt, "dominanceCoeff")
                        )
                assert mt.distribution_type == self.slim_metadata_key0(
                    slim_mt, "distributionType"
                )
                assert len(mt.distribution_args) == len(slim_mt["distributionParams"])
                for a, b in zip(mt.distribution_args, slim_mt["distributionParams"]):
                    assert a == b
            assert ge_index == len(ge["mutationTypes"])

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
            mutation_types=[self.example_mut_types[5][1]],
        )
        dfe1 = stdpopsim.DFE(
            id="dfe",
            description="the second one",
            long_description="I'm different but have the same name! =( =(",
            proportions=[0.2, 0.8],
            mutation_types=[m for _, m in self.example_mut_types[7:9]],
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
                    distribution_type="lp", distribution_args=[0.01, 0.2]
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

    def test_warn_when_dfe_intervals_outside_contig(self):
        contig = get_test_contig()
        dfe = stdpopsim.DFE(
            id="test",
            description="non-neutral",
            long_description="",
            mutation_types=[
                stdpopsim.MutationType(),
                stdpopsim.MutationType(
                    distribution_type="lp", distribution_args=[0.01, 0.2]
                ),
            ],
            proportions=[0.5, 0.5],
        )
        with pytest.warns(UserWarning, match="No intervals remain"):
            contig.add_dfe(np.array([[1e6, 1.1e6]], dtype="int"), dfe)
        assert contig.interval_list[1].shape[0] == 0

    def test_chromosomal_segment(self):
        left = 100101
        right = 201024
        contig = get_test_contig(length_multiplier=1, left=left, right=right)
        L = contig.length
        example_dfes = self.get_example_dfes()
        for j, dfe in enumerate(example_dfes):
            contig.add_dfe(
                left
                + np.array([[10 * (j + 1), 100 * (10 - j)], [L / 2, L]], dtype="int"),
                dfe,
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
        assert int(ts.sequence_length) == contig.length

    def test_dominance_coeff_list(self):
        # test that discretized h-s relationships work
        contig = get_test_contig()
        L = int(contig.length)
        emts = [emt for _, emt in self.example_mut_types[10:12]]
        for emt in emts:
            assert emt.dominance_coeff_list is not None
        dfe0 = stdpopsim.DFE(
            id="dfe",
            description="the first one",
            long_description="hello world",
            proportions=[1.0, 0.0],
            mutation_types=emts,
        )
        emts = [emt for _, emt in self.example_mut_types[12:15]]
        for emt in emts:
            assert emt.dominance_coeff_list is not None
        dfe1 = stdpopsim.DFE(
            id="dfe",
            description="the second one",
            long_description="I'm different but have the same name! =( =(",
            proportions=[0.2, 0.7, 0.1],
            mutation_types=emts,
        )
        contig.add_dfe(np.array([[0, 0.5 * L]], dtype="int"), dfe0)
        contig.add_dfe(np.array([[0.2 * L, 0.4 * L]], dtype="int"), dfe0)
        contig.add_dfe(np.array([[0.45 * L, L]], dtype="int"), dfe1)
        contig.add_dfe(np.array([[0.7 * L, 0.9 * L]], dtype="int"), dfe0)
        assert len(contig.dfe_list) == 5
        engine = stdpopsim.get_engine("slim")
        contig.mutation_rate *= 10
        ts = engine.simulate(
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            verbosity=3,  # to get metadata output
        )
        assert len(ts.metadata["stdpopsim"]["DFEs"]) == len(contig.dfe_list) + 1
        # slim mutation type IDs with dominance coeff lists:
        mut_id_haslist = {}
        for dfe in ts.metadata["stdpopsim"]["DFEs"]:
            for mt in dfe["mutation_types"]:
                haslist = (
                    "dominance_coeff_list" in mt
                    and mt["dominance_coeff_list"] is not None
                )
                for slim_id in mt["slim_mutation_type_id"]:
                    assert slim_id not in mut_id_haslist
                    mut_id_haslist[slim_id] = haslist
        num_target_muts = 0
        for mut in ts.mutations():
            for md in mut.metadata["mutation_list"]:
                if mut_id_haslist[md["mutation_type"]]:
                    num_target_muts += 1
        # the number 20 is not important, just want to make sure we have *some*
        assert num_target_muts > 20
        self.verify_genomic_elements(contig, ts)
        self.verify_mutation_rates(contig, ts)


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
        assert header[0] == "tick"
        assert header[1][:8] == "fitness_"
        assert header[2][:8] == "fitness_"
        # neutral model, should have no fitness variation
        assert np.all(data[:, 1] == 1.0)
        assert np.all(data[:, 2] == 0.0)


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class TestDrawMutation(PiecewiseConstantSizeMixin):
    def test_draw_mutation(self):
        contig = get_test_contig()
        contig.add_single_site(id=self.mut_id, coordinate=100)
        extended_events = [
            stdpopsim.DrawMutation(
                time=self.T_mut,
                single_site_id=self.mut_id,
                population="pop_0",
            ),
        ]
        engine = stdpopsim.get_engine("slim")
        engine.simulate(
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            extended_events=extended_events,
            dry_run=True,
        )

    def test_draw_mutations_at_different_sites(self):
        contig = get_test_contig()
        contig.add_single_site(id="recent", coordinate=100)
        contig.add_single_site(id="older", coordinate=101)
        extended_events = [
            stdpopsim.DrawMutation(
                time=self.T_mut // 2,
                single_site_id="recent",
                population="pop_0",
            ),
            stdpopsim.DrawMutation(
                time=self.T_mut,
                single_site_id="older",
                population="pop_0",
            ),
        ]
        engine = stdpopsim.get_engine("slim")
        engine.simulate(
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            extended_events=extended_events,
            dry_run=True,
        )

    def test_draw_multiple_mutations_at_same_site(self):
        contig = get_test_contig()
        contig.add_single_site(id="mutant", coordinate=100)
        extended_events = [
            stdpopsim.DrawMutation(
                time=self.T_mut // 2,
                single_site_id="mutant",
                population="pop_0",
            ),
            stdpopsim.DrawMutation(
                time=self.T_mut,
                single_site_id="mutant",
                population="pop_0",
            ),
        ]
        engine = stdpopsim.get_engine("slim")
        with pytest.raises(ValueError, match="maximum of one mutation is allowed"):
            engine.simulate(
                demographic_model=self.model,
                contig=contig,
                samples=self.samples,
                extended_events=extended_events,
                dry_run=True,
            )

    def test_invalid_single_site_id(self):
        engine = stdpopsim.get_engine("slim")
        for single_site_id in ["deleterious", "sweep"]:
            extended_events = [
                stdpopsim.DrawMutation(
                    time=self.T_mut,
                    single_site_id=single_site_id,
                    population="pop_0",
                ),
            ]
            with pytest.raises(ValueError, match="must exist and be uniquely labelled"):
                engine.simulate(
                    demographic_model=self.model,
                    contig=self.contig,
                    samples=self.samples,
                    extended_events=extended_events,
                    dry_run=True,
                )

    def test_no_mutation_types_defined(self):
        extended_events = [
            stdpopsim.DrawMutation(
                time=self.T_mut,
                single_site_id=self.mut_id,
                population="pop_0",
            ),
        ]
        contig = get_test_contig()
        contig.add_single_site(id=self.mut_id, coordinate=100)
        contig.dfe_list[1].mutation_types = []
        engine = stdpopsim.get_engine("slim")
        with pytest.raises(ValueError, match="must contain a single mutation type"):
            engine.simulate(
                demographic_model=self.model,
                contig=contig,
                samples=self.samples,
                extended_events=extended_events,
                dry_run=True,
            )

    def test_multiple_mutation_types_defined(self):
        extended_events = [
            stdpopsim.DrawMutation(
                time=self.T_mut,
                single_site_id="test",
                population="pop_0",
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
        with pytest.raises(ValueError, match="must contain a single mutation type"):
            engine.simulate(
                demographic_model=self.model,
                contig=contig,
                samples=self.samples,
                extended_events=extended_events,
                dry_run=True,
            )

    def test_fitness_distribution_not_fixed(self):
        extended_events = [
            stdpopsim.DrawMutation(
                time=self.T_mut,
                single_site_id="test",
                population="pop_0",
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
        with pytest.raises(ValueError, match="instead of a fixed fitness coefficient"):
            engine.simulate(
                demographic_model=self.model,
                contig=contig,
                samples=self.samples,
                extended_events=extended_events,
                dry_run=True,
            )

    def test_multiple_intervals(self):
        extended_events = [
            stdpopsim.DrawMutation(
                time=self.T_mut,
                single_site_id="test",
                population="pop_0",
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
        with pytest.raises(ValueError, match="refers to a DFE with intervals"):
            engine.simulate(
                demographic_model=self.model,
                contig=contig,
                samples=self.samples,
                extended_events=extended_events,
                dry_run=True,
            )

    def test_interval_too_large(self):
        extended_events = [
            stdpopsim.DrawMutation(
                time=self.T_mut,
                single_site_id="test",
                population="pop_0",
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
        with pytest.raises(ValueError, match="refers to a DFE with intervals"):
            engine.simulate(
                demographic_model=self.model,
                contig=contig,
                samples=self.samples,
                extended_events=extended_events,
                dry_run=True,
            )

    def test_duplicate_mutation_ids(self):
        extended_events = [
            stdpopsim.DrawMutation(
                time=self.T_mut,
                single_site_id="test",
                population="pop_0",
            ),
        ]
        contig = get_test_contig()
        contig.add_single_site(id="test", coordinate=100)
        contig.add_single_site(id="test", coordinate=110)
        engine = stdpopsim.get_engine("slim")
        with pytest.raises(ValueError, match="must exist and be uniquely labelled"):
            engine.simulate(
                demographic_model=self.model,
                contig=contig,
                samples=self.samples,
                extended_events=extended_events,
                dry_run=True,
            )

    def test_mutation_has_no_interval(self):
        extended_events = [
            stdpopsim.DrawMutation(
                time=self.T_mut,
                single_site_id="test",
                population="pop_0",
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
                stdpopsim.DrawMutation(
                    time=time,
                    single_site_id="irrelevant",
                    population="pop_0",
                )
        for time in (0, -1):
            with pytest.raises(ValueError):
                stdpopsim.DrawMutation(
                    time=stdpopsim.GenerationAfter(time),
                    single_site_id="irrelevant",
                    population="pop_0",
                )


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class TestAlleleFrequencyConditioning(PiecewiseConstantSizeMixin):
    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_drawn_mutation_not_lost(self):
        engine = stdpopsim.get_engine("slim")
        ct = get_test_contig()
        ct.mutation_rate = 0.0
        ct.add_single_site(id="test", coordinate=100)
        extended_events = [
            stdpopsim.DrawMutation(
                time=self.T_mut,
                single_site_id="test",
                population="pop_0",
            ),
            stdpopsim.ConditionOnAlleleFrequency(
                start_time=0,
                end_time=0,
                single_site_id="test",
                population="pop_0",
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
            stdpopsim.DrawMutation(
                time=self.T_mut,
                single_site_id="test",
                population="pop_0",
            ),
            stdpopsim.ConditionOnAlleleFrequency(
                start_time=0,
                end_time=0,
                single_site_id="test",
                population="pop_0",
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
                stdpopsim.DrawMutation(
                    time=self.T_mut,
                    single_site_id="test",
                    population="pop_0",
                ),
                # Condition on desired AF at end of simulation.
                stdpopsim.ConditionOnAlleleFrequency(
                    start_time=0,
                    end_time=0,
                    single_site_id="test",
                    population="pop_0",
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
                stdpopsim.ConditionOnAlleleFrequency(
                    start_time=0,
                    end_time=0,
                    single_site_id="irrelevant",
                    population="pop_0",
                    op=op,
                    allele_frequency=af,
                )

    def test_bad_times(self):
        for start_time, end_time in [(-1, 0), (0, -1), (1, 100)]:
            with pytest.raises(ValueError):
                stdpopsim.ConditionOnAlleleFrequency(
                    start_time=start_time,
                    end_time=end_time,
                    single_site_id="irrelevant",
                    population="pop_0",
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
                    stdpopsim.DrawMutation(
                        time=self.T_mut,
                        single_site_id=self.mut_id,
                        population="pop_0",
                    ),
                    stdpopsim.ConditionOnAlleleFrequency(
                        start_time=stdpopsim.GenerationAfter(start_time),
                        end_time=end_time,
                        single_site_id=self.mut_id,
                        population="pop_0",
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
        op_types = stdpopsim.ConditionOnAlleleFrequency.op_types
        for op in op_types:
            id = stdpopsim.ConditionOnAlleleFrequency.op_id(op)
            assert 0 <= id < len(op_types)
        for op in ("==", "=", "!=", {}, ""):
            with pytest.raises(ValueError):
                id = stdpopsim.ConditionOnAlleleFrequency.op_id(op)

    def test_no_drawn_mutation(self):
        extended_events = [
            stdpopsim.ConditionOnAlleleFrequency(
                start_time=stdpopsim.GenerationAfter(self.T_mut),
                end_time=0,
                single_site_id=self.mut_id,
                population="pop_0",
                op=">",
                allele_frequency=0,
            ),
        ]
        engine = stdpopsim.get_engine("slim")
        with pytest.raises(ValueError, match="no mutation is drawn at this site"):
            engine.simulate(
                demographic_model=self.model,
                contig=self.contig,
                samples=self.samples,
                extended_events=extended_events,
                dry_run=True,
            )


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class TestChangeMutationFitness(PiecewiseConstantSizeMixin):
    # Testing stdpopsim.ChangeMutationFitness is challenging, because
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
                stdpopsim.DrawMutation(
                    time=self.T_mut,
                    single_site_id=self.mut_id,
                    population="pop_0",
                ),
                stdpopsim.ChangeMutationFitness(
                    start_time=stdpopsim.GenerationAfter(self.T_mut),
                    end_time=0,
                    single_site_id=self.mut_id,
                    population="pop_0",
                    selection_coeff=0.1,
                    dominance_coeff=0.5,
                ),
                # Condition on AF > 0, to restore() immediately if the
                # allele is lost.
                stdpopsim.ConditionOnAlleleFrequency(
                    start_time=stdpopsim.GenerationAfter(self.T_mut),
                    end_time=0,
                    single_site_id=self.mut_id,
                    population="pop_0",
                    op=">",
                    allele_frequency=0,
                ),
                # Condition on desired AF at end of simulation.
                stdpopsim.ConditionOnAlleleFrequency(
                    start_time=0,
                    end_time=0,
                    single_site_id=self.mut_id,
                    population="pop_0",
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

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_referenced_single_site_is_nonneutral(self):
        engine = stdpopsim.get_engine("slim")
        contig = get_test_contig()
        contig.add_single_site("one", coordinate=100)
        contig.add_single_site("two", coordinate=101)
        extended_events = [
            stdpopsim.DrawMutation(
                time=self.T_mut,
                single_site_id="one",
                population="pop_0",
            ),
        ]
        ts = engine.simulate(
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            extended_events=extended_events,
        )
        referenced_dfe = ts.metadata["stdpopsim"]["DFEs"][1]
        assert referenced_dfe["id"] == "one"
        assert referenced_dfe["mutation_types"][0]["is_neutral"] is False
        unreferenced_dfe = ts.metadata["stdpopsim"]["DFEs"][2]
        assert unreferenced_dfe["id"] == "two"
        assert unreferenced_dfe["mutation_types"][0]["is_neutral"] is True

    def test_no_drawn_mutation(self):
        extended_events = [
            stdpopsim.ChangeMutationFitness(
                start_time=stdpopsim.GenerationAfter(self.T_mut),
                end_time=0,
                single_site_id=self.mut_id,
                selection_coeff=0.1,
                dominance_coeff=0.5,
            ),
        ]
        engine = stdpopsim.get_engine("slim")
        with pytest.raises(ValueError, match="no mutation is drawn at this site"):
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
                stdpopsim.ChangeMutationFitness(
                    start_time=start_time,
                    end_time=end_time,
                    single_site_id="irrelevant",
                    population="pop_0",
                    selection_coeff=0.1,
                    dominance_coeff=0.5,
                )

    def test_population_name(self):
        engine = stdpopsim.get_engine("slim")
        for pop in ["not present", 0]:
            extended_events = [
                stdpopsim.ChangeMutationFitness(
                    start_time=stdpopsim.GenerationAfter(self.T_mut),
                    end_time=0,
                    single_site_id=self.mut_id,
                    population=pop,
                    selection_coeff=0.1,
                    dominance_coeff=0.5,
                ),
            ]
            with pytest.raises(ValueError, match="is not in demographic model"):
                engine.simulate(
                    demographic_model=self.model,
                    contig=self.contig,
                    samples=self.samples,
                    extended_events=extended_events,
                    dry_run=True,
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
                    stdpopsim.DrawMutation(
                        time=self.T_mut,
                        single_site_id=self.mut_id,
                        population="pop_0",
                    ),
                    stdpopsim.ChangeMutationFitness(
                        start_time=stdpopsim.GenerationAfter(start_time),
                        end_time=end_time,
                        single_site_id=self.mut_id,
                        population="pop_0",
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


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
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
    def test_sweep(self, tmp_path):
        mutation_generation_ago = 1000
        for start_generation_ago in [
            mutation_generation_ago // 2,
            mutation_generation_ago,
        ]:
            for end_generation_ago in [0, mutation_generation_ago // 4]:
                logfile = tmp_path / "sweep.logfile"
                engine = stdpopsim.get_engine("slim")
                contig = get_test_contig()
                locus_id = "sweep"
                contig.add_single_site(
                    id=locus_id,
                    coordinate=100,
                )
                extended_events = stdpopsim.selective_sweep(
                    single_site_id=locus_id,
                    population="pop_0",
                    mutation_generation_ago=mutation_generation_ago,
                    start_generation_ago=start_generation_ago,
                    end_generation_ago=end_generation_ago,
                    selection_coeff=0.1,
                )
                engine.simulate(
                    demographic_model=self.model,
                    contig=contig,
                    samples=self.samples,
                    extended_events=extended_events,
                    slim_burn_in=1,
                    logfile=logfile,
                    logfile_interval=1,
                )
                in_sweep, outside_sweep, _ = self._fitness_per_generation(
                    logfile=logfile,
                    start_generation_ago=start_generation_ago,
                    end_generation_ago=end_generation_ago,
                    pop=0,
                )
                assert np.all(in_sweep > 1)
                assert np.all(outside_sweep == 1)

    @pytest.mark.usefixtures("tmp_path")
    def test_sweep_meets_min_freq_at_start(self, tmp_path):
        mutation_generation_ago = 10
        start_generation_ago = 8
        s = 0.01
        # condition on a difficult-to-reach frequency, i.e.
        min_freq = 3 * (mutation_generation_ago - start_generation_ago) / (2 * self.N0)
        logfile = tmp_path / "sweep_start_af.logfile"
        engine = stdpopsim.get_engine("slim")
        contig = get_test_contig()
        locus_id = "sweep"
        contig.add_single_site(
            id=locus_id,
            coordinate=100,
        )
        extended_events = stdpopsim.selective_sweep(
            single_site_id=locus_id,
            population="pop_0",
            mutation_generation_ago=mutation_generation_ago,
            start_generation_ago=start_generation_ago,
            min_freq_at_start=min_freq,
            selection_coeff=s,
            dominance_coeff=0.5,
        )
        while True:
            engine.simulate(
                demographic_model=self.model,
                contig=contig,
                samples=self.samples,
                extended_events=extended_events,
                slim_burn_in=1,
                logfile=logfile,
                logfile_interval=1,
            )
            in_sweep, outside_sweep, rejections = self._fitness_per_generation(
                logfile=logfile,
                start_generation_ago=start_generation_ago,
                end_generation_ago=0,
                pop=0,
            )
            # ensure that rejections are occuring in the generation of the AF
            # condition
            if start_generation_ago in rejections.keys():
                break
        assert np.all(outside_sweep == 1.0)
        assert in_sweep[0] >= s * min_freq + 1

    @pytest.mark.usefixtures("tmp_path")
    def test_sweep_meets_min_freq_at_end(self, tmp_path):
        mutation_generation_ago = 10
        end_generation_ago = 8
        s = 0.01
        # condition on a difficult-to-reach frequency, i.e.
        min_freq = 3 * (mutation_generation_ago - end_generation_ago) / (2 * self.N0)
        logfile = tmp_path / "sweep_end_af.logfile"
        engine = stdpopsim.get_engine("slim")
        contig = get_test_contig()
        locus_id = "sweep"
        contig.add_single_site(
            id=locus_id,
            coordinate=100,
        )
        extended_events = stdpopsim.selective_sweep(
            single_site_id=locus_id,
            population="pop_0",
            mutation_generation_ago=mutation_generation_ago,
            end_generation_ago=end_generation_ago,
            min_freq_at_end=min_freq,
            selection_coeff=s,
            dominance_coeff=0.5,
        )
        while True:
            engine.simulate(
                demographic_model=self.model,
                contig=contig,
                samples=self.samples,
                extended_events=extended_events,
                slim_burn_in=1,
                logfile=logfile,
                logfile_interval=1,
            )
            in_sweep, outside_sweep, rejections = self._fitness_per_generation(
                logfile=logfile,
                start_generation_ago=mutation_generation_ago,
                end_generation_ago=end_generation_ago,
                pop=0,
            )
            # ensure that rejections are occuring in the generation of the AF
            # condition
            if end_generation_ago in rejections.keys():
                break
        assert np.all(outside_sweep == 1.0)
        assert in_sweep[-1] >= s * min_freq + 1

    def test_sweep_with_negative_selection_coeff(self):
        with pytest.raises(ValueError, match="coefficient must be"):
            stdpopsim.selective_sweep(
                single_site_id="irrelevant",
                population="irrelevant",
                mutation_generation_ago=1000,
                selection_coeff=-0.1,
            )

    def test_sweep_with_bad_AF_conditions(self):
        for start_freq, end_freq in zip([-0.1, 0.1], [0.1, -0.1]):
            with pytest.raises(ValueError, match="of the sweep must be in"):
                stdpopsim.selective_sweep(
                    single_site_id="irrelevant",
                    population="irrelevant",
                    mutation_generation_ago=1000,
                    selection_coeff=0.1,
                    start_generation_ago=100,
                    min_freq_at_start=start_freq,
                    min_freq_at_end=end_freq,
                )
        with pytest.raises(
            ValueError, match="coincides with the introduction of the mutation"
        ):
            stdpopsim.selective_sweep(
                single_site_id="irrelevant",
                population="irrelevant",
                mutation_generation_ago=1000,
                start_generation_ago=1000,
                selection_coeff=0.1,
                min_freq_at_start=0.1,
            )

    @pytest.mark.usefixtures("tmp_path")
    def test_global_sweep(self, tmp_path):
        mutation_generation_ago = 100
        start_generation_ago = 50
        end_generation_ago = 30
        logfile = tmp_path / "sweep_not_restricted.logfile"
        engine = stdpopsim.get_engine("slim")
        model = self._get_island_model()
        contig = get_test_contig()
        locus_id = "sweep"
        contig.add_single_site(
            id=locus_id,
            coordinate=100,
        )
        extended_events = stdpopsim.selective_sweep(
            single_site_id=locus_id,
            population="pop_1",
            mutation_generation_ago=mutation_generation_ago,
            start_generation_ago=start_generation_ago,
            end_generation_ago=end_generation_ago,
            selection_coeff=1000.0,  # so that mutation stays in non-source pop
            globally_adaptive=True,
        )
        while True:
            # The while loop is to avoid rare stochastic failure due to
            # mutation never migrating into the non-source population
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=self.samples,
                extended_events=extended_events,
                slim_burn_in=1,
                logfile=logfile,
                logfile_interval=1,
            )
            p0_in_sweep, p0_outside_sweep, _ = self._fitness_per_generation(
                logfile=logfile,
                start_generation_ago=start_generation_ago,
                end_generation_ago=end_generation_ago,
                pop=0,
            )
            # Mutation is in non-source pop at start of sweep
            if p0_in_sweep[0] > 1:
                break
        p1_in_sweep, p1_outside_sweep, _ = self._fitness_per_generation(
            logfile=logfile,
            start_generation_ago=start_generation_ago,
            end_generation_ago=end_generation_ago,
            pop=1,
        )
        assert np.all(p0_in_sweep > 1)
        assert np.all(p0_outside_sweep == 1)
        assert np.all(p1_in_sweep > 1)
        assert np.all(p1_outside_sweep == 1)

    @pytest.mark.usefixtures("tmp_path")
    def test_local_sweep(self, tmp_path):
        mutation_generation_ago = 100
        start_generation_ago = 50
        end_generation_ago = 30
        logfile = tmp_path / "sweep_restricted.logfile"
        engine = stdpopsim.get_engine("slim")
        model = self._get_island_model()
        contig = get_test_contig()
        locus_id = "sweep"
        contig.add_single_site(
            id=locus_id,
            coordinate=100,
        )
        extended_events = stdpopsim.selective_sweep(
            single_site_id=locus_id,
            population="pop_1",
            mutation_generation_ago=mutation_generation_ago,
            start_generation_ago=start_generation_ago,
            end_generation_ago=end_generation_ago,
            selection_coeff=0.1,
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
        )
        p0_in_sweep, p0_outside_sweep, _ = self._fitness_per_generation(
            logfile=logfile,
            start_generation_ago=start_generation_ago,
            end_generation_ago=end_generation_ago,
            pop=0,
        )
        p1_in_sweep, p1_outside_sweep, _ = self._fitness_per_generation(
            logfile=logfile,
            start_generation_ago=start_generation_ago,
            end_generation_ago=end_generation_ago,
            pop=1,
        )
        assert np.all(p0_in_sweep == 1)
        assert np.all(p0_outside_sweep == 1)
        assert np.all(p1_in_sweep > 1)
        assert np.all(p1_outside_sweep == 1)

    @pytest.mark.usefixtures("tmp_path")
    def test_sweeps_at_multiple_sites(self, tmp_path):
        logfile = tmp_path / "sweep_multi.logfile"
        engine = stdpopsim.get_engine("slim")
        contig = get_test_contig()
        model = self._get_island_model()
        mutation_generation_ago = 100
        ids = ["sweep1", "sweep2"]
        coords = [100, contig.length - 100]
        start_generation_agos = [mutation_generation_ago, mutation_generation_ago // 2]
        end_generation_agos = [0, 0]
        pop_ids = ["pop_0", "pop_1"]
        extended_events = []
        for locus_id, coord, start_generation_ago, end_generation_ago, pop in zip(
            ids, coords, start_generation_agos, end_generation_agos, pop_ids
        ):
            contig.add_single_site(
                id=locus_id,
                coordinate=coord,
            )
            extended_events += stdpopsim.selective_sweep(
                single_site_id=locus_id,
                population=pop,
                mutation_generation_ago=mutation_generation_ago,
                start_generation_ago=start_generation_ago,
                end_generation_ago=end_generation_ago,
                selection_coeff=1.0,
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
        )
        for i, _ in enumerate(pop_ids):
            in_sweep, outside_sweep, _ = self._fitness_per_generation(
                logfile=logfile,
                start_generation_ago=start_generation_agos[i],
                end_generation_ago=end_generation_agos[i],
                pop=i,
            )
            assert np.all(in_sweep > 1)
            assert np.all(outside_sweep == 1)

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


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class TestSelectionCoeffFromMutation:

    species = stdpopsim.get_species("HomSap")
    model = stdpopsim.PiecewiseConstantSize(100)
    samples = {"pop_0": 5}
    dfe = stdpopsim.DFE(
        id="test",
        description="",
        long_description="",
        mutation_types=[
            stdpopsim.MutationType(distribution_type="f", distribution_args=[-0.01]),
            stdpopsim.MutationType(distribution_type="f", distribution_args=[0.0]),
        ],
        proportions=[0.5, 0.5],
    )

    def test_stacked(self):
        engine = stdpopsim.get_engine("slim")
        contig = self.species.get_contig(length=20, mutation_rate=1e-2)
        contig.add_dfe(np.array([[0, contig.length // 2]]), self.dfe)
        while True:
            ts = engine.simulate(
                demographic_model=self.model,
                contig=contig,
                samples=self.samples,
                slim_burn_in=10,
            )
            is_stacked = [len(m.metadata["mutation_list"]) > 1 for m in ts.mutations()]
            if any(is_stacked):
                break
        selection_coeffs = [
            stdpopsim.selection_coeff_from_mutation(ts, m) for m in ts.mutations()
        ]
        assert np.all(
            np.logical_or(
                np.isclose(selection_coeffs, -0.01),
                np.isclose(selection_coeffs, 0.0),
            )
        )

    def test_msprime(self):
        engine = stdpopsim.get_engine("msprime")
        contig = self.species.get_contig(length=20, mutation_rate=1e-2)
        while True:
            ts = engine.simulate(
                demographic_model=self.model,
                contig=contig,
                samples=self.samples,
            )
            if ts.num_mutations > 0:
                break
        selection_coeffs = [
            stdpopsim.selection_coeff_from_mutation(ts, m) for m in ts.mutations()
        ]
        assert np.allclose(selection_coeffs, 0.0)

    def test_errors(self):
        engine = stdpopsim.get_engine("msprime")
        contig = self.species.get_contig(length=20, mutation_rate=1e-2)
        while True:
            ts = engine.simulate(
                demographic_model=self.model,
                contig=contig,
                samples=self.samples,
            )
            if ts.num_mutations > 0:
                break
        with pytest.raises(ValueError, match="must be a"):
            stdpopsim.selection_coeff_from_mutation("foo", next(ts.mutations()))
        with pytest.raises(ValueError, match="must be a"):
            stdpopsim.selection_coeff_from_mutation(ts, "bar")


@pytest.mark.skipif(IS_WINDOWS, reason="SLiM not available on windows")
class TestPloidy:
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
            engine.simulate(model, contig, samples={"pop_0": 2}, verbosity=2)
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
            ts = engine.simulate(model, contig, samples={"pop_0": 2})
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
        ts = engine.simulate(model, contig, samples={"pop_0": 3})
        ts_hap = stdpopsim.utils.haploidize_individuals(ts)
        assert ts_hap.num_individuals == ts.num_individuals * 2
        for i, j in zip(ts.samples(), ts_hap.samples()):
            assert i == j
            ind_i = ts.nodes_individual[i]
            ind_j = ts_hap.nodes_individual[j]
            assert ts.tables.individuals[ind_i] == ts_hap.tables.individuals[ind_j]
