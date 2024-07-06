"""
Tests for SLiM simulation engine.
"""
import os
import io
import tempfile
from unittest import mock
import numpy as np
import collections
import re

import pytest
import tskit
import msprime

import stdpopsim
import stdpopsim.cli
from .test_cli import capture_output

slim_path = os.environ.get("SLIM", "slim")


def count_mut_types(ts):
    selection_coeffs = [
        stdpopsim.ext.selection_coeff_from_mutation(ts, mut) for mut in ts.mutations()
    ]
    num_neutral = sum([s == 0 for s in selection_coeffs])
    return [num_neutral, abs(len(selection_coeffs) - num_neutral)]


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

    @pytest.mark.filterwarnings("ignore:Recombination map has length:UserWarning")
    def test_recombination_map(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr21", genetic_map="HapMapII_GRCh37")
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
        assert all(tree.num_roots == 1 for tree in ts.trees())

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
            assert all(tree.num_roots == 1 for tree in ts.trees())

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
        out = re.sub(
            r'defineConstant\("trees_file.+;',
            rf'defineConstant("trees_file", "{treefile}");',
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
        assert all(tree.num_roots == 1 for tree in ts.trees())
        n_mut_types = count_mut_types(ts)
        assert n_mut_types[0] > 0
        assert n_mut_types[1] > 0

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    @pytest.mark.usefixtures("tmp_path")
    def test_dfe_no_demography(self, tmp_path):
        fname = tmp_path / "test_dfe_no_demography.trees"
        cmd = (
            f"-q -e slim --slim-scaling-factor 50 --slim-path {slim_path} "
            f"HomSap -c chr22 -l 0.02 -o {fname} --dfe Gamma_K17 -s 24 "
            f"pop_0:5"
        ).split()
        capture_output(stdpopsim.cli.stdpopsim_main, cmd)
        ts = tskit.load(fname)
        self.verify_slim_sim(ts, num_samples=10)

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    @pytest.mark.usefixtures("tmp_path")
    def test_dfe_interval(self, tmp_path):
        fname = tmp_path / "test_dfe_interval.trees"
        cmd = (
            f"-q -e slim --slim-scaling-factor 40 --slim-path {slim_path} "
            f"HomSap -c chr21 -l 0.01 -o {fname} --dfe Gamma_K17 -s 984 "
            f"--dfe-interval 1000,100000 pop_0:5"
        ).split()
        capture_output(stdpopsim.cli.stdpopsim_main, cmd)
        ts = tskit.load(fname)
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
    @pytest.mark.usefixtures("tmp_path")
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
            f"-q -e slim --slim-scaling-factor 30 --slim-path {slim_path} "
            f"HomSap -c chr20 -s 1234 -l 0.01 -o {fname} --dfe Gamma_K17 -s 183 "
            f"--dfe-bed-file {tmp_path / 'ex.bed'} pop_0:5"
        ).split()
        capture_output(stdpopsim.cli.stdpopsim_main, cmd)
        ts = tskit.load(fname)
        self.verify_slim_sim(ts, num_samples=10)

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_chromosomal_segment(self):
        left = 101024
        right = 200111
        with tempfile.NamedTemporaryFile(mode="w") as f:
            cmd = (
                f"-q -e slim --slim-scaling-factor 30 --slim-path {slim_path} "
                f"HomSap -c chr22 --left {left} --right {right} -o {f.name} "
                f"--dfe Gamma_K17 -s 24 pop_0:5"
            ).split()
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)
            ts = tskit.load(f.name)
        assert ts.sequence_length == right - left
        self.verify_slim_sim(ts, num_samples=10)

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_chromosomal_segment_with_dfe_interval(self):
        left = 101024
        right = 200111
        dfe_left = left - 5000
        dfe_right = left + 90000
        with tempfile.NamedTemporaryFile(mode="w") as f:
            cmd = (
                f"-q -e slim --slim-scaling-factor 30 --slim-path {slim_path} "
                f"HomSap -c chr22 --left {left} --right {right} -o {f.name} "
                f"--dfe Gamma_K17 -s 984 --dfe-interval {dfe_left},{dfe_right} "
                f"pop_0:5"
            ).split()
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)
            ts = tskit.load(f.name)
        assert ts.sequence_length == right - left
        self.verify_slim_sim(ts, num_samples=10)

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_chromosomal_segment_with_dfe_bed_file(self, tmp_path):
        left = 101024
        right = 300111
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
            f"-q -e slim --slim-scaling-factor 30 --slim-path {slim_path} "
            f"HomSap -c chr22 -s 1234 --left {left} --right {right} -o {fname} "
            f"--dfe Gamma_K17 -s 183 --dfe-bed-file {tmp_path / 'ex.bed'} "
            f"pop_0:5"
        ).split()
        capture_output(stdpopsim.cli.stdpopsim_main, cmd)
        ts = tskit.load(fname)
        assert ts.sequence_length == right - left
        self.verify_slim_sim(ts, num_samples=10)

    @pytest.mark.filterwarnings("ignore::stdpopsim.SLiMScalingFactorWarning")
    def test_chromosomal_segment_with_dfe_annotation(self):
        left = 37500000
        right = 38000000
        with tempfile.NamedTemporaryFile(mode="w") as f:
            cmd = (
                f"-q -e slim --slim-scaling-factor 20 --slim-path {slim_path} "
                f"HomSap -c chr22 --left {left} --right {right} -o {f.name} "
                f"--dfe Gamma_K17 -s 913 --dfe-annotation ensembl_havana_104_CDS "
                f"pop_0:5"
            ).split()
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)
            ts = tskit.load(f.name)
        assert ts.sequence_length == right - left
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
                extended_events = stdpopsim.ext.selective_sweep(
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
        extended_events = stdpopsim.ext.selective_sweep(
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
        extended_events = stdpopsim.ext.selective_sweep(
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
            stdpopsim.ext.selective_sweep(
                single_site_id="irrelevant",
                population="irrelevant",
                mutation_generation_ago=1000,
                selection_coeff=-0.1,
            )

    def test_sweep_with_bad_AF_conditions(self):
        for start_freq, end_freq in zip([-0.1, 0.1], [0.1, -0.1]):
            with pytest.raises(ValueError, match="of the sweep must be in"):
                stdpopsim.ext.selective_sweep(
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
            stdpopsim.ext.selective_sweep(
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
        extended_events = stdpopsim.ext.selective_sweep(
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
        extended_events = stdpopsim.ext.selective_sweep(
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
            extended_events += stdpopsim.ext.selective_sweep(
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
        extended_events = stdpopsim.ext.selective_sweep(
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
