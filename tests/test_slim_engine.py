"""
Tests for SLiM simulation engine.
"""
import os
import re
import io
import sys
import itertools
import unittest
import tempfile
import math
from unittest import mock
import numpy as np

import tskit
import pyslim
import msprime

import stdpopsim
import stdpopsim.cli
from .test_cli import capture_output

IS_WINDOWS = sys.platform.startswith("win")


def slim_simulate_no_recap(seed=1234, **kwargs):
    """
    Return the tree sequence produced by SLiM, without recapitation, etc.
    """
    kwargs.update(slim_script=True)
    engine = stdpopsim.get_engine("slim")
    out, _ = capture_output(engine.simulate, **kwargs)

    # Find the name of the temporary trees_file in the script.
    match = re.search(r'"trees_file",\s*"([^"]*)"', out)
    assert match is not None
    tmp_trees_file = match.group(1)

    with tempfile.TemporaryDirectory() as tmpdir:
        script_file = os.path.join(tmpdir, "script.slim")
        trees_file = os.path.join(tmpdir, "out.trees")
        # Write out the script with a new location for the trees_file.
        out = out.replace(tmp_trees_file, trees_file)
        with open(script_file, "w") as f:
            f.write(out)
        engine._run_slim(script_file, seed=seed)
        ts = pyslim.load(trees_file)
    return ts


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
                    demographic_model=model,
                    contig=contig,
                    samples=samples,
                    slim_scaling_factor=scaling_factor,
                    dry_run=True,
                )

        for burn_in in (-1, -1e-6):
            with self.assertRaises(ValueError):
                engine.simulate(
                    demographic_model=model,
                    contig=contig,
                    samples=samples,
                    slim_burn_in=burn_in,
                    dry_run=True,
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
        self.assertTrue("sim.registerLateEvent" in out)

        model = species.get_demographic_model("AncientEurasia_9K19")
        samples = model.get_samples(10, 20, 30, 40, 50, 60, 70)
        out, _ = capture_output(
            engine.simulate,
            demographic_model=model,
            contig=contig,
            samples=samples,
            slim_script=True,
        )
        self.assertTrue("sim.registerLateEvent" in out)

        model = species.get_demographic_model("AmericanAdmixture_4B11")
        samples = model.get_samples(10, 10, 10)
        out, _ = capture_output(
            engine.simulate,
            demographic_model=model,
            contig=contig,
            samples=samples,
            slim_script=True,
        )
        self.assertTrue("sim.registerLateEvent" in out)

    def test_no_recombination_map(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr1", genetic_map="HapMapII_GRCh37")
        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        contig.recombination_map = None
        with self.assertRaises(ValueError):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=model.get_samples(10),
                dry_run=True,
            )

    def test_recombination_map(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr1", genetic_map="HapMapII_GRCh37")
        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        samples = model.get_samples(10)
        out, _ = capture_output(
            engine.simulate,
            demographic_model=model,
            contig=contig,
            samples=samples,
            dry_run=True,
        )

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
        self.assertEqual(ts.num_samples, 10)
        self.assertTrue(all(tree.num_roots == 1 for tree in ts.trees()))

    def test_recap_and_rescale(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        model = species.get_demographic_model("OutOfAfrica_3G09")
        samples = model.get_samples(10, 10, 10)
        for proportion, seed in zip((0, 1), (1234, 2345)):
            contig = species.get_contig("chr22", length_multiplier=0.001)
            if proportion:
                contig.fully_neutral(slim_mutations=True)
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
            ts2_headless = slim_simulate_no_recap(
                demographic_model=model,
                contig=contig,
                samples=samples,
                extended_events=extended_events,
                slim_scaling_factor=10,
                slim_burn_in=0,
                seed=seed,
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

            self.assertEqual(tables1.nodes, tables2.nodes)
            self.assertEqual(tables1.edges, tables2.edges)
            self.assertEqual(tables1.mutations, tables2.mutations)


@unittest.skipIf(IS_WINDOWS, "SLiM not available on windows")
class TestCLI(unittest.TestCase):
    def docmd(self, _cmd):
        cmd = (
            "-q -e slim --slim-scaling-factor 20 --slim-burn-in 0 "
            f"{_cmd} -l 0.001 -c chr1 -s 1234 10"
        ).split()
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
        slim_path = os.environ.get("SLIM", "slim")
        with tempfile.NamedTemporaryFile(mode="w") as f:
            self.docmd(f"--slim-path {slim_path} HomSap -o {f.name}")
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
            cmd = (
                "-q -e slim --slim-scaling-factor 20 --slim-burn-in 0 "
                f"HomSap -o {f.name} -l 0.001 -c chr1 -s 1234 "
                "-d OutOfAfrica_3G09 0 0 8"
            ).split()
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
        slim_path = os.environ.get("SLIM", "slim")
        self.assertTrue(slim_path in mocked_popen.call_args[0][0])
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

    # this is an expected failure, because no warning should be emitted
    @unittest.expectedFailure
    def test_odd_sample_warning_for_even_samples(self):
        cmd = "-q -e slim --slim-script HomSap -d OutOfAfrica_2T12 4 6".split()
        with self.assertWarns(stdpopsim.SLiMOddSampleWarning):
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)

    def test_odd_sample_warning(self):
        cmd = "-q -e slim --slim-script HomSap -d OutOfAfrica_2T12 -L 100 4 5".split()
        with self.assertWarns(stdpopsim.SLiMOddSampleWarning):
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)

        cmd = "-q -e slim --slim-script HomSap -d OutOfAfrica_2T12 -L 100 3 5".split()
        with self.assertWarns(stdpopsim.SLiMOddSampleWarning):
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)

    def triplet(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr22", length_multiplier=0.001)
        return engine, species, contig

    def test_bad_population_size_addSubPop(self):
        engine, species, contig = self.triplet()
        model = stdpopsim.PiecewiseConstantSize(100)
        samples = model.get_samples(2)

        with self.assertWarns(stdpopsim.UnspecifiedSLiMWarning):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                slim_scaling_factor=10,
                dry_run=True,
            )

    def test_no_populations_in_generation_1(self):
        engine, species, contig = self.triplet()
        model = stdpopsim.PiecewiseConstantSize(100)
        samples = model.get_samples(2)

        with self.assertRaises(stdpopsim.SLiMException):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                slim_scaling_factor=200,
                dry_run=True,
            )

    def test_bad_population_size_addSubpopSplit(self):
        engine, species, contig = self.triplet()
        model = stdpopsim.IsolationWithMigration(
            NA=1000, N1=100, N2=1000, T=1000, M12=0, M21=0
        )
        samples = model.get_samples(2)

        with self.assertWarns(stdpopsim.UnspecifiedSLiMWarning):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                slim_scaling_factor=10,
                dry_run=True,
            )

        with self.assertRaises(stdpopsim.SLiMException):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                slim_scaling_factor=200,
                dry_run=True,
            )

    def test_bad_population_size_setSubpopulationSize(self):
        engine, species, contig = self.triplet()
        model = stdpopsim.PiecewiseConstantSize(100, (1000, 1000))
        samples = model.get_samples(2)

        with self.assertWarns(stdpopsim.UnspecifiedSLiMWarning):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                slim_scaling_factor=10,
                dry_run=True,
            )

        with self.assertRaises(stdpopsim.SLiMException):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                slim_scaling_factor=200,
                dry_run=True,
            )

    def test_sample_size_too_big(self):
        engine, species, contig = self.triplet()
        model = stdpopsim.PiecewiseConstantSize(1000)
        samples = model.get_samples(300)

        with self.assertRaises(stdpopsim.SLiMException):
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
        return stdpopsim.DemographicModel(
            id="exp_decline",
            description="exp_decline",
            long_description="exp_decline",
            populations=[stdpopsim.models._pop0],
            generation_time=1,
            population_configurations=[
                msprime.PopulationConfiguration(
                    initial_size=N0,
                    growth_rate=r,
                    metadata=stdpopsim.models._pop0.asdict(),
                )
            ],
            demographic_events=[
                msprime.PopulationParametersChange(
                    time=T, initial_size=N1, growth_rate=0, population_id=0
                ),
            ],
        )

    def test_bad_population_size_exp_decline(self):
        engine, species, contig = self.triplet()
        model = self.exp_decline()
        samples = model.get_samples(2)

        with self.assertWarns(stdpopsim.UnspecifiedSLiMWarning):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                slim_scaling_factor=10,
                dry_run=True,
            )

        with self.assertRaises(stdpopsim.SLiMException):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                slim_scaling_factor=200,
                dry_run=False,
            )

    def test_sample_size_too_big_exp_decline(self):
        engine, species, contig = self.triplet()
        model = self.exp_decline()
        samples = model.get_samples(30)

        with self.assertRaises(stdpopsim.SLiMException):
            engine.simulate(
                demographic_model=model,
                contig=contig,
                samples=samples,
                slim_scaling_factor=10,
                dry_run=True,
            )

    # this is an expected failure, because no warning should be emitted
    @unittest.expectedFailure
    def test_warning_when_not_scaling(self):
        with self.assertWarns(stdpopsim.SLiMScalingFactorWarning):
            for cmd in [
                "HomSap 100 -D",
                "-e slim HomSap 100 -D",
                "-e slim --slim-scaling-factor 1 HomSap 100 -D",
                "-e slim --slim-scaling-factor 1.0 HomSap 100 -D",
            ]:
                capture_output(stdpopsim.cli.stdpopsim_main, cmd.split())

    def test_warning_when_scaling(self):
        for cmd in [
            "-e slim --slim-scaling-factor 2 HomSap 100 -D -L 1",
            "-e slim --slim-scaling-factor 1000 EscCol 100 -D -L 1",
        ]:
            with self.assertWarns(stdpopsim.SLiMScalingFactorWarning):
                capture_output(stdpopsim.cli.stdpopsim_main, cmd.split())


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


def get_test_contig(spp="HomSap", chrom="chr22", length_multiplier=0.001):
    species = stdpopsim.get_species(spp)
    contig = species.get_contig(chrom, length_multiplier=length_multiplier)
    contig.fully_neutral(convert_to_substitution=False)
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
    mut_id = 0

    def allele_frequency(self, ts):
        """
        Get the allele frequency of the drawn mutation.
        """
        # surely there's a simpler way!
        assert ts.num_mutations == 1
        alive = list(
            itertools.chain.from_iterable(
                ts.individual(i).nodes for i in ts.individuals_alive_at(0)
            )
        )
        mut = next(ts.mutations())
        tree = ts.at(ts.site(mut.site).position)
        have_mut = [u for u in alive if tree.is_descendant(u, mut.node)]
        af = len(have_mut) / len(alive)
        return af


@unittest.skipIf(IS_WINDOWS, "SLiM not available on windows")
class TestGenomicElementTypes(unittest.TestCase, PiecewiseConstantSizeMixin):
    def test_single_genomic_element_type_in_script(self):
        contig = get_test_contig()
        engine = stdpopsim.get_engine("slim")
        out, _ = capture_output(
            engine.simulate,
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            slim_script=True,
        )
        self.assertEqual(out.count("initializeGenomicElementType"), 1)
        self.assertEqual(out.count("initializeGenomicElement("), 1)

    def test_multiple_genomic_element_types_in_script(self):
        contig = get_test_contig()
        engine = stdpopsim.get_engine("slim")
        contig.clear_genomic_mutation_types()
        contig.add_genomic_element_type(
            intervals=np.array([[0, 100]]),
            mutation_types=[stdpopsim.ext.MutationType()],
            proportions=[0],
        )
        contig.add_genomic_element_type(
            intervals=np.array([[100, 200]]),
            mutation_types=[stdpopsim.ext.MutationType()],
            proportions=[0],
        )
        contig.add_genomic_element_type(
            intervals=np.array([[200, 300]]),
            mutation_types=[stdpopsim.ext.MutationType()],
            proportions=[0],
        )
        out, _ = capture_output(
            engine.simulate,
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            slim_script=True,
        )
        self.assertEqual(out.count("initializeGenomicElementType"), 3)
        self.assertEqual(out.count("initializeGenomicElement("), 3)


@unittest.skipIf(IS_WINDOWS, "SLiM not available on windows")
class TestMutationTypes(unittest.TestCase, PiecewiseConstantSizeMixin):
    def test_single_mutation_type_in_script(self):
        contig = get_test_contig()
        engine = stdpopsim.get_engine("slim")
        out, _ = capture_output(
            engine.simulate,
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            slim_script=True,
        )
        self.assertEqual(out.count("initializeMutationType"), 1)

        out, _ = capture_output(
            engine.simulate,
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            slim_script=True,
        )
        self.assertEqual(out.count("initializeMutationType"), 1)

    def test_multiple_mutation_types_in_script(self):
        contig = get_test_contig()
        engine = stdpopsim.get_engine("slim")
        contig.mutation_types = [
            stdpopsim.ext.MutationType(),
            stdpopsim.ext.MutationType(),
        ]
        out, _ = capture_output(
            engine.simulate,
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            slim_script=True,
        )
        self.assertEqual(out.count("initializeMutationType"), 2)

        positive = stdpopsim.ext.MutationType(convert_to_substitution=False)
        contig.mutation_types = [stdpopsim.ext.MutationType() for i in range(10)] + [
            positive
        ]
        contig.genomic_element_types[0].proportions = [1 / 11 for i in range(11)]
        contig.genomic_element_types[0].mutation_type_ids = [i for i in range(11)]
        out, _ = capture_output(
            engine.simulate,
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            slim_script=True,
        )
        self.assertEqual(out.count("initializeMutationType"), 11)

    def test_unweighted_mutations_are_not_simulated_by_slim(self):
        contig = get_test_contig()
        contig.mutation_types = [
            stdpopsim.ext.MutationType(convert_to_substitution=True),
            stdpopsim.ext.MutationType(convert_to_substitution=False),
        ]
        contig.genomic_element_types[0].mutation_type_ids = [0, 1]
        contig.genomic_element_types[0].proportions = [0, 0]
        ts = slim_simulate_no_recap(
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            slim_scaling_factor=10,
            slim_burn_in=0.1,
        )
        self.assertEqual(ts.num_sites, 0)

        contig.mutation_types = [
            stdpopsim.ext.MutationType(),
            stdpopsim.ext.MutationType(
                distribution_type="g", distribution_args=[-0.01, 0.2]
            ),
        ]
        contig.genomic_element_types[0].proportions = [0, 0]
        ts = slim_simulate_no_recap(
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            slim_scaling_factor=10,
            slim_burn_in=0.1,
        )
        self.assertEqual(ts.num_sites, 0)

    def test_weighted_mutations_are_simulated_by_slim(self):
        contig = get_test_contig()
        contig.mutation_types = [
            stdpopsim.ext.MutationType(convert_to_substitution=True)
        ]
        contig.genomic_element_types[0].proportions = [1]
        ts = slim_simulate_no_recap(
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            slim_scaling_factor=10,
            slim_burn_in=0.1,
        )
        self.assertTrue(ts.num_sites > 0)

        contig.mutation_types = [
            stdpopsim.ext.MutationType(convert_to_substitution=False)
        ]
        contig.genomic_element_types[0].proportions = [1]
        ts = slim_simulate_no_recap(
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            slim_scaling_factor=10,
            slim_burn_in=0.1,
        )
        self.assertTrue(ts.num_sites > 0)

        contig.mutation_types = [stdpopsim.ext.MutationType() for i in range(10)]
        contig.genomic_element_types[0].proportions = [1 / 10 for i in range(10)]
        contig.genomic_element_types[0].mutation_type_ids = list(range(10))
        ts = slim_simulate_no_recap(
            demographic_model=self.model,
            contig=contig,
            samples=self.samples,
            slim_scaling_factor=10,
            slim_burn_in=0.1,
        )
        self.assertTrue(ts.num_sites > 0)

    def test_dominance_coeff(self):
        for dominance_coeff in (0, 0.5, 1, 50):
            stdpopsim.ext.MutationType(dominance_coeff=dominance_coeff)

    def test_bad_dominance_coeff(self):
        for dominance_coeff in (-1,):
            with self.assertRaises(ValueError):
                stdpopsim.ext.MutationType(dominance_coeff=dominance_coeff)

    def test_bad_distribution_type(self):
        for distribution_type in (1, {}, None, "~", "!", "F"):
            with self.assertRaises(ValueError):
                stdpopsim.ext.MutationType(distribution_type=distribution_type)

    def test_distribution_type_f(self):
        for distribution_args in ([-0.1], [0], [0.1], [50]):
            stdpopsim.ext.MutationType(
                distribution_type="f", distribution_args=distribution_args
            )

    def test_bad_distribution_args_f(self):
        for distribution_args in ([0.1, 0.2], []):
            with self.assertRaises(ValueError):
                stdpopsim.ext.MutationType(
                    distribution_type="f", distribution_args=distribution_args
                )

    def test_distribution_type_g(self):
        for distribution_args in ([-0.1, 0.1], [0.1, 0.1], [50, 50]):
            stdpopsim.ext.MutationType(
                distribution_type="g", distribution_args=distribution_args
            )

    def test_bad_distribution_args_g(self):
        for distribution_args in ([], [0.1, 0], [0.1, -0.1], [0.1, 0.4, 0.5]):
            with self.assertRaises(ValueError):
                stdpopsim.ext.MutationType(
                    distribution_type="g", distribution_args=distribution_args
                )


@unittest.skipIf(IS_WINDOWS, "SLiM not available on windows")
class TestDrawMutation(unittest.TestCase, PiecewiseConstantSizeMixin):
    def test_draw_mutation_save(self):
        extended_events = [
            stdpopsim.ext.DrawMutation(
                time=self.T_mut,
                mutation_type_id=self.mut_id,
                population_id=0,
                coordinate=100,
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
                mutation_type_id=self.mut_id,
                population_id=0,
                coordinate=100,
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

    def test_invalid_mutation_type_id(self):
        engine = stdpopsim.get_engine("slim")
        for mut_type_id in [-1, 10]:
            extended_events = [
                stdpopsim.ext.DrawMutation(
                    time=self.T_mut,
                    mutation_type_id=mut_type_id,
                    population_id=0,
                    coordinate=100,
                ),
            ]
            with self.assertRaises(ValueError):
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
                mutation_type_id=self.mut_id,
                population_id=0,
                coordinate=100,
            ),
        ]
        contig = self.contig
        contig.mutation_types = []
        engine = stdpopsim.get_engine("slim")
        with self.assertRaises(ValueError):
            engine.simulate(
                demographic_model=self.model,
                contig=contig,
                samples=self.samples,
                extended_events=extended_events,
                dry_run=True,
            )

    def test_bad_time(self):
        for time in (-1,):
            with self.assertRaises(ValueError):
                stdpopsim.ext.DrawMutation(
                    time=time, mutation_type_id=0, population_id=0, coordinate=0
                )
        for time in (0, -1):
            with self.assertRaises(ValueError):
                stdpopsim.ext.DrawMutation(
                    time=stdpopsim.ext.GenerationAfter(time),
                    mutation_type_id=0,
                    population_id=0,
                    coordinate=0,
                )


@unittest.skipIf(IS_WINDOWS, "SLiM not available on windows")
class TestAlleleFrequencyConditioning(unittest.TestCase, PiecewiseConstantSizeMixin):
    def test_save_point_creation(self):
        extended_events = [
            stdpopsim.ext.DrawMutation(
                time=self.T_mut,
                mutation_type_id=self.mut_id,
                population_id=0,
                coordinate=100,
                save=True,
            ),
            stdpopsim.ext.ConditionOnAlleleFrequency(
                start_time=stdpopsim.ext.GenerationAfter(self.T_mut),
                end_time=0,
                mutation_type_id=self.mut_id,
                population_id=0,
                op=">",
                allele_frequency=0,
                save=True,
            ),
            stdpopsim.ext.ConditionOnAlleleFrequency(
                start_time=self.T_mut // 2,
                end_time=self.T_mut // 2,
                mutation_type_id=self.mut_id,
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

    def test_drawn_mutation_not_lost(self):
        extended_events = [
            stdpopsim.ext.DrawMutation(
                time=self.T_mut,
                mutation_type_id=self.mut_id,
                population_id=0,
                coordinate=100,
                save=True,
            ),
            stdpopsim.ext.ConditionOnAlleleFrequency(
                start_time=0,
                end_time=0,
                mutation_type_id=self.mut_id,
                population_id=0,
                op=">",
                allele_frequency=0,
            ),
        ]
        ts = slim_simulate_no_recap(
            demographic_model=self.model,
            contig=self.contig,
            samples=self.samples,
            extended_events=extended_events,
            slim_scaling_factor=10,
            slim_burn_in=0.1,
        )
        self.assertEqual(ts.num_mutations, 1)

    def test_drawn_mutation_is_lost(self):
        extended_events = [
            stdpopsim.ext.DrawMutation(
                time=self.T_mut,
                mutation_type_id=self.mut_id,
                population_id=0,
                coordinate=100,
                save=True,
            ),
            stdpopsim.ext.ConditionOnAlleleFrequency(
                start_time=0,
                end_time=0,
                mutation_type_id=self.mut_id,
                population_id=0,
                op="<=",
                allele_frequency=0,
            ),
        ]
        ts = slim_simulate_no_recap(
            demographic_model=self.model,
            contig=self.contig,
            samples=self.samples,
            extended_events=extended_events,
            slim_scaling_factor=10,
            slim_burn_in=0.1,
        )
        self.assertEqual(ts.num_mutations, 0)

    def test_drawn_mutation_meets_AF_threshold(self):
        for af_threshold, seed in zip((0.01, 0.1, 0.2), (1, 2, 3)):
            extended_events = [
                stdpopsim.ext.DrawMutation(
                    time=self.T_mut,
                    mutation_type_id=self.mut_id,
                    population_id=0,
                    coordinate=100,
                    save=True,
                ),
                # Condition on desired AF at end of simulation.
                stdpopsim.ext.ConditionOnAlleleFrequency(
                    start_time=0,
                    end_time=0,
                    mutation_type_id=self.mut_id,
                    population_id=0,
                    op=">=",
                    allele_frequency=af_threshold,
                ),
            ]
            ts = slim_simulate_no_recap(
                demographic_model=self.model,
                contig=self.contig,
                samples=self.samples,
                extended_events=extended_events,
                slim_scaling_factor=10,
                slim_burn_in=0.1,
                seed=seed,
            )
            self.assertEqual(ts.num_mutations, 1)
            self.assertTrue(self.allele_frequency(ts) >= af_threshold)

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
            with self.assertRaises(ValueError):
                stdpopsim.ext.ConditionOnAlleleFrequency(
                    start_time=0,
                    end_time=0,
                    mutation_type_id=self.mut_id,
                    population_id=0,
                    op=op,
                    allele_frequency=af,
                )

    def test_bad_times(self):
        for start_time, end_time in [(-1, 0), (0, -1), (1, 100)]:
            with self.assertRaises(ValueError):
                stdpopsim.ext.ConditionOnAlleleFrequency(
                    start_time=start_time,
                    end_time=end_time,
                    mutation_type_id=0,
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
            with self.assertRaises(ValueError):
                extended_events = [
                    stdpopsim.ext.DrawMutation(
                        time=self.T_mut,
                        mutation_type_id=self.mut_id,
                        population_id=0,
                        coordinate=100,
                        save=True,
                    ),
                    stdpopsim.ext.ConditionOnAlleleFrequency(
                        start_time=stdpopsim.ext.GenerationAfter(start_time),
                        end_time=end_time,
                        mutation_type_id=self.mut_id,
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
            self.assertTrue(0 <= id < len(op_types))
        for op in ("==", "=", "!=", {}, ""):
            with self.assertRaises(ValueError):
                id = stdpopsim.ext.ConditionOnAlleleFrequency.op_id(op)

    def test_conditioning_without_save(self):
        extended_events = [
            stdpopsim.ext.DrawMutation(
                time=self.T_mut,
                mutation_type_id=self.mut_id,
                population_id=0,
                coordinate=100,
            ),
            stdpopsim.ext.ConditionOnAlleleFrequency(
                start_time=stdpopsim.ext.GenerationAfter(self.T_mut),
                end_time=0,
                mutation_type_id=self.mut_id,
                population_id=0,
                op=">=",
                allele_frequency=1,
            ),
        ]
        engine = stdpopsim.get_engine("slim")
        with self.assertRaises(stdpopsim.SLiMException):
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
                mutation_type_id=self.mut_id,
                population_id=0,
                op=">",
                allele_frequency=0,
                save=True,
            ),
        ]
        engine = stdpopsim.get_engine("slim")
        with self.assertRaises(ValueError):
            engine.simulate(
                demographic_model=self.model,
                contig=self.contig,
                samples=self.samples,
                extended_events=extended_events,
                dry_run=True,
            )


@unittest.skipIf(IS_WINDOWS, "SLiM not available on windows")
class TestChangeMutationFitness(unittest.TestCase, PiecewiseConstantSizeMixin):
    # Testing stdpopsim.ext.ChangeMutationFitness is challenging, because
    # the side-effects are not deterministic. But if we condition on fixation
    # of a drawn mutation, such a simulation will be very slow without strong
    # positive selection (because we effectively do rejection sampling on the
    # simulation until we get one that meets the allele frequency condition).
    # So if this test takes more than a few seconds to run, that's a good
    # indication that selection is not acting.
    def test_positive_mutation_meets_AF_threshold(self):
        for af_threshold, seed in zip((0.5, 1), (1, 2)):
            extended_events = [
                stdpopsim.ext.DrawMutation(
                    time=self.T_mut,
                    mutation_type_id=self.mut_id,
                    population_id=0,
                    coordinate=100,
                    save=True,
                ),
                stdpopsim.ext.ChangeMutationFitness(
                    start_time=stdpopsim.ext.GenerationAfter(self.T_mut),
                    end_time=0,
                    mutation_type_id=self.mut_id,
                    population_id=0,
                    selection_coeff=0.1,
                    dominance_coeff=0.5,
                ),
                # Condition on AF > 0, to restore() immediately if the
                # allele is lost.
                stdpopsim.ext.ConditionOnAlleleFrequency(
                    start_time=0,
                    end_time=0,
                    mutation_type_id=self.mut_id,
                    population_id=0,
                    op=">",
                    allele_frequency=0,
                ),
                # Condition on desired AF at end of simulation.
                stdpopsim.ext.ConditionOnAlleleFrequency(
                    start_time=0,
                    end_time=0,
                    mutation_type_id=self.mut_id,
                    population_id=0,
                    op=">=",
                    allele_frequency=af_threshold,
                ),
            ]
            ts = slim_simulate_no_recap(
                demographic_model=self.model,
                contig=self.contig,
                samples=self.samples,
                extended_events=extended_events,
                slim_scaling_factor=10,
                slim_burn_in=0.1,
                seed=seed,
            )
            self.assertEqual(ts.num_mutations, 1)
            self.assertTrue(self.allele_frequency(ts) >= af_threshold)

    def test_no_drawn_mutation(self):
        extended_events = [
            stdpopsim.ext.ChangeMutationFitness(
                start_time=stdpopsim.ext.GenerationAfter(self.T_mut),
                end_time=0,
                mutation_type_id=self.mut_id,
                population_id=0,
                selection_coeff=0.1,
                dominance_coeff=0.5,
            ),
        ]
        engine = stdpopsim.get_engine("slim")
        with self.assertRaises(ValueError):
            engine.simulate(
                demographic_model=self.model,
                contig=self.contig,
                samples=self.samples,
                extended_events=extended_events,
                dry_run=True,
            )

    def test_bad_times(self):
        for start_time, end_time in [(-1, 0), (0, -1), (1, 100)]:
            with self.assertRaises(ValueError):
                stdpopsim.ext.ChangeMutationFitness(
                    start_time=start_time,
                    end_time=end_time,
                    mutation_type_id=0,
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
            with self.assertRaises(ValueError):
                extended_events = [
                    stdpopsim.ext.DrawMutation(
                        time=self.T_mut,
                        mutation_type_id=self.mut_id,
                        population_id=0,
                        coordinate=100,
                        save=True,
                    ),
                    stdpopsim.ext.ChangeMutationFitness(
                        start_time=stdpopsim.ext.GenerationAfter(start_time),
                        end_time=end_time,
                        mutation_type_id=self.mut_id,
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


@unittest.skipIf(IS_WINDOWS, "SLiM not available on windows")
class TestExtendedEvents(unittest.TestCase, PiecewiseConstantSizeMixin):
    def test_bad_extended_events(self):
        engine = stdpopsim.get_engine("slim")
        for bad_ee in [
            msprime.PopulationParametersChange(time=0, initial_size=100),
            None,
            {},
            "",
        ]:
            with self.assertRaises(ValueError):
                engine.simulate(
                    demographic_model=self.model,
                    contig=self.contig,
                    samples=self.samples,
                    extended_events=[bad_ee],
                    dry_run=True,
                )
