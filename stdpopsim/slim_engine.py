"""
SLiM simulation engine.

This is a translation of the msprime API into SLiM's Eidos language, which
resembles R. The generated SLiM script is designed differently to the recipes
described in the SLiM reference manual. In our generated SLiM script, all the
demographic model parameters are defined in multi-dimensional arrays at the
top of the script, in the `initialize()` block. These arrays define the event
generations, and event blocks are subsequently constructed programmatically
using `sim.registerLateEvent()`, rather than writing out the blocks verbatim.
This design is intended to permit modification of demographic parameters in
the generated SLiM script, without needing to directly convert event times in
the past into forwards-time generations.

How backwards-time demographic events are mapped to forwards-time SLiM code:

 * `msprime.DemographyDebugger()` does much of the hard work by extracting
   epochs from the given model's `demographic_events`, and calculating a
   migration_matrix for each epoch from the `msprime.MigrationRateChange`
   events. The epoch boundaries defined here are indirectly translated into
   "late events" in SLiM.

 * `msprime.PopulationParametersChange` events are translated into SLiM as
   `pop.setSubpopulationSize()`. If `growth_rate` is not None, the population
   size is changed in every generation to match the specified rate.

 * `msprime.MassMigration` events with proportion=1 are population splits
   in forwards time. In SLiM, these are `sim.addSubpopSplit()`.

 * `msprime.MassMigration` events with proportion<1 indicate an admixture
   pulse at a single point in time. In SLiM, we call `pop.setMigrationRates()`
   in the relevant generation, and turn off migrations in the next generation.
   When multiple MassMigration events correspond to a single SLiM generation,
   the migration proportions multiply, following the msprime behaviour and
   event ordering.

 * The migration_matrix for each epoch describes continuous migrations that
   occur over long time periods. In SLiM, we call `pop.setMigrationRates()`.
"""

import os
import sys
import string
import tempfile
import subprocess
import functools
import itertools
import collections
import contextlib
import random

import stdpopsim
import numpy as np
import msprime
import pyslim

_slim_upper = """
initialize() {
    defineConstant("verbosity", $verbosity);

    // Scaling factor to speed up simulation.
    // See SLiM manual:
    // `5.5 Rescaling population sizes to improve simulation performance`.
    defineConstant("Q", $scaling_factor);

    defineConstant("generation_time", $generation_time);
    defineConstant("mutation_rate", Q * $mutation_rate);
    defineConstant("recombination_rate", (1-(1-2*$recombination_rate)^Q)/2);
    defineConstant("chromosome_length", $chromosome_length);
    defineConstant("trees_file", "$trees_file");
    defineConstant("check_coalescence", $check_coalescence);
"""

_slim_lower = """
    defineConstant("N", asInteger(_N/Q));

    initializeTreeSeq(checkCoalescence=check_coalescence);
    initializeMutationRate(mutation_rate);
    initializeMutationType("m1", 0.5, "f", 0);
    initializeGenomicElementType("g1", m1, 1.0);
    initializeGenomicElement(g1, 0, chromosome_length-1);
    initializeRecombinationRate(recombination_rate);
}

function (void)dbg(string$ s, [integer$ debug_level = 2]) {
    if (verbosity >= debug_level)
        catn(sim.generation + ": " + s);
}

// Return the number of generations that separate t0 and t1.
function (integer)gdiff(numeric$ t0, numeric t1) {
    return asInteger(round((t0-t1)/generation_time/Q));
}

// Output tree sequence file and end the simulation.
function (void)end(void) {
    sim.treeSeqOutput(trees_file);
    sim.simulationFinished();
}

// Create initial populations and migration rates.
1 {
    // Create initial populations.
    for (i in 0:(num_populations-1)) {
        if (N[0,i] > 0) {
            dbg("sim.addSubpop("+i+", "+N[0,i]+");");
            sim.addSubpop(i, N[0,i]);
        }
    }

    // Migration rates.
    i = 0;
    for (j in 0:(num_populations-1)) {
        for (k in 0:(num_populations-1)) {
            if (j==k | N[i,j] == 0 | N[i,k] == 0)
                next;

            m = migration_matrices[k,j,i];
            p = sim.subpopulations[j];
            dbg("p"+j+".setMigrationRates("+k+", "+m+");");
            p.setMigrationRates(k, m);
        }
    }
}

// Check if burn in has completed.
1 {
    if (!check_coalescence) {
        setup();
        return;
    }

    if (sim.treeSeqCoalesced()) {
        /*
         * All current generation individuals now have a common ancestor
         * born after the start of our simulation.  But the genealogy at this
         * first coalescence is not a good representative of an average
         * genealogy, as the TMRCA is biased low.  To understand why,
         * consider a haploid population with constant N=2.  At first
         * coalescence, the TMRCA is always one generation, yet clearly
         * genealogies with longer TMRCA's are possible for this population.
         *
         * To obtain a burn-in genealogy drawn from the full distribution
         * of possible genealogies, we may continue the simulation forwards
         * until this ergodic process converges.  Unfortunately, it's not
         * clear how many more generations are required for convergence,
         * and this is likely dependent on demography, and non-neutral
         * processes (if any).  10*N appears to be reasonable for a single
         * neutrally evolving population with constant N.
         */
        N_max = max(N[0,0:(num_populations-1)]);
        g = sim.generation + 10 * N_max;
        sim.registerEarlyEvent(NULL, "{setup();}", g, g);
    } else {
        if (sim.generation == 1)
            dbg("Waiting for burn in...");
        // Reschedule the current script block 10 generations hence.
        // XXX: find a less arbitrary generation interval.
        g = sim.generation + 10;
        sim.rescheduleScriptBlock(self, g, g);
    }
}

// Register events occurring at time _T[0] or more recently.
function (void)setup(void) {

    dbg("setup()");

    // Once burn in is complete, we know the starting generation (which
    // corresponds to T_0) and can thus calculate the generation
    // for each remaining event.
    G_start = sim.generation;
    T_0 = max(_T);
    G = G_start + gdiff(T_0, _T);
    G_end = max(G)+1;

    // Split events.
    if (length(subpopulation_splits) > 0 ) {
        for (i in 0:(ncol(subpopulation_splits)-1)) {
            g = G_start + gdiff(T_0, subpopulation_splits[0,i]);
            newpop = subpopulation_splits[1,i];
            size = asInteger(subpopulation_splits[2,i] / Q);
            oldpop = subpopulation_splits[3,i];
            sim.registerLateEvent(NULL,
                "{dbg(self.source); " +
                "sim.addSubpopSplit("+newpop+","+size+","+oldpop+");}",
                g, g);
        }
    }

    // Population size changes.
    if (num_epochs > 1) {
        for (i in 1:(num_epochs-1)) {
            g = G[i-1];
            for (j in 0:(num_populations-1)) {
                if (N[i,j] != N[i-1,j]) {
                    sim.registerLateEvent(NULL,
                        "{dbg(self.source); " +
                        "p"+j+".setSubpopulationSize("+N[i,j]+");}",
                        g, g);
                }

                if (growth_rates[i,j] != 0) {
                    growth_phase_start = g+1;
                    if (i == num_epochs-1)
                        growth_phase_end = G[i];
                    else {
                        // We already registered a size change at generation G[i].
                        growth_phase_end = G[i] - 1;
                    }

                    if (growth_phase_start >= growth_phase_end) {
                        // Some demographic models have duplicate epoch times,
                        // which should be ignored.
                        next;
                    }

                    N0 = N[i,j];
                    r = Q * growth_rates[i,j];

                    sim.registerLateEvent(NULL,
                        "{" +
                            "dbg(self.source); " +
                            "gx=sim.generation-"+g+"; " +
                            "size=asInteger(round("+N0+"*exp("+r+"*gx))); " +
                            "p"+j+".setSubpopulationSize(size);" +
                        "}",
                        growth_phase_start, growth_phase_end);
                }
            }
        }

        // Migration rates.
        for (i in 1:(num_epochs-1)) {
            for (j in 0:(num_populations-1)) {
                for (k in 0:(num_populations-1)) {
                    if (j==k | N[i,j] == 0 | N[i,k] == 0)
                        next;

                    m_last = migration_matrices[k,j,i-1];
                    m = migration_matrices[k,j,i];
                    if (m == m_last) {
                        // Do nothing if the migration rate hasn't changed.
                        next;
                    }
                    g = G[i-1];
                    sim.registerLateEvent(NULL,
                        "{dbg(self.source); " +
                        "p"+j+".setMigrationRates("+k+", "+m+");}",
                        g, g);
                }
            }
        }
    }

    // Admixture pulses.
    if (length(admixture_pulses) > 0 ) {
        for (i in 0:(ncol(admixture_pulses)-1)) {
            g = G_start + gdiff(T_0, admixture_pulses[0,i]);
            dest = admixture_pulses[1,i];
            src = admixture_pulses[2,i];
            rate = admixture_pulses[3,i];
            sim.registerLateEvent(NULL,
                "{dbg(self.source); " +
                "p"+dest+".setMigrationRates("+src+", "+rate+");}",
                g, g);
            sim.registerLateEvent(NULL,
                "{dbg(self.source); " +
                "p"+dest+".setMigrationRates("+src+", 0);}",
                g+1, g+1);
        }
    }

    // Sample individuals.
    for (i in 0:(ncol(sampling_episodes)-1)) {
        pop = sampling_episodes[0,i];
        n = sampling_episodes[1,i];
        g = G_start + gdiff(T_0, sampling_episodes[2,i]);
        sim.registerLateEvent(NULL,
            "{dbg(self.source); " +
            "inds=p"+pop+".sampleIndividuals("+n+"); " +
            "sim.treeSeqRememberIndividuals(inds);}",
            g, g);
    }

    sim.registerLateEvent(NULL, "{dbg(self.source); end();}", G_end, G_end);
}
"""


def slim_makescript(
        script_file, trees_file,
        demographic_model, contig, samples,
        scaling_factor, check_coalescence, verbosity):

    recombination_map = contig.recombination_map
    if len(recombination_map.get_positions()) > 2:
        raise NotImplementedError("recombination_map not supported")

    pop_names = [pc.metadata["id"] for pc in demographic_model.population_configurations]

    # Reassign event times according to integral SLiM generations.
    # This collapses the time deltas used in HomSap/AmericanAdmixture_4B11.
    for event in demographic_model.demographic_events:
        event.time = round(event.time / scaling_factor) * scaling_factor

    # The demography debugger constructs event epochs, which we use
    # to define the forwards-time events.
    dd = msprime.DemographyDebugger(
            population_configurations=demographic_model.population_configurations,
            migration_matrix=demographic_model.migration_matrix,
            demographic_events=demographic_model.demographic_events)

    epochs = sorted(dd.epochs, key=lambda e: e.start_time, reverse=True)
    T = [round(e.start_time*demographic_model.generation_time) for e in epochs]
    migration_matrices = [e.migration_matrix for e in epochs]

    N = np.empty(shape=(dd.num_populations, len(epochs)), dtype=int)
    growth_rates = np.empty(shape=(dd.num_populations, len(epochs)),
                            dtype=float)
    for j, epoch in enumerate(epochs):
        for i, pop in enumerate(epoch.populations):
            N[i, j] = int(pop.end_size)
            growth_rates[i, j] = pop.growth_rate

    admixture_pulses = []
    subpopulation_splits = []
    for i, epoch in enumerate(epochs):
        for de in epoch.demographic_events:
            if isinstance(de, msprime.MassMigration):

                if de.proportion < 1:
                    # Calculate remainder of population after previous
                    # MassMigration events in this epoch.
                    rem = 1 - np.sum([ap[3] for ap in admixture_pulses
                                     if ap[0] == i and ap[1] == de.source])
                    admixture_pulses.append((
                        i,
                        de.source,  # forwards-time dest
                        de.dest,    # forwards-time source
                        rem*de.proportion))
                    continue

                # Backwards: de.source is being merged into de.dest.
                # Forwards: de.source is being created, taking individuals
                #           from de.dest.
                #
                # If the proportion==1, we can use SLiM function:
                #       sim.addSubpopSplit(newpop, size, oldpop),
                # which we trigger by adding a row to subpopulation_splits.
                # This SLiM function creates newpop (=de.source), under the
                # assumption that it doesn't already exist.

                subpopulation_splits.append((
                    f"_T[{i}]",
                    de.source,
                    f"_N[{i+1},{de.source}]",
                    de.dest))

                # Zero out the population size for generations before this
                # epoch, to avoid simulating invididuals that contribute no
                # genealogy.
                N[de.source, 0:(i+1)] = 0
                growth_rates[de.source, 0:(i+1)] = 0

                # Ensure there are no migrations to or from de.source before
                # this epoch.
                for j in range(i+1):
                    for k in range(dd.num_populations):
                        migration_matrices[j][k][de.source] = 0
                        migration_matrices[j][de.source][k] = 0

    printsc = functools.partial(print, file=script_file)

    # Header
    printsc('/*')
    printsc(' * stdpopsim ' + stdpopsim.__version__)
    printsc(' *')
    printsc(' * Demographic model: ' + demographic_model.id)
    printsc(' * ' + "\n * ".join(
        [line.strip() for line in demographic_model.description.split('\n')]))
    for citation in demographic_model.citations:
        printsc(' * ' + str(citation))
    printsc(' */')

    printsc(string.Template(_slim_upper).substitute(
                scaling_factor=scaling_factor if scaling_factor is not None else 1,
                chromosome_length=int(recombination_map.get_length()),
                recombination_rate=recombination_map.mean_recombination_rate,
                mutation_rate=contig.mutation_rate,
                generation_time=demographic_model.generation_time,
                trees_file=trees_file,
                verbosity=verbosity,
                check_coalescence="T" if check_coalescence else "F",
                ))

    def matrix2str(matrix, row_comments=None, col_comment=None, indent=2,
                   fmt="", dim=(None, None)):
        """
        Return an Eidos representation of the matrix as a string.
        """
        if row_comments is not None:
            assert len(matrix) == len(row_comments)

        if len(matrix) == 0:
            return "c()"

        s = ["array(c(\n"]
        if col_comment is not None:
            s.append(indent*4*' ' + '// ' + col_comment + '\n')

        for i in range(len(matrix)):
            s.append(indent*4*" ")
            s.append('c({})'.format(", ".join(
                [format(x, fmt) for x in matrix[i]])))
            if i != len(matrix)-1:
                s.append(",")
            if row_comments is not None:
                s.append(" // " + row_comments[i])
            s.append("\n")

        s.append((indent-1)*4*" ")

        if dim[0] is None:
            dim = (len(matrix[0]), dim[1])
        if dim[1] is None:
            dim = (dim[0], len(matrix))
        s.append(f'), c({dim[0]}, {dim[1]}))')

        return "".join(s)

    printsc('    // Time of epoch boundaries, in years before present.')
    printsc('    // The first epoch spans from INF to _T[0].')
    printsc('    defineConstant("_T", c({}));'.format(", ".join(map(str, T))))
    printsc()

    # Population sizes.
    printsc('    // Population sizes in each epoch.')
    printsc('    _N = ' +
            matrix2str(
                N,
                row_comments=pop_names,
                col_comment="INF:_T[0], _T[0]:_T[1], etc.") +
            ';')
    printsc()

    printsc('    defineConstant("num_epochs", length(_T));')
    printsc('    defineConstant("num_populations", ncol(_N));')
    printsc()

    # Growth rates.
    printsc('    // Population growth rates for each epoch.')
    printsc('    defineConstant("growth_rates", ' +
            matrix2str(
                growth_rates,
                row_comments=pop_names,
                col_comment="INF:_T[0], _T[0]:_T[1], etc.",
                dim=("num_epochs", "num_populations")) +
            ');')
    printsc()

    printsc('    no_migration = rep(0, num_populations*num_populations);')
    printsc()

    # Migration rates.
    printsc('    // Migration rates for each epoch.')
    printsc('    // Migrations involving a population with size=0 are ignored.')
    printsc('    // XXX: document what the rows & cols correspond to.')
    printsc('    defineConstant("migration_matrices", array(c(')
    for i in range(len(migration_matrices)):
        epoch_str = f"INF:_T[{i}]" if i == 0 else f"_T[{i}]:_T[{i+1}]"
        printsc()
        printsc(2*4*' ' + '// ' + epoch_str)

        end = ",\n" if i != len(migration_matrices)-1 else "\n"
        if np.all(np.array(migration_matrices[i]) == 0):
            printsc(2*4*' ' + 'no_migration', end=end)
        else:
            printsc(2*4*' ' +
                    matrix2str(
                        migration_matrices[i],
                        indent=3,
                        fmt="g",
                        dim=("num_populations", "num_populations")),
                    end=end)
    printsc()
    printsc(4*' '+'), c(num_populations, num_populations, num_epochs)));')
    printsc()

    # Population splits.
    printsc('    // Population splits, one row for each event.')
    printsc('    defineConstant("subpopulation_splits", ' +
            matrix2str(
                subpopulation_splits,
                col_comment="time, newpop, size, oldpop") +
            ');')
    printsc()

    # Admixture pulses.
    # Output _T[...] variable rather than an index.
    admixture_pulses = [(f"_T[{ap[0]}]", *ap[1:]) for ap in admixture_pulses]
    printsc('    // Admixture pulses, one row for each pulse.')
    printsc('    defineConstant("admixture_pulses", ' +
            matrix2str(
                admixture_pulses,
                col_comment="time, dest, source, rate") +
            ');')
    printsc()

    # Sampling episodes.
    s_counts = collections.Counter([(s.population, s.time) for s in samples])
    sampling_episodes = []
    for i, ((pop, time), count) in enumerate(s_counts.items()):
        # XXX: SLiM can only sample individuals, which we assume are diploid.
        n_inds = (count+1) // 2
        sampling_episodes.append((pop, n_inds, int(time)))

    printsc('    // One row for each sampling episode.')
    printsc('    defineConstant("sampling_episodes", ' +
            matrix2str(
                sampling_episodes,
                col_comment='pop, n_inds, time') +
            ');')

    printsc(_slim_lower)


def simplify_remembered(ts):
    """
    Remove all samples except those individuals that were explicity
    sampled in SLiM with sim.treeSeqRememberIndividuals().
    """
    nodes = itertools.chain.from_iterable(
                i.nodes for i in ts.individuals()
                if i.flags & pyslim.INDIVIDUAL_REMEMBERED)
    return ts.simplify(samples=list(nodes))


class _SLiMEngine(stdpopsim.Engine):
    id = "slim"  #:
    description = "SLiM forward-time Wright-Fisher simulator"  #:
    citations = [
            stdpopsim.Citation(
                doi="https://doi.org/10.1111/1755-0998.12968",
                year=2019,
                author="Haller et al.",
                reasons={stdpopsim.CiteReason.ENGINE}),
            ]

    def slim_path(self):
        return os.environ.get("SLIM", "slim")

    def get_version(self):
        s = subprocess.check_output([self.slim_path(), "-v"])
        return s.split()[2].decode("ascii").rstrip(",")

    def simulate(
            self, demographic_model=None, contig=None, samples=None, seed=None,
            verbosity=0, slim_path=None, slim_script=False, slim_scaling_factor=10,
            slim_no_recapitation=False, slim_no_burnin=False, **kwargs):
        """
        Simulate the demographic model using SLiM.
        See :meth:`.Engine.simulate()` for definitions of the
        ``demographic_model``, ``contig``, and ``samples`` parameters.

        :param seed: The seed for the random number generator.
        :type seed: int
        :param slim_path: The full path to the slim executable, or the name of
            a command in the current PATH.
        :type slim_path: str
        :param slim_script: If true, the simulation will not be executed.
            Instead the generated SLiM script will be printed to stdout.
        :type slim_script: bool
        :param slim_scaling_factor: Rescale model parameters by the given value,
            to speed up simulation. Population sizes and generation times are
            divided by this factor, whereas the mutation rate, recombination
            rate, and growth rates are multiplied by the factor.
            See SLiM manual: `5.5 Rescaling population sizes to improve
            simulation performance.`
        :type slim_scaling_factor: float
        :param slim_no_recapitation: Do an explicit burn in, and add
            mutations, within the SLiM simulation. This may be much slower than
            the defaults (recapitation and neutral mutation overlay with
            msprime). The burn in behaviour is to wait until all individuals in
            the ancestral populations have a common ancestor within their
            respective population, and then wait another 10*N generations.
        :type slim_no_recapitation: bool
        :param slim_no_burnin: Do not perform a burn in at the start of the
            simulation.  This option is only relevant when
            ``slim_no_recapitation=True``.
        :type slim_no_burnin: bool
        """

        run_slim = not slim_script
        do_recap = not slim_no_recapitation
        check_coalescence = slim_no_recapitation and not slim_no_burnin

        if slim_path is None:
            slim_path = self.slim_path()

        if do_recap:
            mutation_rate = contig.mutation_rate
            # Ensure no mutations are introduced by SLiM.
            contig = stdpopsim.Contig(
                    recombination_map=contig.recombination_map,
                    mutation_rate=0,
                    genetic_map=contig.genetic_map)

        slim_cmd = [slim_path]
        if seed is not None:
            slim_cmd.extend(["-s", f"{seed}"])

        mktemp = functools.partial(tempfile.NamedTemporaryFile, mode="w")

        @contextlib.contextmanager
        def script_file_f():
            f = mktemp(suffix=".slim") if not slim_script else sys.stdout
            yield f
            # Don't close sys.stdout.
            if not slim_script:
                f.close()

        with script_file_f() as script_file, mktemp(suffix=".ts") as ts_file:

            slim_makescript(
                    script_file, ts_file.name,
                    demographic_model, contig, samples,
                    slim_scaling_factor, check_coalescence, verbosity)

            script_file.flush()

            if not run_slim:
                return None

            slim_cmd.append(script_file.name)
            stdout = subprocess.DEVNULL if verbosity == 0 else None
            subprocess.check_call(slim_cmd, stdout=stdout)

            ts = pyslim.load(ts_file.name)

        if do_recap:
            rng = random.Random(seed)
            s1, s2 = rng.randrange(1, 2**32), rng.randrange(1, 2**32)

            # Recapitation.
            r_map = contig.recombination_map
            pop_configs = demographic_model.population_configurations
            ts = ts.recapitate(
                    recombination_rate=r_map.mean_recombination_rate,
                    population_configurations=pop_configs,
                    random_seed=s1)

        ts = simplify_remembered(ts)

        if do_recap:
            # Add neutral mutations.
            ts = pyslim.SlimTreeSequence(msprime.mutate(
                ts, rate=mutation_rate, keep=True, random_seed=s2))

        return ts

    def add_arguments(self, parser):
        def slim_exec(path):
            # Hack to set the SLIM environment variable at parse time,
            # before get_version() can be called.
            os.environ["SLIM"] = path
            return path
        parser.add_argument(
                "--slim-path", metavar="PATH", type=slim_exec, default=None,
                help="Full path to `slim' executable.")
        parser.add_argument(
                "--slim-script", action="store_true", default=False,
                help="Write script to stdout and exit without running SLiM.")
        parser.add_argument(
                "--slim-scaling-factor", metavar="Q", default=10, type=float,
                help="Rescale model parameters by Q to speed up simulation. "
                     "See SLiM manual: `5.5 Rescaling population sizes to "
                     "improve simulation performance`. "
                     "[default=%(default)s].")
        parser.add_argument(
                "--slim-no-recapitation", action="store_true", default=False,
                help="Explicitly wait for coalescence, and add "
                     "mutations, within the SLiM simulation. This may be much "
                     "slower than the defaults (recapitation and neutral mutation "
                     "overlay with msprime).")
        parser.add_argument(
                "--slim-no-burnin", action="store_true", default=False,
                help="Don't wait for coalescence in SLiM before proceeding. "
                     "This option is only relevant in combination with "
                     "--slim-no-recapitation.")


stdpopsim.register_engine(_SLiMEngine())
