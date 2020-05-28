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
import textwrap
import logging
import warnings

import stdpopsim
import numpy as np
import msprime
import pyslim

logger = logging.getLogger(__name__)

_slim_upper = """
initialize() {
    if (!exists("dry_run"))
        defineConstant("dry_run", F);
    if (!exists("verbosity"))
        defineConstant("verbosity", 2);

    // Scaling factor to speed up simulation.
    // See SLiM manual:
    // `5.5 Rescaling population sizes to improve simulation performance`.
    defineConstant("Q", $scaling_factor);

    defineConstant("burn_in", $burn_in);
    defineConstant("generation_time", $generation_time);
    defineConstant("mutation_rate", Q * $mutation_rate);
    defineConstant("chromosome_length", $chromosome_length);
    defineConstant("trees_file", "$trees_file");
    defineConstant("pop_names", $pop_names);

    _recombination_rates = $recombination_rates;
    _recombination_ends = $recombination_ends;
    defineConstant("recombination_rates", (1-(1-2*_recombination_rates)^Q)/2);
    defineConstant("recombination_ends", _recombination_ends);
"""

_slim_lower = """
    defineConstant("N", asInteger(_N/Q));

    initializeTreeSeq();
    initializeMutationRate(mutation_rate);
    initializeMutationType("m1", 0.5, "f", 0);
    initializeGenomicElementType("g1", m1, 1.0);
    initializeGenomicElement(g1, 0, chromosome_length-1);
    initializeRecombinationRate(recombination_rates, recombination_ends);
}

function (void)err(string$ s) {
    stop("ERROR: " + s);
}

function (void)warn(string$ s) {
    catn("WARNING: " + s);
}

function (void)dbg(string$ s, [integer$ debug_level = 2]) {
    if (verbosity >= debug_level) {
        catn(sim.generation + ": " + s);
    }
}

// Check that sizes aren't dangerously low or zero (e.g. due to scaling).
function (void)check_size(integer$ pop, integer$ size, integer$ g) {
    if (size == 0) {
        err("The population size of p"+pop+" ("+pop_names[pop]+") is zero " +
            "at generation "+g+".");
    } else if (size < 50) {
        warn("p"+pop+" ("+pop_names[pop]+") has only "+size+" individuals " +
             "alive at generation "+g+".");
    }
}

// Return the epoch index for generation g.
function (integer)epoch(integer G, integer $g) {
    for (i in 0:(num_epochs-1)) {
        if (g < G[i]) {
            return i;
        }
    }
    return num_epochs - 1;
}

// Return the population size of pop at generation g.
function (integer)pop_size_at(integer G, integer$ pop, integer$ g) {
    e = epoch(G, g);
    N0 = N[e,pop];
    r = Q * growth_rates[e,pop];
    if (r == 0) {
        N_g = N0;
    } else {
        g_diff = g - G[e-1];
        N_g = asInteger(round(N0*exp(r*g_diff)));
    }
    return N_g;
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

1 {
    /*
     * Create initial populations and migration rates.
     */

    // Initial populations.
    for (i in 0:(num_populations-1)) {
        if (N[0,i] > 0) {
            check_size(i, N[0,i], sim.generation);
            dbg("sim.addSubpop("+i+", "+N[0,i]+");");
            sim.addSubpop(i, N[0,i]);
        }
    }

    if (length(sim.subpopulations) == 0) {
        err("No populations with non-zero size in generation 1.");
    }

    // Initial migration rates.
    i = 0;
    for (j in 0:(num_populations-1)) {
        for (k in 0:(num_populations-1)) {
            if (j==k | N[i,j] == 0 | N[i,k] == 0) {
                next;
            }

            m = Q * migration_matrices[k,j,i];
            p = sim.subpopulations[j];
            dbg("p"+j+".setMigrationRates("+k+", "+m+");");
            p.setMigrationRates(k, m);
        }
    }


    // The end of the burn-in is the starting generation, and corresponds to
    // time T_start. All remaining events are relative to this generation.
    N_max = max(N[0,0:(num_populations-1)]);
    G_start = sim.generation + asInteger(round(burn_in * N_max));
    T_start = max(_T);
    G = G_start + gdiff(T_start, _T);
    G_end = max(G);

    /*
     * Register events occurring at time T_start or more recently.
     */

    // Split events.
    if (length(subpopulation_splits) > 0 ) {
        for (i in 0:(ncol(subpopulation_splits)-1)) {
            g = G_start + gdiff(T_start, subpopulation_splits[0,i]);
            newpop = drop(subpopulation_splits[1,i]);
            size = asInteger(subpopulation_splits[2,i] / Q);
            oldpop = subpopulation_splits[3,i];
            check_size(newpop, size, g);
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
                // Change population size if this hasn't already been taken
                // care of by sim.addSubpop() or sim.addSubpopSplit().
                if (N[i,j] != N[i-1,j] & N[i-1,j] != 0) {
                    check_size(j, N[i,j], g);
                    sim.registerLateEvent(NULL,
                        "{dbg(self.source); " +
                        "p"+j+".setSubpopulationSize("+N[i,j]+");}",
                        g, g);
                }

                if (growth_rates[i,j] != 0) {
                    growth_phase_start = g+1;
                    if (i == num_epochs-1) {
                        growth_phase_end = G[i];
                    } else {
                        // We already registered a size change at generation G[i].
                        growth_phase_end = G[i] - 1;
                    }

                    if (growth_phase_start >= growth_phase_end) {
                        // Some demographic models have duplicate epoch times,
                        // which should be ignored.
                        next;
                    }

                    N_growth_phase_end = pop_size_at(G, j, growth_phase_end);
                    check_size(j, N_growth_phase_end, growth_phase_end);

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
                    if (j==k | N[i,j] == 0 | N[i,k] == 0) {
                        next;
                    }

                    m_last = Q * migration_matrices[k,j,i-1];
                    m = Q * migration_matrices[k,j,i];
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
            g = G_start + gdiff(T_start, admixture_pulses[0,i]);
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
        pop = drop(sampling_episodes[0,i]);
        n = sampling_episodes[1,i];
        g = G_start + gdiff(T_start, sampling_episodes[2,i]);

        // Check that there will be at least n individuals for sampling.
        N_g = pop_size_at(G, pop, g);
        if (n > N_g) {
            err("Request to sample "+n+" individuals from p"+pop+
                " ("+pop_names[pop]+") at generation "+g+", but only "+
                N_g+" individuals will be alive.");
        }

        sim.registerLateEvent(NULL,
            "{dbg(self.source); " +
            "inds=p"+pop+".sampleIndividuals("+n+"); " +
            "sim.treeSeqRememberIndividuals(inds);}",
            g, g);
    }

    sim.registerLateEvent(NULL, "{dbg(self.source); end();}", G_end, G_end);

    if (G_start > sim.generation) {
        dbg("Starting burn-in...");
    }

    if (dry_run) {
        sim.simulationFinished();
    }
}
"""


def msprime_rm_to_slim_rm(recombination_map):
    """
    Convert recombination map from start position coords to end position coords.
    """
    rates = recombination_map.get_rates()
    ends = [int(pos)-1 for pos in recombination_map.get_positions()]
    return rates[:-1], ends[1:]


def slim_makescript(
        script_file, trees_file,
        demographic_model, contig, samples,
        scaling_factor, burn_in):

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

    recomb_rates, recomb_ends = msprime_rm_to_slim_rm(contig.recombination_map)
    indent = 8*" "
    recomb_rates_str = (
            "c(\n" +
            textwrap.fill(
                    ", ".join(map(str, recomb_rates)),
                    width=80,
                    initial_indent=indent,
                    subsequent_indent=indent) +
            ")")
    recomb_ends_str = (
            "c(\n" +
            textwrap.fill(
                    ", ".join(map(str, recomb_ends)),
                    width=80,
                    initial_indent=indent,
                    subsequent_indent=indent) +
            ")")

    pop_names_str = ', '.join(map(lambda x: f'"{x}"', pop_names))

    printsc(string.Template(_slim_upper).substitute(
                scaling_factor=scaling_factor,
                burn_in=float(burn_in),
                chromosome_length=int(contig.recombination_map.get_length()),
                recombination_rates=recomb_rates_str,
                recombination_ends=recomb_ends_str,
                mutation_rate=contig.mutation_rate,
                generation_time=demographic_model.generation_time,
                trees_file=trees_file,
                pop_names=f"c({pop_names_str})"
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
    sample_counts = collections.Counter([
        (sample.population, round(sample.time * demographic_model.generation_time))
        for sample in samples])
    sampling_episodes = []
    for (pop, time), count in sample_counts.items():
        # SLiM can only sample individuals, which we assume are diploid.
        n_inds = (count+1) // 2
        if count % 2 != 0:
            pop_id = pop_names[pop]
            gen = time / demographic_model.generation_time
            warnings.warn(stdpopsim.SLiMOddSampleWarning(
                    f"SLiM simulates diploid individuals, so {n_inds} "
                    f"individuals will be sampled for the {count} haploids "
                    f"requested from population {pop_id} at time {gen}. "
                    "See #464."))
        sampling_episodes.append((pop, n_inds, time))

    printsc('    // One row for each sampling episode.')
    printsc('    defineConstant("sampling_episodes", ' +
            matrix2str(
                sampling_episodes,
                col_comment='pop, n_inds, time') +
            ');')

    printsc(_slim_lower)

    return epochs[0]


class SLiMException(Exception):
    pass


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
            slim_path=None, slim_script=False, slim_scaling_factor=1.0,
            slim_burn_in=10.0, dry_run=False, **kwargs):
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
        :param slim_burn_in: Length of the burn-in phase, in units of N
            generations.
        :type slim_burn_in: float
        :param dry_run: If True, run the first generation setup and then end the
            simulation.
        :type dry_run: bool
        """

        if slim_scaling_factor <= 0:
            raise ValueError("slim_scaling_factor must be positive")
        if slim_burn_in < 0:
            raise ValueError("slim_burn_in must be non-negative")

        if slim_scaling_factor != 1:
            warnings.warn(stdpopsim.SLiMScalingFactorWarning(
                f"You're using a scaling factor ({slim_scaling_factor}). "
                "This should give similar results for many situations, "
                "but is not equivalent, especially in the presence of selection. "
                "When using rescaling, you should be careful---do checks and "
                "compare results across different values of the scaling factor."))

        run_slim = not slim_script

        mutation_rate = contig.mutation_rate
        # Ensure no mutations are introduced by SLiM.
        contig = stdpopsim.Contig(
                recombination_map=contig.recombination_map,
                mutation_rate=0,
                genetic_map=contig.genetic_map)

        mktemp = functools.partial(tempfile.NamedTemporaryFile, mode="w")

        @contextlib.contextmanager
        def script_file_f():
            f = mktemp(suffix=".slim") if not slim_script else sys.stdout
            yield f
            # Don't close sys.stdout.
            if not slim_script:
                f.close()

        with script_file_f() as script_file, mktemp(suffix=".ts") as ts_file:

            recap_epoch = slim_makescript(
                    script_file, ts_file.name,
                    demographic_model, contig, samples,
                    slim_scaling_factor, slim_burn_in)

            script_file.flush()

            if not run_slim:
                return None

            self._run_slim(
                    script_file.name, slim_path=slim_path, seed=seed,
                    dry_run=dry_run)

            if dry_run:
                return None

            ts = pyslim.load(ts_file.name)

        ts = self._recap_and_rescale(
                ts, seed, recap_epoch, contig, mutation_rate, slim_scaling_factor)
        return ts

    def _run_slim(self, script_file, slim_path=None, seed=None, dry_run=False):
        """
        Run SLiM.

        We capture the output using Popen's line-oriented text buffering
        (bufsize=1, universal_newlines=True) and redirect all messages to
        Python's logging module.
        By convention, messages from SLiM prefixed with "ERROR: " or
        "WARNING: " are treated as ERROR or WARN loglevels respectively.
        All other output on stdout is given the DEBUG loglevel.
        ERROR messages, and any output from SLiM on stderr, will raise a
        SLiMException here.
        """
        if slim_path is None:
            slim_path = self.slim_path()
        slim_cmd = [slim_path]
        if seed is not None:
            slim_cmd.extend(["-s", f"{seed}"])
        if dry_run:
            slim_cmd.extend(["-d", "dry_run=T"])
        slim_cmd.append(script_file)

        with subprocess.Popen(
                slim_cmd, bufsize=1, universal_newlines=True,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE) as proc:
            for line in proc.stdout:
                line = line.rstrip()
                if line.startswith("ERROR: "):
                    logger.error(line[len("ERROR: "):])
                elif line.startswith("WARNING: "):
                    warnings.warn(stdpopsim.UnspecifiedSLiMWarning(
                        line[len("WARNING: "):]))
                else:
                    # filter `dbg` function calls that generate output
                    line = line.replace("dbg(self.source); ", "")
                    logger.debug(line)
            stderr = proc.stderr.read()

        if proc.returncode != 0 or stderr:
            raise SLiMException(
                    f"{slim_path} exited with code {proc.returncode}.\n"
                    f"{stderr}")

    def _simplify_remembered(self, ts):
        """
        Remove all samples except those individuals that were explicity
        sampled in SLiM with sim.treeSeqRememberIndividuals().
        """
        nodes = itertools.chain.from_iterable(
                    i.nodes for i in ts.individuals()
                    if i.flags & pyslim.INDIVIDUAL_REMEMBERED)
        return ts.simplify(samples=list(nodes), filter_populations=False)

    def _recap_and_rescale(
            self, ts, seed, recap_epoch, contig, mutation_rate, slim_scaling_factor):
        """
        Apply post-SLiM transformations to ``ts``. This rescales node times,
        does recapitation, simplification, and adds neutral mutations.
        """
        # Node times come from SLiM generation numbers, which may have been
        # divided by a scaling factor for computational tractability.
        tables = ts.dump_tables()
        for table in (tables.nodes, tables.migrations):
            table.time *= slim_scaling_factor
        ts = pyslim.SlimTreeSequence.load_tables(tables)
        ts.slim_generation *= slim_scaling_factor

        rng = random.Random(seed)
        s1, s2 = rng.randrange(1, 2**32), rng.randrange(1, 2**32)

        population_configurations = [
                msprime.PopulationConfiguration(
                    initial_size=pop.start_size,
                    growth_rate=pop.growth_rate)
                for pop in recap_epoch.populations]
        ts = ts.recapitate(
                recombination_rate=contig.recombination_map.mean_recombination_rate,
                population_configurations=population_configurations,
                migration_matrix=recap_epoch.migration_matrix,
                random_seed=s1)

        ts = self._simplify_remembered(ts)

        # Add neutral mutations.
        ts = pyslim.SlimTreeSequence(msprime.mutate(
            ts, rate=mutation_rate, keep=True, random_seed=s2))

        return ts

    def recap_and_rescale(
            self, ts, demographic_model, contig, samples,
            slim_scaling_factor=1.0, seed=None, **kwargs):
        """
        Apply post-SLiM transformations to ``ts``. This rescales node times,
        does recapitation, simplification, and adds neutral mutations.

        If the SLiM engine was used to output a SLiM script, and the script was
        run outside of stdpopsim, this function can be used to transform the
        SLiM tree sequence following the procedure that would have been used
        if stdpopsim had run SLiM itself.
        The parameters after ``ts`` have the same meaning as for :func:`simulate`,
        and the values for ``demographic_model``, ``contig``, ``samples``,
        and ``slim_scaling_factor`` should match those that were used to
        generate the SLiM script with :func:`simulate`.

        :param ts: The tree sequence output by SLiM.
        :type ts: :class:`pyslim.SlimTreeSequence`

        .. warning::
            The :func:`recap_and_rescale` function is provided in the hope that
            it will be useful. But as we can't anticipate what changes you'll
            make to the SLiM code before using it, the stdpopsim source code
            should be consulted to determine if it's behaviour is appropriate
            for your case.
        """
        with open(os.devnull, "w") as script_file:
            recap_epoch = slim_makescript(
                    script_file, "unused.trees",
                    demographic_model, contig, samples,
                    slim_scaling_factor, 1)

        ts = self._recap_and_rescale(
                ts, seed, recap_epoch, contig, contig.mutation_rate, slim_scaling_factor)
        return ts


# SLiM does not currently work on Windows.
if sys.platform != "win32":
    stdpopsim.register_engine(_SLiMEngine())
