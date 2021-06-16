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
import copy
import string
import tempfile
import subprocess
import functools
import itertools
import contextlib
import random
import textwrap
import logging
import warnings

import stdpopsim
import numpy as np
import msprime
import pyslim
import tskit

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
    // save/restore bookkeeping
    sim.setValue("n_restores", 0);
    sim.setValue("n_saves", 0);
    sim.setValue("restore_function", F);

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

    // Save/restore events. These should come before all other events.
    if (length(drawn_mutations) > 0) {
        for (i in 0:(ncol(drawn_mutations)-1)) {
            save = drawn_mutations[4,i] == 1;
            if (!save) {
                next;
            }
            g = G_start + gdiff(T_start, drawn_mutations[0,i]);
            // Unconditionally save the state before the mutation is drawn.
            sim.registerLateEvent(NULL, "{save();}", g, g);
        }
    }
    if (length(condition_on_allele_frequency) > 0) {
        for (i in 0:(ncol(condition_on_allele_frequency)-1)) {
            g_start = G_start + gdiff(T_start, condition_on_allele_frequency[0,i]);
            g_end = G_start + gdiff(T_start, condition_on_allele_frequency[1,i]);
            mut_type = asInteger(condition_on_allele_frequency[2,i]);
            pop_id = asInteger(condition_on_allele_frequency[3,i]);
            op = op_types[asInteger(drop(condition_on_allele_frequency[4,i]))];
            af = condition_on_allele_frequency[5,i];
            save = condition_on_allele_frequency[6,i] == 1;

            if (g_start > g_end) {
                err("Attempt to register AF conditioning callback with g_start="+
                    g_start+" > g_end="+g_end);
            }

            if (save) {
                // Save the state conditional on the allele frequency.
                // If the condition isn't met, we restore.
                sim.registerLateEvent(NULL,
                    "{if (af(m"+mut_type+", p"+pop_id+") "+op+" "+af+")" +
                    " save(); else restore();}",
                    g_start, g_start);
                g_start = g_start + 1;
            }

            if (g_start <= g_end) {
                sim.registerLateEvent(NULL,
                    "{if (!(af(m"+mut_type+", p"+pop_id+") "+op+" "+af+"))" +
                    " restore();}",
                    g_start, g_end);
            }
        }
    }

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
                        // Demographic models could have duplicate epoch times,
                        // which should be fixed.
                        warn("growth_phase_start="+growth_phase_start+
                             " >= growth_phase_end="+growth_phase_end);
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

    // Draw mutations.
    if (length(drawn_mutations) > 0) {
        for (i in 0:(ncol(drawn_mutations)-1)) {
            g = G_start + gdiff(T_start, drawn_mutations[0,i]);
            mut_type = drawn_mutations[1,i];
            pop_id = drawn_mutations[2,i];
            coordinate = drawn_mutations[3,i];
            sim.registerLateEvent(NULL,
                "{dbg(self.source); " +
                "add_mut(m"+mut_type+", p"+pop_id+", "+coordinate+");}",
                g, g);
        }
    }

    // Setup fitness callbacks.
    if (length(fitness_callbacks) > 0) {
        for (i in 0:(ncol(fitness_callbacks)-1)) {
            g_start = G_start + gdiff(T_start, fitness_callbacks[0,i]);
            g_end = G_start + gdiff(T_start, fitness_callbacks[1,i]);
            mut_type = asInteger(fitness_callbacks[2,i]);
            pop_id = asInteger(fitness_callbacks[3,i]);
            selection_coeff = Q * fitness_callbacks[4,i];
            dominance_coeff = fitness_callbacks[5,i];

            if (g_start > g_end) {
                err("Attempt to register fitness callback with g_start="+
                    g_start+" > g_end="+g_end);
            }

            sim.registerLateEvent(NULL,
                "{dbg('s="+selection_coeff+", h="+dominance_coeff+
                " for m"+mut_type+" in p"+pop_id+"');}",
                g_start, g_start);
            sim.registerLateEvent(NULL,
                "{dbg('s, h defaults for m"+mut_type+" in p"+pop_id+"');}",
                g_end, g_end);
            /* We explicitly format() here to prevent integral-valued floats
             * from getting converted to integers during string interpolation
             * (this triggers a type error when the fitness callback runs). */
            f_hom = format("%e", 1 + selection_coeff);
            f_het = format("%e", 1 + selection_coeff * dominance_coeff);
            sim.registerFitnessCallback(NULL,
                "{if (homozygous) return "+f_hom+"; else return "+f_het+";}",
                mut_type, pop_id, g_start, g_end);
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

// Add `mut_type` mutation at `pos`, to a single individual in `pop`.
function (void)add_mut(object$ mut_type, object$ pop, integer$ pos) {
   targets = sample(pop.genomes, 1);
   targets.addNewDrawnMutation(mut_type, pos);
}

// Return the allele frequency of a drawn mutation in the specified population.
// Assumes there's only one mutation of the given type.
function (float$)af(object$ mut_type, object$ pop) {
    mut = sim.mutationsOfType(mut_type);
    if (length(mut) == 0) {
        return 0.0;
    }
    return sim.mutationFrequencies(pop, mut);
}

// Save the state of the simulation.
function (void)save(void) {
    if (sim.getValue("restore_function")) {
        // Don't save if we're in the restore() function.
        return;
    }
    n_saves = 1 + sim.getValue("n_saves");
    sim.setValue("n_saves", n_saves);
    dbg("save() "+n_saves);
    sim.treeSeqOutput(trees_file);
}

// Restore the simulation state.
function (void)restore(void) {
    g_restore = sim.generation;
    n_restores = 1 + sim.getValue("n_restores");
    sim.setValue("n_restores", n_restores);
    n_saves = sim.getValue("n_saves");
    if (n_saves == 0) {
        err("restore() in generation "+g_restore+", but nothing is saved.");
    }
    sim.readFromPopulationFile(trees_file);
    dbg("restore() "+n_restores+" from generation "+g_restore+", returning "+
        "to state at save() "+n_saves);

    /*
     * The generation counter sim.generation has now been reset to the
     * value it had when save() was called. There are two issues relating
     * to event scheduling which must now be dealt with.
     *
     * 1. There may be additional late events for the generation in which
     * restore() was called, and they are still scheduled to run.
     * So we deactivate all script blocks to avoid unexpected problems.
     * They will be automatically reactivated at the start of the next
     * generation (see SLiM manual section 23.10).
     */
    sim.scriptBlocks.active = F;

    /*
     * 2. The late events below were run in the save() generation,
     * but after the save() call. We execute these again here, because
     * the next late events to run will be for sim.generation + 1.
     * Note that the save() event is indistinguishable from the other
     * late events in this generation, so we set a flag `restore_function`
     * to signal the save() function not to save again.
     */
    g = sim.generation;
    sim.setValue("restore_function", T);
    for (sb in sim.scriptBlocks) {
        if (sb.type == "late" & g >= sb.start & g <= sb.end) {
            self = sb;
            executeLambda(sb.source);
        }
    }
    sim.setValue("restore_function", F);
}

"""


def msprime_rm_to_slim_rm(recombination_map):
    """
    Convert recombination map from start position coords to end position coords.
    """
    rates = recombination_map.rate.copy()
    # replace missing values with 0 recombination rate
    rates[recombination_map.missing] = 0
    ends = [int(pos) - 1 for pos in recombination_map.position]
    return rates, ends[1:]


def get_slim_fractions(contig):
    return np.array([sum(g.proportions) for g in contig.genomic_element_types])


def compute_msp_mutation_rate_map(contig, slim_fractions):
    """
    breaks = [0,10,20] means [0,10) and [10,20)
    """
    length = contig.recombination_map.sequence_length
    breaks = [0]
    rates = []
    for (start, end, t) in contig.all_intervals_array:
        if start not in breaks:
            breaks.append(start)
            rates.append(contig.mutation_rate)
        breaks.append(end)
        rates.append(
            np.around((1 - slim_fractions[t]) * contig.mutation_rate, decimals=16)
        )
    if not np.isclose(breaks[-1], int(length)):
        breaks.append(int(length))
        rates.append(contig.mutation_rate)
    assert len(breaks) == len(rates) + 1
    return msprime.RateMap(position=breaks, rate=rates)


def compute_slim_mutation_breaks_and_rates(contig, msp_mutation_rate_map):
    """
    breaks = [10,20] means [0,10] and (10,20]
    """
    # making sure the mut types are still concordant with the ge types
    breaks = msp_mutation_rate_map.position
    rates = msp_mutation_rate_map.rate
    assert breaks[0] == 0
    slim_breaks = [b - 1 for b in breaks[1:]]
    slim_rates = [contig.mutation_rate - r for r in rates]
    return (slim_breaks, slim_rates)


def get_msp_and_slim_mutation_rate_maps(contig):
    """
    Returns a :class:`msprime.RateMap` with the mutation rate map for overlay
    of neutral mutations and a tuple with the breakpoints and rates of mutations
    for SLiM.
    """
    slim_fractions = get_slim_fractions(contig)
    msp_mutation_rate_map = compute_msp_mutation_rate_map(contig, slim_fractions)

    slim_breaks, slim_rates = compute_slim_mutation_breaks_and_rates(
        contig, msp_mutation_rate_map
    )
    return (msp_mutation_rate_map, (slim_breaks, slim_rates))


def slim_array_string(iterable, indent, width=80):
    """
    Format an array as a SLiM c() array and return as a line-wrapped string.
    """
    return (
        "c(\n"
        + textwrap.fill(
            ", ".join(map(str, iterable)),
            width=width,
            initial_indent=indent,
            subsequent_indent=indent,
        )
        + ")"
    )


def slim_makescript(
    script_file,
    trees_file,
    demographic_model,
    contig,
    samples,
    extended_events,
    scaling_factor,
    burn_in,
    slim_rate_map,
):
    mutation_types = contig.mutation_types
    pop_names = [pop.name for pop in demographic_model.model.populations]
    # Use copies of these so that the time frobbing below doesn't have
    # side-effects in the caller's model.
    demographic_events = copy.deepcopy(demographic_model.model.events)
    if extended_events is None:
        extended_events = []
    else:
        extended_events = copy.deepcopy(extended_events)

    # Reassign event times according to integral SLiM generations.
    # This collapses the time deltas used in HomSap/AmericanAdmixture_4B11,
    # and calculates times for GenerationAfter objects.
    def fix_time(event):
        for attr in ("time", "start_time", "end_time"):
            if not hasattr(event, attr):
                continue
            t = getattr(event, attr)
            t_rounded = round(float(t) / scaling_factor) * scaling_factor
            if isinstance(t, stdpopsim.ext.GenerationAfter):
                t_rounded -= scaling_factor
            if t_rounded < 0:
                raise ValueError(f"Bad {attr}: {getattr(event, attr)}")
            setattr(event, attr, t_rounded)

    for event in demographic_events:
        fix_time(event)
    for event in extended_events:
        fix_time(event)

    # The demography debugger constructs event epochs, which we use
    # to define the forwards-time events.
    dd = demographic_model.model.debug()
    # msprime.DemographyDebugger(
    #     population_configurations=demographic_model.population_configurations,
    #     migration_matrix=demographic_model.migration_matrix,
    #     demographic_events=demographic_events,
    # )

    epochs = sorted(dd.epochs, key=lambda e: e.start_time, reverse=True)
    T = [round(e.start_time * demographic_model.generation_time) for e in epochs]
    migration_matrices = [e.migration_matrix for e in epochs]

    N = np.empty(shape=(dd.num_populations, len(epochs)), dtype=int)
    growth_rates = np.empty(shape=(dd.num_populations, len(epochs)), dtype=float)
    for j, epoch in enumerate(epochs):
        for i, pop in enumerate(epoch.populations):
            N[i, j] = int(pop.end_size)
            growth_rates[i, j] = pop.growth_rate

    admixture_pulses = []
    subpopulation_splits = []
    for i, epoch in enumerate(epochs):
        for de in epoch.demographic_events:
            if isinstance(de, msprime.demography.LineageMovementEvent):
                # This is using internal msprime APIs here, but it's not worth
                # updating before we change over to Demes.
                for lm in de._as_lineage_movements():
                    if lm.proportion < 1:
                        # Calculate remainder of population after previous
                        # MassMigration events in this epoch.
                        rem = 1 - np.sum(
                            [
                                ap[3]
                                for ap in admixture_pulses
                                if ap[0] == i and ap[1] == lm.source
                            ]
                        )
                        admixture_pulses.append(
                            (
                                i,
                                lm.source,  # forwards-time dest
                                lm.dest,  # forwards-time source
                                rem * lm.proportion,
                            )
                        )
                        continue

                    # Backwards: lm.source is being merged into lm.dest.
                    # Forwards: lm.source is being created, taking individuals
                    #           from lm.dest.
                    #
                    # If the proportion==1, we can use SLiM function:
                    #       sim.addSubpopSplit(newpop, size, oldpop),
                    # which we trigger by adding a row to subpopulation_splits.
                    # This SLiM function creates newpop (=lm.source), under the
                    # assumption that it doesn't already exist.

                    subpopulation_splits.append(
                        (f"_T[{i}]", lm.source, f"_N[{i+1},{lm.source}]", lm.dest)
                    )

                    # Zero out the population size for generations before this
                    # epoch, to avoid simulating invididuals that contribute no
                    # genealogy.
                    N[lm.source, 0 : (i + 1)] = 0
                    growth_rates[lm.source, 0 : (i + 1)] = 0

                    # Ensure there are no migrations to or from lm.source before
                    # this epoch.
                    for j in range(i + 1):
                        for k in range(dd.num_populations):
                            migration_matrices[j][k][lm.source] = 0
                            migration_matrices[j][lm.source][k] = 0

    drawn_mutations = []
    fitness_callbacks = []
    condition_on_allele_frequency = []
    op_id = stdpopsim.ext.ConditionOnAlleleFrequency.op_id
    for ee in extended_events:
        if hasattr(ee, "mutation_type_id"):
            mt_id = getattr(ee, "mutation_type_id")
            cls_name = ee.__class__.__name__
            if not (0 <= mt_id < len(mutation_types)):
                raise ValueError(
                    f"Invalid {cls_name} {ee.mutation_type_id} {len(mutation_types)} "
                    f"event with mutation type id {mt_id}.{contig}"
                )
        if hasattr(ee, "start_time") and hasattr(ee, "end_time"):
            # Now that GenerationAfter times have been accounted for, we can
            # properly catch invalid start/end times.
            start_time = getattr(ee, "start_time")
            end_time = getattr(ee, "end_time")
            stdpopsim.ext.validate_time_range(start_time, end_time)

        if isinstance(ee, stdpopsim.ext.DrawMutation):
            # TODO: warning if mutation type id points to a neutral mut.
            time = ee.time * demographic_model.generation_time
            save = 1 if ee.save else 0
            drawn_mutations.append(
                (time, ee.mutation_type_id, ee.population_id, ee.coordinate, save)
            )
        elif isinstance(ee, stdpopsim.ext.ChangeMutationFitness):
            start_time = ee.start_time * demographic_model.generation_time
            end_time = ee.end_time * demographic_model.generation_time
            fitness_callbacks.append(
                (
                    start_time,
                    end_time,
                    ee.mutation_type_id,
                    ee.population_id,
                    ee.selection_coeff,
                    ee.dominance_coeff,
                )
            )
        elif isinstance(ee, stdpopsim.ext.ConditionOnAlleleFrequency):
            start_time = ee.start_time * demographic_model.generation_time
            end_time = ee.end_time * demographic_model.generation_time
            save = 1 if ee.save else 0
            condition_on_allele_frequency.append(
                (
                    start_time,
                    end_time,
                    ee.mutation_type_id,
                    ee.population_id,
                    op_id(ee.op),
                    ee.allele_frequency,
                    save,
                )
            )
        else:
            raise ValueError(f"Unknown extended event type {type(ee)}")

    # Check that drawn mutations exist for extended events that need them.
    drawn_mut_type_ids = {mt_id for _, mt_id, _, _, _ in drawn_mutations}
    for ee in extended_events:
        if isinstance(ee, stdpopsim.ext.ChangeMutationFitness) or isinstance(
            ee, stdpopsim.ext.ConditionOnAlleleFrequency
        ):
            if ee.mutation_type_id not in drawn_mut_type_ids:
                cls_name = ee.__class__.__name__
                raise ValueError(
                    f"Invalid {cls_name} event. No drawn mutation for "
                    f"mutation type id {ee.mutation_type_id}"
                )

    printsc = functools.partial(print, file=script_file)

    # Header
    printsc("/*")
    printsc(" * stdpopsim " + stdpopsim.__version__)
    printsc(" *")
    printsc(" * Demographic model: " + demographic_model.id)
    printsc(
        " * "
        + "\n * ".join(
            [line.strip() for line in demographic_model.description.split("\n")]
        )
    )
    for citation in demographic_model.citations:
        printsc(" * " + str(citation))
    printsc(" */")

    recomb_rates, recomb_ends = msprime_rm_to_slim_rm(contig.recombination_map)
    indent = 8 * " "
    recomb_rates_str = slim_array_string(recomb_rates, indent)
    recomb_ends_str = slim_array_string(recomb_ends, indent)

    pop_names_str = ", ".join(map(lambda x: f'"{x}"', pop_names))

    printsc(
        string.Template(_slim_upper).substitute(
            scaling_factor=scaling_factor,
            burn_in=float(burn_in),
            recombination_rates=recomb_rates_str,
            recombination_ends=recomb_ends_str,
            generation_time=demographic_model.generation_time,
            trees_file=trees_file,
            pop_names=f"c({pop_names_str})",
        )
    )

    def matrix2str(
        matrix, row_comments=None, col_comment=None, indent=2, fmt="", dim=(None, None)
    ):
        """
        Return an Eidos representation of the matrix as a string.
        """
        if row_comments is not None:
            assert len(matrix) == len(row_comments)

        if len(matrix) == 0:
            return "c()"

        s = ["array(c(\n"]
        if col_comment is not None:
            s.append(indent * 4 * " " + "// " + col_comment + "\n")

        for i in range(len(matrix)):
            s.append(indent * 4 * " ")
            s.append("c({})".format(", ".join([format(x, fmt) for x in matrix[i]])))
            if i != len(matrix) - 1:
                s.append(",")
            if row_comments is not None:
                s.append(" // " + row_comments[i])
            s.append("\n")

        s.append((indent - 1) * 4 * " ")

        if dim[0] is None:
            dim = (len(matrix[0]), dim[1])
        if dim[1] is None:
            dim = (dim[0], len(matrix))
        s.append(f"), c({dim[0]}, {dim[1]}))")

        return "".join(s)

    # Mutation types.
    for i, m in enumerate(mutation_types):
        distrib_args = [str(arg) for arg in m.distribution_args]
        distrib_args[m.Q_scaled_index] = "Q * " + distrib_args[m.Q_scaled_index]
        distrib_args = ", ".join(distrib_args)
        printsc(
            f"    initializeMutationType({i}, {m.dominance_coeff}, "
            + f'"{m.distribution_type}", {distrib_args});'
        )
        if not m.convert_to_substitution:
            # T is the default for WF simulations.
            printsc(f"    m{i}.convertToSubstitution = F;")
    # Genomic element types.
    genomic_element_types = contig.genomic_element_types
    for j, g in enumerate(genomic_element_types):
        mut_props = ", ".join([str(prop) for prop in g.proportions])
        mut_types = ", ".join([str(mid) for mid in g.mutation_type_ids])
        element_starts = slim_array_string(g.intervals[:, 0], indent)
        # stdpopsim intervals are 0-based left inclusive, right exclusive, but
        # SLiM intervals are right inclusive
        element_ends = slim_array_string(g.intervals[:, 1] - 1, indent)
        printsc(
            f"    initializeGenomicElementType({j}, c({mut_types}), c({mut_props}));"
        )
        printsc(f"    initializeGenomicElement({j}, {element_starts}, {element_ends});")
    # Mutation rate map.
    mut_breaks = slim_array_string(map(int, slim_rate_map[0]), indent)
    mut_rates = slim_array_string(slim_rate_map[1], indent)
    printsc(f"    initializeMutationRate(Q*{mut_rates}, {mut_breaks});")
    printsc()

    # Epoch times.
    printsc("    // Time of epoch boundaries, in years before present.")
    printsc("    // The first epoch spans from INF to _T[0].")
    printsc('    defineConstant("_T", c({}));'.format(", ".join(map(str, T))))
    printsc()

    # Population sizes.
    printsc("    // Population sizes in each epoch.")
    printsc(
        "    _N = "
        + matrix2str(
            N, row_comments=pop_names, col_comment="INF:_T[0], _T[0]:_T[1], etc."
        )
        + ";"
    )
    printsc()

    printsc('    defineConstant("num_epochs", length(_T));')
    printsc('    defineConstant("num_populations", ncol(_N));')
    printsc()

    # Growth rates.
    printsc("    // Population growth rates for each epoch.")
    printsc(
        '    defineConstant("growth_rates", '
        + matrix2str(
            growth_rates,
            row_comments=pop_names,
            col_comment="INF:_T[0], _T[0]:_T[1], etc.",
            dim=("num_epochs", "num_populations"),
        )
        + ");"
    )
    printsc()

    printsc("    no_migration = rep(0, num_populations*num_populations);")
    printsc()

    # Migration rates.
    printsc("    // Migration rates for each epoch.")
    printsc("    // Migrations involving a population with size=0 are ignored.")
    printsc("    // XXX: document what the rows & cols correspond to.")
    printsc('    defineConstant("migration_matrices", array(c(')
    for i in range(len(migration_matrices)):
        epoch_str = f"INF:_T[{i}]" if i == 0 else f"_T[{i}]:_T[{i+1}]"
        printsc()
        printsc(2 * 4 * " " + "// " + epoch_str)

        end = ",\n" if i != len(migration_matrices) - 1 else "\n"
        if np.all(np.array(migration_matrices[i]) == 0):
            printsc(2 * 4 * " " + "no_migration", end=end)
        else:
            printsc(
                2 * 4 * " "
                + matrix2str(
                    migration_matrices[i],
                    indent=3,
                    fmt="g",
                    dim=("num_populations", "num_populations"),
                ),
                end=end,
            )
    printsc()
    printsc(4 * " " + "), c(num_populations, num_populations, num_epochs)));")
    printsc()

    # Population splits.
    printsc("    // Population splits, one row for each event.")
    printsc(
        '    defineConstant("subpopulation_splits", '
        + matrix2str(subpopulation_splits, col_comment="time, newpop, size, oldpop")
        + ");"
    )
    printsc()

    # Admixture pulses.
    # Output _T[...] variable rather than an index.
    admixture_pulses = [(f"_T[{ap[0]}]", *ap[1:]) for ap in admixture_pulses]
    printsc("    // Admixture pulses, one row for each pulse.")
    printsc(
        '    defineConstant("admixture_pulses", '
        + matrix2str(admixture_pulses, col_comment="time, dest, source, rate")
        + ");"
    )
    printsc()

    # Drawn mutations.
    printsc("    // Drawn mutations, one row for each mutation.")
    printsc(
        '    defineConstant("drawn_mutations", '
        + matrix2str(
            drawn_mutations,
            col_comment="time, mut_type, pop_id, genomic_coordinate, save",
        )
        + ");"
    )
    printsc()

    # Fitness callbacks.
    printsc("    // Fitness callbacks, one row for each callback.")
    printsc(
        '    defineConstant("fitness_callbacks", '
        + matrix2str(
            fitness_callbacks,
            col_comment="start_time, end_time, mut_type, pop_id, "
            "selection_coeff, dominance_coeff",
        )
        + ");"
    )
    printsc()

    # Allele frequency conditioning
    op_types = ", ".join(
        f'"{op}"' for op in stdpopsim.ext.ConditionOnAlleleFrequency.op_types
    )
    printsc(f'    defineConstant("op_types", c({op_types}));')
    printsc("    // Allele frequency conditioning, one row for each.")
    printsc(
        '    defineConstant("condition_on_allele_frequency", '
        + matrix2str(
            condition_on_allele_frequency,
            col_comment="start_time, end_time, mut_type, pop_id, "
            "op, allele_frequency, save",
        )
        + ");"
    )
    printsc()

    # Sampling episodes.
    sampling_episodes = []
    for sample_set in samples:
        pop = demographic_model.model[sample_set.population]
        # We're currently sampling genomes, and set the ploidy to 1 to make
        # sure it all fits together. This is probably confusing though and
        # maybe we should change it.
        assert sample_set.ploidy == 1
        # SLiM can only sample individuals, which we assume are diploid.
        count = sample_set.num_samples
        time = 0 if sample_set.time is None else sample_set.time
        time = round(time * demographic_model.generation_time)
        n_inds = (count + 1) // 2
        if count % 2 != 0:
            gen = time / demographic_model.generation_time
            warnings.warn(
                stdpopsim.SLiMOddSampleWarning(
                    f"SLiM simulates diploid individuals, so {n_inds} "
                    f"individuals will be sampled for the {count} haploids "
                    f"requested from population {pop.name} at time {gen}. "
                    "See #464."
                )
            )
        sampling_episodes.append((pop.id, n_inds, time))

    printsc("    // One row for each sampling episode.")
    printsc(
        '    defineConstant("sampling_episodes", '
        + matrix2str(sampling_episodes, col_comment="pop, n_inds, time")
        + ");"
    )

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
            reasons={stdpopsim.CiteReason.ENGINE},
        ),
    ]

    def slim_path(self):
        return os.environ.get("SLIM", "slim")

    def get_version(self, slim_path=None):
        if slim_path is None:
            slim_path = self.slim_path()
        s = subprocess.check_output([slim_path, "-v"])
        return s.split()[2].decode("ascii").rstrip(",")

    def _assert_min_version(self, min_required_version, slim_path):
        def version_split(version):
            return [int(v) for v in version.split(".")]

        current_version = self.get_version(slim_path)
        if version_split(current_version) < version_split(min_required_version):
            raise RuntimeError(
                f"Minimum supported SLiM version is {min_required_version}, "
                f"but only found version {current_version}"
            )

    def simulate(
        self,
        demographic_model,
        contig,
        samples,
        *,
        seed=None,
        extended_events=None,
        slim_path=None,
        slim_script=False,
        slim_scaling_factor=1.0,
        slim_burn_in=10.0,
        dry_run=False,
    ):
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
            warnings.warn(
                stdpopsim.SLiMScalingFactorWarning(
                    f"You're using a scaling factor ({slim_scaling_factor}). "
                    "This should give similar results for many situations, "
                    "but is not equivalent, especially in the presence of selection. "
                    "When using rescaling, you should be careful---do checks and "
                    "compare results across different values of the scaling factor."
                )
            )

        # figuring out mutations in SLiM and neutral mutations to be added
        # by msprime later
        msp_rate_map, slim_rate_map = get_msp_and_slim_mutation_rate_maps(contig)

        # TODO: remove this after a release or two. See #745.
        self._warn_zigzag(demographic_model)
        self._warn_mutation_rate_mismatch(contig, demographic_model)

        run_slim = not slim_script

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
                script_file,
                ts_file.name,
                demographic_model,
                contig,
                samples,
                extended_events,
                slim_scaling_factor,
                slim_burn_in,
                slim_rate_map,
            )

            script_file.flush()

            if not run_slim:
                return None

            self._run_slim(
                script_file.name, slim_path=slim_path, seed=seed, dry_run=dry_run
            )

            if dry_run:
                return None

            ts = tskit.load(ts_file.name)

        ts = self._recap_and_rescale(
            ts, seed, recap_epoch, contig, slim_scaling_factor, msp_rate_map
        )

        if contig.inclusion_mask is not None:
            ts = stdpopsim.utils.mask_tree_sequence(ts, contig.inclusion_mask, False)
        if contig.exclusion_mask is not None:
            ts = stdpopsim.utils.mask_tree_sequence(ts, contig.exclusion_mask, True)

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
        ERROR messages will raise a SLiMException here too, because
        they are always generated by the `stop()` eidos function which
        makes SLiM exit with a non-zero return code.
        """
        if slim_path is None:
            slim_path = self.slim_path()

        # SLiM v3.6 sends `stop()` output to stderr, which we rely upon.
        self._assert_min_version("3.6", slim_path)

        slim_cmd = [slim_path]
        if seed is not None:
            slim_cmd.extend(["-s", f"{seed}"])
        if dry_run:
            slim_cmd.extend(["-d", "dry_run=T"])
        slim_cmd.append(script_file)

        with subprocess.Popen(
            slim_cmd,
            bufsize=1,
            universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ) as proc:
            for line in proc.stdout:
                line = line.rstrip()
                if line.startswith("WARNING: "):
                    warnings.warn(
                        stdpopsim.UnspecifiedSLiMWarning(line[len("WARNING: ") :])
                    )
                else:
                    # filter `dbg` function calls that generate output
                    line = line.replace("dbg(self.source); ", "")
                    logger.debug(line)

            stderr = proc.stderr.read()
            for line in stderr.splitlines():
                if line.startswith("ERROR: "):
                    logger.error(line[len("ERROR: ") :])

        if proc.returncode != 0:
            raise SLiMException(
                f"{slim_path} exited with code {proc.returncode}.\n{stderr}"
            )

    def _simplify_remembered(self, ts):
        """
        Remove all samples except those individuals that were explicity
        sampled in SLiM with sim.treeSeqRememberIndividuals().
        """
        nodes = itertools.chain.from_iterable(
            i.nodes for i in ts.individuals() if i.flags & pyslim.INDIVIDUAL_REMEMBERED
        )
        return ts.simplify(samples=list(nodes), filter_populations=False)

    def _recap_and_rescale(
        self, ts, seed, recap_epoch, contig, slim_scaling_factor, msp_rate_map
    ):
        """
        Apply post-SLiM transformations to ``ts``. This rescales node times,
        does recapitation, simplification, and adds neutral mutations.
        """
        # Node times come from SLiM generation numbers, which may have been
        # divided by a scaling factor for computational tractability.
        tables = ts.dump_tables()
        for table in (tables.nodes, tables.migrations, tables.mutations):
            table.time *= slim_scaling_factor
        ts_metadata = tables.metadata
        ts_metadata["SLiM"]["generation"] *= slim_scaling_factor
        tables.metadata = ts_metadata
        ts = tables.tree_sequence()
        slim_gen = ts.metadata["SLiM"]["generation"]

        rng = random.Random(seed)
        s1, s2 = rng.randrange(1, 2 ** 32), rng.randrange(1, 2 ** 32)

        population_configurations = [
            msprime.PopulationConfiguration(
                initial_size=pop.start_size, growth_rate=pop.growth_rate
            )
            for pop in recap_epoch.populations
        ]

        demography = msprime.Demography.from_old_style(
            population_configurations=population_configurations,
            migration_matrix=recap_epoch.migration_matrix,
        )

        # TODO: using recombination map in recapitate. maps should be discrete.
        ts = msprime.sim_ancestry(
            initial_state=ts,
            demography=demography,
            recombination_rate=contig.recombination_map.mean_rate,
            random_seed=s1,
        )

        ts = self._simplify_remembered(ts)

        # Figuring out SLiM mutation id metadata
        max_id = -1
        for mut in ts.mutations():
            for d in mut.derived_state.split(","):
                max_id = max(max_id, int(d))

        # Creating a SLiM mutation model
        model = msprime.SLiMMutationModel(
            type=len(contig.mutation_types), next_id=max_id + 1
        )

        # Add mutations to recapitated part of trees.
        ts = msprime.sim_mutations(
            ts,
            rate=contig.mutation_rate,
            model=model,
            start_time=slim_gen,
            keep=True,
            random_seed=s2,
        )

        # Adding neutral mutations to the SLiM part
        s3 = rng.randrange(1, 2 ** 32)
        ts = msprime.sim_mutations(
            ts,
            rate=msp_rate_map,
            model=model,
            end_time=slim_gen,
            keep=True,
            random_seed=s3,
        )

        return ts

    def recap_and_rescale(
        self,
        ts,
        demographic_model,
        contig,
        samples,
        extended_events=None,
        slim_scaling_factor=1.0,
        seed=None,
        **kwargs,
    ):
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
        :type ts: :class:`tskit.TreeSequence`

        .. warning::
            The :func:`recap_and_rescale` function is provided in the hope that
            it will be useful. But as we can't anticipate what changes you'll
            make to the SLiM code before using it, the stdpopsim source code
            should be consulted to determine if the behaviour is appropriate
            for your case.
        """
        msp_rate_map, slim_rate_map = get_msp_and_slim_mutation_rate_maps(contig)

        with open(os.devnull, "w") as script_file:
            recap_epoch = slim_makescript(
                script_file,
                "unused.trees",
                demographic_model,
                contig,
                samples,
                extended_events,
                slim_scaling_factor,
                1,
                slim_rate_map,
            )

        ts = self._recap_and_rescale(
            ts, seed, recap_epoch, contig, slim_scaling_factor, msp_rate_map
        )
        return ts


# SLiM does not currently work on Windows.
if sys.platform != "win32":
    stdpopsim.register_engine(_SLiMEngine())
