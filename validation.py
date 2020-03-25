#!/usr/bin/env python3

import os
import sys
import pathlib
import argparse
import random
import functools
import itertools
import concurrent.futures
import time

import numpy as np
import tskit
import msprime
# Work around an issue on systems with large numbers of cores.
# https://github.com/cggh/scikit-allel/issues/285
os.environ["NUMEXPR_MAX_THREADS"] = f"{os.cpu_count()}"  # NOQA
import allel

import matplotlib
matplotlib.use('Agg')  # NOQA # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import stdpopsim
import stdpopsim.cli


def warning(msg):
    """
    Print a warning, with output less ugly than that of warnings.warn().
    """
    print(f"WARNING: {msg}", file=sys.stderr)


#
# Simulation functions.
#


def _onepop_PC(engine_id, out_dir, seed, N0=1000, *size_changes, **sim_kwargs):
    species = stdpopsim.get_species("CanFam")
    contig = species.get_contig("chr35", length_multiplier=0.01)  # ~265 kb
    model = stdpopsim.PiecewiseConstantSize(N0, *size_changes)
    model.generation_time = species.generation_time
    samples = model.get_samples(100)
    engine = stdpopsim.get_engine(engine_id)
    t0 = time.perf_counter()
    ts = engine.simulate(model, contig, samples, seed=seed, **sim_kwargs)
    t1 = time.perf_counter()
    out_file = out_dir / f"{seed}.trees"
    ts.dump(out_file)
    return out_file, t1 - t0


def onepop_constantN_msprime1(out_dir, seed):
    """
    Single population with constant population size.
    """
    return _onepop_PC("msprime", out_dir, seed)


def onepop_constantN_slim1(out_dir, seed):
    """
    Single population with constant population size.
    There are no demographic_events, so SLiM exits immediately, and
    tree sequences are constructed via recapitation.
    """
    return _onepop_PC("slim", out_dir, seed)


def onepop_constantN_slim2(out_dir, seed):
    """
    Single population with constant population size.
    No recapitation. No scaling.
    """
    return _onepop_PC(
            "slim", out_dir, seed,
            slim_no_recapitation=True, slim_scaling_factor=1)


def onepop_constantN_slim3(out_dir, seed):
    """
    Single population with constant population size.
    No recapitation. Time is rescaled by a factor of 10.
    """
    return _onepop_PC(
            "slim", out_dir, seed,
            slim_no_recapitation=True, slim_scaling_factor=10)


def onepop_bottleneck_msprime1(out_dir, seed):
    """
    Single population with bottleneck and recovery.
    """
    return _onepop_PC("msprime", out_dir, seed, 5000, (800, 100), (1000, 1000))


def onepop_bottleneck_slim1(out_dir, seed):
    """
    Single population with bottleneck and recovery.
    """
    return _onepop_PC("slim", out_dir, seed, 5000, (800, 100), (1000, 1000))


def onepop_bottleneck_slim2(out_dir, seed):
    """
    Single population with bottleneck and recovery.
    No recapitation. No scaling.
    """
    return _onepop_PC(
            "slim", out_dir, seed, 5000, (800, 100), (1000, 1000),
            slim_no_recapitation=True, slim_scaling_factor=1)


def onepop_bottleneck_slim3(out_dir, seed):
    """
    Single population with bottleneck and recovery.
    No recapitation. Time is rescaled by a factor of 10.
    """
    return _onepop_PC(
            "slim", out_dir, seed, 5000, (800, 100), (1000, 1000),
            slim_no_recapitation=True, slim_scaling_factor=10)


class _PiecewiseSize(stdpopsim.DemographicModel):
    """
    A copy of stdpopsim.PiecewiseConstantSize that permits growth rates.
    """
    id = "Piecewise"
    description = "Piecewise size population model over multiple epochs."
    citations = []
    populations = [stdpopsim.Population(id="pop0", description="Population 0")]
    author = None
    year = None
    doi = None

    def __init__(self, N0, growth_rate, *args):
        self.population_configurations = [
            msprime.PopulationConfiguration(
                initial_size=N0, growth_rate=growth_rate,
                metadata=self.populations[0].asdict())
        ]
        self.migration_matrix = [[0]]
        self.demographic_events = []
        for t, initial_size, growth_rate in args:
            self.demographic_events.append(msprime.PopulationParametersChange(
                time=t, initial_size=initial_size, growth_rate=growth_rate,
                population_id=0))


def _onepop_expgrowth(
        engine_id, out_dir, seed, N0=5000, N1=500, T=1000, **sim_kwargs):
    growth_rate = - np.log(N1 / N0) / T
    species = stdpopsim.get_species("DroMel")
    contig = species.get_contig("chr2R", length_multiplier=0.01)  # ~250 kb
    model = _PiecewiseSize(N0, growth_rate, (T, N1, 0))
    model.generation_time = species.generation_time
    samples = model.get_samples(100)
    engine = stdpopsim.get_engine(engine_id)
    t0 = time.perf_counter()
    ts = engine.simulate(model, contig, samples, seed=seed, **sim_kwargs)
    t1 = time.perf_counter()
    out_file = out_dir / f"{seed}.trees"
    ts.dump(out_file)
    return out_file, t1 - t0


def onepop_expgrowth_msprime1(out_dir, seed):
    """
    Single population with exponential population size growth.
    """
    return _onepop_expgrowth("msprime", out_dir, seed)


def onepop_expgrowth_slim1(out_dir, seed):
    """
    Single population with exponential population size growth.
    There are no demographic_events, so SLiM exits immediately, and
    tree sequences are constructed via recapitation.
    """
    return _onepop_expgrowth("slim", out_dir, seed)


def onepop_expgrowth_slim2(out_dir, seed):
    """
    Single population with exponential population size growth.
    No recapitation. No scaling.
    """
    return _onepop_expgrowth(
            "slim", out_dir, seed,
            slim_no_recapitation=True, slim_scaling_factor=1)


def onepop_expgrowth_slim3(out_dir, seed):
    """
    Single population with exponential population size growth.
    No recapitation. Time is rescaled by a factor of 10.
    """
    return _onepop_expgrowth(
            "slim", out_dir, seed,
            slim_no_recapitation=True, slim_scaling_factor=10)


def _twopop_IM(
        engine_id, out_dir, seed,
        NA=1000, N1=500, N2=5000, T=1000, M12=0, M21=0, pulse=None,
        **sim_kwargs):
    species = stdpopsim.get_species("AraTha")
    contig = species.get_contig("chr5", length_multiplier=0.01)  # ~270 kb
    model = stdpopsim.IsolationWithMigration(
            NA=NA, N1=N1, N2=N2, T=T, M12=M12, M21=M21)
    if pulse is not None:
        model.demographic_events.append(pulse)
        model.demographic_events.sort(key=lambda x: x.time)
    model.generation_time = species.generation_time
    samples = model.get_samples(50, 50, 0)
    engine = stdpopsim.get_engine(engine_id)
    t0 = time.perf_counter()
    ts = engine.simulate(model, contig, samples, seed=seed, **sim_kwargs)
    t1 = time.perf_counter()
    out_file = out_dir / f"{seed}.trees"
    ts.dump(out_file)
    return out_file, t1 - t0


def twopop_no_migration_msprime1(out_dir, seed):
    """
    Two populations with different sizes and no migrations.
    """
    return _twopop_IM("msprime", out_dir, seed)


def twopop_no_migration_slim1(out_dir, seed):
    """
    Two populations with different sizes and no migrations.
    Recapitation. Default scaling factor of 10 is applied.
    """
    return _twopop_IM("slim", out_dir, seed)


def twopop_no_migration_slim2(out_dir, seed):
    """
    Two populations with different sizes and no migrations.
    No recapitation. No scaling.
    """
    return _twopop_IM(
            "slim", out_dir, seed,
            slim_no_recapitation=True, slim_scaling_factor=1)


def twopop_asymmetric_migration_msprime1(out_dir, seed):
    """
    Two populations with different sizes and migrations from pop2 to pop1.
    """
    return _twopop_IM("msprime", out_dir, seed, M12=0, M21=0.001)


def twopop_asymmetric_migration_slim1(out_dir, seed):
    """
    Two populations with different sizes and migrations from pop2 to pop1.
    Recapitation. Default scaling factor of 10 is applied.
    """
    return _twopop_IM("slim", out_dir, seed, M12=0, M21=0.001)


def twopop_asymmetric_migration_slim2(out_dir, seed):
    """
    Two populations with different sizes and migrations from pop2 to pop1.
    No recapitation. No scaling.
    """
    return _twopop_IM(
            "slim", out_dir, seed, M12=0, M21=0.001,
            slim_no_recapitation=True, slim_scaling_factor=1)


_pulse_m21 = msprime.MassMigration(
        time=20, proportion=0.1, source=1, destination=0)


def twopop_pulse_migration_msprime1(out_dir, seed):
    """
    Two populations with different sizes and introgression from pop2 to pop1.
    """
    return _twopop_IM("msprime", out_dir, seed, pulse=_pulse_m21)


def twopop_pulse_migration_slim1(out_dir, seed):
    """
    Two populations with different sizes and introgression from pop2 to pop1.
    Recapitation. Default scaling factor of 10 is applied.
    """
    return _twopop_IM("slim", out_dir, seed, pulse=_pulse_m21)


def twopop_pulse_migration_slim2(out_dir, seed):
    """
    Two populations with different sizes and introgression from pop2 to pop1.
    No recapitation. No scaling.
    """
    return _twopop_IM(
            "slim", out_dir, seed, pulse=_pulse_m21,
            slim_no_recapitation=True, slim_scaling_factor=1)


def do_cmd(cmd, out_dir, seed):
    cmd = cmd.split()
    assert "-o" not in cmd and "--output" not in cmd
    assert "-s" not in cmd and "--seed" not in cmd
    out_file = out_dir / f"{seed}.trees"
    full_cmd = cmd + f" -q -o {out_file} -s {seed}".split()
    t0 = time.perf_counter()
    stdpopsim.cli.stdpopsim_main(full_cmd)
    t1 = time.perf_counter()
    assert os.path.exists(out_file)
    return out_file, t1 - t0


_homsap_250k = " HomSap -c chr1 -l 0.001 "


def Africa_1T12_msprime1(out_dir, seed):
    cmd = "-e msprime" + _homsap_250k + "-d Africa_1T12 100"
    return do_cmd(cmd, out_dir, seed)


def Africa_1T12_slim1(out_dir, seed):
    cmd = "-e slim" + _homsap_250k + "-d Africa_1T12 100"
    return do_cmd(cmd, out_dir, seed)


def OutOfAfrica_3G09_msprime1(out_dir, seed):
    samples = 3 * " 33"
    cmd = "-e msprime" + _homsap_250k + "-d OutOfAfrica_3G09" + samples
    return do_cmd(cmd, out_dir, seed)


def OutOfAfrica_3G09_slim1(out_dir, seed):
    samples = 3 * " 33"
    cmd = "-e slim" + _homsap_250k + "-d OutOfAfrica_3G09" + samples
    return do_cmd(cmd, out_dir, seed)


def AmericanAdmixture_4B11_msprime1(out_dir, seed):
    samples = 4 * " 25"
    cmd = "-e msprime" + _homsap_250k + "-d AmericanAdmixture_4B11" + samples
    return do_cmd(cmd, out_dir, seed)


def AmericanAdmixture_4B11_slim1(out_dir, seed):
    samples = 4 * " 25"
    cmd = "-e slim" + _homsap_250k + "-d AmericanAdmixture_4B11" + samples
    return do_cmd(cmd, out_dir, seed)


def AncientEurasia_9K19_msprime1(out_dir, seed):
    samples = 8 * " 12"
    cmd = "-e msprime" + _homsap_250k + "-d AncientEurasia_9K19" + samples
    return do_cmd(cmd, out_dir, seed)


def AncientEurasia_9K19_slim1(out_dir, seed):
    samples = 8 * " 12"
    cmd = "-e slim" + _homsap_250k + "-d AncientEurasia_9K19" + samples
    return do_cmd(cmd, out_dir, seed)


#
# Stats functions.
#


def tmrca(ts):
    """
    Time to most recent common ancestor of sample, aka tree height.
    """
    tmrcas = [tree.time(tree.root) for tree in ts.trees()]
    min_, median, max_ = np.quantile(tmrcas, (0, 0.5, 1))
    return {"min(tmrca)": min_,
            "median(tmrca)": median,
            "max(tmrca)": max_,
            }


def ts_properties(ts):
    """
    TreeSequence properties.
    """
    return {"num_trees": ts.num_trees,
            "num_edges": ts.num_edges,
            "num_nodes": ts.num_nodes,
            "num_sites": ts.num_sites,
            }


def pooled_pop_stats(ts):
    """
    Population statistics, with samples pooled from all populations.
    """
    n = ts.num_samples // 2
    sample_sets = [ts.samples()[:n], ts.samples()[n:]]
    return {"diversity": ts.diversity(),
            "Tajimas_D": ts.Tajimas_D(),
            "$f_2$": ts.f2(sample_sets),
            "$Y_2$": ts.Y2(sample_sets),
            "segregating_sites": ts.segregating_sites(),
            }


def pairwise_pop_stats(ts):
    """
    Pairwise population statistics, calculated for all pairs of populations.
    """
    pops = [i for i in range(ts.num_populations) if len(ts.samples(i)) > 0]
    if len(pops) < 2:
        return None
    stats = dict()
    for j, k in itertools.combinations(pops, 2):
        sample_sets = [ts.samples(j), ts.samples(k)]
        stats[f"$f_2$[{j},{k}]"] = ts.f2(sample_sets)
        stats[f"$Y_2$[{j},{k}]"] = ts.Y2(sample_sets)
    return stats


def linkage_disequilibrium(ts, span=2*10**5, bins=50, min_obs_per_bin=8):
    """
    Average R^2 in `bins` bins over the first `span` bases of ts.
    """
    span = min(ts.sequence_length, span)
    ts = ts.keep_intervals([(0, span)], record_provenance=False)
    position = [site.position for site in ts.sites()]
    num_sites = len(position)
    assert num_sites == int(ts.num_sites)

    nans = np.full(bins, np.nan)
    if num_sites >= min_obs_per_bin:
        gts = np.expand_dims(ts.genotype_matrix(), axis=-1)
        gn = allel.GenotypeArray(gts, dtype='i1').to_n_alt()
        ld = allel.rogers_huff_r(gn)**2
        assert len(ld) == num_sites * (num_sites - 1) // 2

        # Bin the pairwise site R^2 in `ld` by site separation distance.
        r2 = np.zeros(bins)
        n = np.zeros(bins)
        i = 0
        for j in range(num_sites):
            for k in range(j+1, num_sites):
                distance = position[k] - position[j]
                index = int(distance * bins / span)
                if not np.isnan(ld[i]):
                    r2[index] += ld[i]
                    n[index] += 1
                i += 1
        # divide `r2` by `n`, but return NaN where n has insufficient observations.
        r2 = np.divide(r2, n, out=nans, where=n >= min_obs_per_bin)
    else:
        # Too few segregating sites to do anything meaningful.
        # LD plots may be blank.
        r2 = nans

    a = f"{span//bins//1000}k"  # width of one bin, in kb
    b = f"{span//8//1000}k"
    c = f"{span//4//1000}k"
    d = f"{span//2//1000}k"

    return {f"$R^2$[<{a}]": r2[0],
            f"$R^2$[{b}]": r2[bins//8],
            f"$R^2$[{c}]": r2[bins//4],
            f"$R^2$[{d}]": r2[bins//2]
            }


_simulation_functions = [
    onepop_constantN_msprime1,
    onepop_constantN_slim1,
    onepop_constantN_slim2,
    onepop_constantN_slim3,
    onepop_bottleneck_msprime1,
    onepop_bottleneck_slim1,
    onepop_bottleneck_slim2,
    onepop_bottleneck_slim3,
    onepop_expgrowth_msprime1,
    onepop_expgrowth_slim1,
    onepop_expgrowth_slim2,
    onepop_expgrowth_slim3,

    twopop_no_migration_msprime1,
    twopop_no_migration_slim1,
    twopop_no_migration_slim2,
    twopop_asymmetric_migration_msprime1,
    twopop_asymmetric_migration_slim1,
    twopop_asymmetric_migration_slim2,
    twopop_pulse_migration_msprime1,
    twopop_pulse_migration_slim1,
    twopop_pulse_migration_slim2,

    Africa_1T12_msprime1,
    Africa_1T12_slim1,
    OutOfAfrica_3G09_msprime1,
    OutOfAfrica_3G09_slim1,
    AmericanAdmixture_4B11_msprime1,
    AmericanAdmixture_4B11_slim1,
    AncientEurasia_9K19_msprime1,
    AncientEurasia_9K19_slim1,
]

_stats_functions = [
    ts_properties,
    tmrca,
    pooled_pop_stats,
    pairwise_pop_stats,
    linkage_disequilibrium,
]

_default_comparisons = [
    (onepop_constantN_msprime1, onepop_constantN_slim1),
    (onepop_constantN_msprime1, onepop_constantN_slim2),
    (onepop_constantN_msprime1, onepop_constantN_slim3),
    (onepop_bottleneck_msprime1, onepop_bottleneck_slim1),
    (onepop_bottleneck_msprime1, onepop_bottleneck_slim2),
    (onepop_expgrowth_msprime1, onepop_expgrowth_slim1),
    (onepop_expgrowth_msprime1, onepop_expgrowth_slim2),
    (onepop_expgrowth_msprime1, onepop_expgrowth_slim3),

    (twopop_no_migration_msprime1, twopop_no_migration_slim1),
    (twopop_no_migration_msprime1, twopop_no_migration_slim2),
    (twopop_asymmetric_migration_msprime1, twopop_asymmetric_migration_slim1),
    (twopop_asymmetric_migration_msprime1, twopop_asymmetric_migration_slim2),
    (twopop_pulse_migration_msprime1, twopop_pulse_migration_slim1),
    (twopop_pulse_migration_msprime1, twopop_pulse_migration_slim2),

    (Africa_1T12_msprime1, Africa_1T12_slim1),
    (OutOfAfrica_3G09_msprime1, OutOfAfrica_3G09_slim1),
    (AmericanAdmixture_4B11_msprime1, AmericanAdmixture_4B11_slim1),
    (AncientEurasia_9K19_msprime1, AncientEurasia_9K19_slim1),
]

stats_functions = {f.__name__: f for f in _stats_functions}
simulation_functions = {f.__name__: f for f in _simulation_functions}
default_comparisons = [(t[0].__name__, t[1].__name__)
                       for t in _default_comparisons]


def do_simulations(rng, path, num_replicates, executor, key):
    out_dir = path / "trees" / key
    out_dir.mkdir(parents=True, exist_ok=True)
    func = functools.partial(simulation_functions[key], out_dir)
    seeds = (rng.randrange(1, 2**32) for _ in range(num_replicates))
    res = list(executor.map(func, seeds))
    files, times = zip(*res)
    # dump timing info to a file
    np.savetxt(out_dir / "times.txt", times)
    return files, times


def find_simulations(path, key):
    out_dir = path / "trees" / key
    files = list(out_dir.glob("*.trees"))
    if len(files) == 0:
        raise RuntimeError(f"{out_dir}: no *.trees found.")

    times_file = out_dir / "times.txt"
    if times_file.exists():
        times = np.loadtxt(times_file)
    else:
        warning(f"No times.txt found for {key}")

    return files, times


def compute_stats(ts_file):
    st = dict()
    ts = tskit.load(ts_file)
    for key, func in stats_functions.items():
        try:
            res = func(ts)
        except Exception:
            # Print the filename so it's easier to trace problems.
            warning(f"{ts_file} triggered exception")
            raise
        if res is not None:
            st[key] = res
    return st


def do_plots(path, sim_key1, sim_key2, times, stats):
    plotdir = path / "plots"
    plotdir.mkdir(parents=True, exist_ok=True)

    cmap = plt.get_cmap("tab10")
    markers = "oXdPvp*"
    scale = 1.25
    fig_w, fig_h = plt.figaspect(9.0/16.0)
    figsize = (scale*fig_w, scale*fig_h)

    times1, times2 = times[sim_key1], times[sim_key2]
    stats1, stats2 = stats[sim_key1], stats[sim_key2]

    pdf = PdfPages(plotdir / f"{sim_key1}__{sim_key2}.pdf")

    # plot run times
    assert len(times1) > 0 and len(times2) > 0
    fig, ax = plt.subplots(figsize=figsize)
    ax.violinplot([times1, times2])
    ax.set_xticklabels([sim_key1, sim_key2])
    ax.set_title(f"{sim_key1} vs. {sim_key2}: run time")
    ax.set_ylabel("time (seconds)")
    fig.tight_layout()
    pdf.savefig(figure=fig)

    # plot stats
    quantiles = np.linspace(0, 1, 100)
    for stat_key in stats_functions.keys():
        if stat_key not in stats1[0]:
            continue
        inner_keys = stats1[0][stat_key].keys()
        assert inner_keys == stats2[0][stat_key].keys()
        nrows = int(np.ceil(np.sqrt(len(inner_keys))))
        ncols = int(np.ceil(len(inner_keys) / nrows))
        fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
        axs = np.array(axs).reshape(-1)
        assert len(axs) >= len(inner_keys)
        imarkers = itertools.cycle(markers)
        icolour = itertools.cycle(cmap.colors)
        save_fig = False
        for ax, inner_key in zip(axs, inner_keys):
            x = [d[stat_key][inner_key] for d in stats1]
            y = [d[stat_key][inner_key] for d in stats2]
            assert len(x) > 0 and len(y) > 0
            if np.all(np.isnan(x)) or np.all(np.isnan(y)):
                continue

            xq = np.nanquantile(x, quantiles)
            yq = np.nanquantile(y, quantiles)

            ax.scatter(xq, yq, ec=next(icolour), fc="none", marker=next(imarkers))
            ax.set_title(inner_key)

            # draw a diagonal line
            min_ = min(np.min(xq), np.min(yq))
            max_ = max(np.max(xq), np.max(yq))
            ax.plot([min_, max_], [min_, max_],
                    c="lightgray", ls="--", lw=1, zorder=-10)

            save_fig = True

        if save_fig:
            # use a full-figure subplot for labels that span the other subplots
            ax = fig.add_subplot(111, frameon=False)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(f"{sim_key1} vs. {sim_key2}: {stat_key}", pad=30)
            ax.set_xlabel(sim_key1, labelpad=30)
            ax.set_ylabel(sim_key2, labelpad=50)
            fig.tight_layout()
            pdf.savefig(figure=fig, bbox_inches='tight')
        plt.close(fig)

    pdf.close()


def parse_args():
    parser = argparse.ArgumentParser(
            description="Do validation simulations and make QQ plots.")
    parser.add_argument(
            "-o", "--output-folder", metavar="DIR",
            type=pathlib.Path, default=pathlib.Path("validation"),
            help="Folder to store validation plots and tree sequences "
                 "[%(default)s].")

    mutex_group = parser.add_mutually_exclusive_group()
    mutex_group.add_argument(
            "-n", "--no-plots", action="store_true", default=False,
            help="Don't make plots, just do the simulations [%(default)s].")
    mutex_group.add_argument(
            "-p", "--plot-only", action="store_true", default=False,
            help="Don't simulate, just make QQ plots from preexisting files "
                 "[%(default)s].")

    parser.add_argument(
            "-j", "--num-procs", metavar="NPROCS", type=int, default=1,
            help="Number of simulations to run simultaneously [%(default)s].")
    parser.add_argument(
            "-r", "--num-replicates", metavar="NREPS", type=int, default=100,
            help="Number of replicates for each simulation key [%(default)s].")
    parser.add_argument(
            "-s", "--seed", metavar="SEED", type=int, default=1234,
            help="Seed for the random number generator [%(default)s].")
    parser.add_argument(
            "keys", nargs="*",
            help="One or more scenarios to simulate and/or compare.")

    args = parser.parse_args()

    if len(args.keys) == 0:
        args.comparisons = default_comparisons
        args.keys = list(set(itertools.chain(*args.comparisons)))
    else:
        args.comparisons = itertools.combinations(args.keys, 2)

    # sort keys to get deterministic ordering from random number generator
    args.keys.sort()

    for key in args.keys:
        if key not in simulation_functions:
            if args.plot_only:
                # Might be a mistake, but continue anyway to allow validation
                # using arbitrary folders that are in the right place.
                warning(f"unknown scenario key ``{key}''")
            else:
                parser.error(f"unknown scenario key ``{key}''")

    return args


if __name__ == "__main__":
    args = parse_args()
    rng = random.Random(args.seed)

    files = dict()
    times = dict()
    stats = dict()

    with concurrent.futures.ProcessPoolExecutor(args.num_procs) as executor:
        for sim_keys in args.comparisons:
            j, k = sim_keys
            assert j != k
            print(f"{j} / {k}.", end="")

            for key in sim_keys:
                if key in files:
                    assert key in times
                    assert key in stats
                    continue

                if not args.plot_only:
                    files[key], times[key] = do_simulations(
                            rng, args.output_folder, args.num_replicates,
                            executor, key)
                else:
                    files[key], times[key] = find_simulations(
                            args.output_folder, key)

                print(".", end="")

                if not args.no_plots:
                    stats[key] = list(executor.map(compute_stats, files[key]))
                    print(".", end="")

            if not args.no_plots:
                do_plots(args.output_folder, j, k, times, stats)

            print("done")
