"""
Code for liftOver between assemblies
"""

import sys
import os
import multiprocessing as mp
import subprocess as sp
import numpy as np
import matplotlib.pyplot as plt
from datetime import date
from matplotlib import gridspec


def assign_task(tmpFile, task_q, nProcs):
    """
    puts tasks into the queue
    """
    siteID = tmpFile
    c, i, nth_job = 0, 0, 1
    while (i + 1) * nProcs <= len(siteID):
        i += 1
    nP1 = nProcs - (len(siteID) % nProcs)
    for j in range(nP1):
        task_q.put((siteID[c : c + i], nth_job))
        nth_job += 1
        c = c + i
    for j in range(nProcs - nP1):
        task_q.put((siteID[c : c + i + 1], nth_job))
        nth_job += 1
        c = c + i + 1


def create_proc(nProcs, task_q, result_q, params):
    """
    starts a new process
    """
    for _ in range(nProcs):
        p = mp.Process(target=worker, args=(task_q, result_q, params))
        p.daemon = True
        p.start()


def hapmap_to_bed(inFile, outFile):
    """
    convert the hapmap-formatted file into a bed formatted file
    """
    chromName = None
    outFile += ".tmp1"
    oldRates = []
    weights = []
    oldRateDict = {}
    ct = 0
    try:
        with open(inFile, "r") as fIN, open(outFile, "w") as fOUT:
            fIN.readline()
            outline = []
            for line in fIN:
                ar = line.split()
                if len(outline) == 4:
                    outline.insert(2, ar[1])
                    oldRates.append([float(outline[1]), float(outline[3])])
                    weights.append(int(outline[2]) - int(outline[1]))
                    oldRateDict[outline[4]] = [
                        int(outline[1]),
                        int(outline[2]),
                        float(outline[3]),
                    ]
                    fOUT.write("\t".join([str(x) for x in outline]) + "\n")
                    outline = [x for x in ar[:3]]
                    outline.append("int%s" % (ct))
                    ct += 1
                else:
                    outline = [x for x in ar[:3]]
                    outline.append("int%s" % (ct))
                    ct += 1
        chromName = outline[0]
        old_average = np.average(
            np.take(np.array(oldRates), 1, 1), weights=np.array(weights)
        )
        return (
            outFile,
            chromName,
            np.array(oldRates),
            np.array(weights),
            oldRateDict,
            old_average,
        )
    except FileNotFoundError:
        return (
            outFile,
            chromName,
            np.array(oldRates),
            np.array(weights),
            oldRateDict,
            None,
        )


def lift_over(inFile, outFile, chainFile):
    """
    run the liftOver tool on the infile
    """
    newFile = outFile + ".tmp2"
    unmappedFile = outFile + ".tmp2.unmapped"
    cmd = [str(x) for x in ["liftOver", inFile, chainFile, newFile, unmappedFile]]
    sp.run(cmd)
    return newFile, unmappedFile


def base_overlap(target, query):
    """
    Return the number of bases that overlap between the target and query
    """
    return len(range(max(target[0], query[0]), min(target[1], query[1])))


def rescale_rate(oldRateDict, ar):
    """
    Return the rescaled rate from the new line
    """
    old_interval = oldRateDict[ar[4]]
    return (
        float(ar[3]) * (old_interval[1] - old_interval[0]) / (int(ar[2]) - int(ar[1]))
    )


def bed_to_hapmap(
    inFile, outFile, chromName, winLen, gapThresh, jobN, oldRateDict, adjAvg
):
    """
    Convert the lifted-over bed file to a hapmap-formatted
    file using weighted average
    """
    max_pos = 0
    rates = []
    intervalLens = []
    intervals = []
    with open(inFile, "r") as fIN:
        for line in fIN:
            ar = line.split()
            if ar[0] == chromName:
                rescaled_rate = rescale_rate(oldRateDict, ar)
                intervals.append([int(ar[1]), int(ar[2]), rescaled_rate, ar[0], ar[4]])
                rates.append(rescaled_rate)
                intervalLens.append(int(ar[2]) - int(ar[1]))
                max_pos = max(max_pos, int(ar[2]))
    rates = np.array(rates)
    intervalLens = np.array(intervalLens)
    sorted_intervals = sorted(intervals)
    new_chrom_average = np.average(rates, weights=intervalLens)

    # fill in gaps between intervals with the chromosome mean rate
    gap_ct, over_ct = 0, 0
    outLines = []
    skip_idx = []
    for i in range(len(sorted_intervals) - 1):
        if i in skip_idx:
            continue
        gap = sorted_intervals[i + 1][0] - sorted_intervals[i][1]
        if gap == 0:
            outLines.append(
                [
                    sorted_intervals[i][3],
                    sorted_intervals[i][0],
                    sorted_intervals[i][1],
                    sorted_intervals[i][2],
                ]
            )
        elif gapThresh >= gap > 0:
            gap_ct += 1
            if adjAvg:
                new_average = (
                    sorted_intervals[i][2] + sorted_intervals[i + 1][2]
                ) / 2.0
            else:
                new_average = new_chrom_average
            outLines.append(
                [
                    sorted_intervals[i][3],
                    sorted_intervals[i][0],
                    sorted_intervals[i][1],
                    sorted_intervals[i][2],
                ]
            )
            outLines.append(
                [
                    sorted_intervals[i][3],
                    sorted_intervals[i][1],
                    sorted_intervals[i + 1][0],
                    new_average,
                ]
            )
        elif gap > gapThresh:
            gap_ct += 1
            outLines.append(
                [
                    sorted_intervals[i][3],
                    sorted_intervals[i][0],
                    sorted_intervals[i][1],
                    sorted_intervals[i][2],
                ]
            )
            outLines.append(
                [
                    sorted_intervals[i][3],
                    sorted_intervals[i][1],
                    sorted_intervals[i + 1][0],
                    0.000000,
                ]
            )
        else:
            tmp_idx = [i]
            tmp_rate = [sorted_intervals[i][2]]
            tmp_weight = [sorted_intervals[i][1] - sorted_intervals[i][0]]
            tmp_max = sorted_intervals[i][1]
            j = 1
            while True:
                if sorted_intervals[i + j][0] < tmp_max:
                    tmp_max = max(tmp_max, sorted_intervals[i + j][1])
                    tmp_idx.append(i + j)
                    tmp_rate.append(sorted_intervals[i + j][2])
                    tmp_weight.append(
                        sorted_intervals[i + j][1] - sorted_intervals[i + j][0]
                    )
                    j += 1
                else:
                    break
            overlap_avg = np.average(np.array(tmp_rate), weights=np.array(tmp_weight))
            skip_idx.extend(tmp_idx)
            outLines.append(
                [sorted_intervals[i][3], sorted_intervals[i][0], tmp_max, overlap_avg]
            )

            # check for a gap after the overlapping windows
            gap = sorted_intervals[i + j][0] - tmp_max
            if gapThresh >= gap > 0:
                if adjAvg:
                    new_average = (sorted_intervals[i + j][2] + overlap_avg) / 2.0
                else:
                    new_average = new_chrom_average
                outLines.append(
                    [
                        sorted_intervals[i][3],
                        tmp_max,
                        sorted_intervals[i + j][0],
                        new_average,
                    ]
                )
            elif gap > gapThresh:
                gap_ct += 1
                outLines.append(
                    [
                        sorted_intervals[i][3],
                        tmp_max,
                        sorted_intervals[i + j][0],
                        0.000000,
                    ]
                )
            over_ct += 1
    # last interval
    if len(sorted_intervals) - 1 not in skip_idx:
        outLines.append(
            [
                sorted_intervals[-1][3],
                sorted_intervals[-1][0],
                sorted_intervals[-1][1],
                sorted_intervals[-1][2],
            ]
        )

    joinedRates = []
    with open(outFile, "w") as fOUT:
        fOUT.write(
            "\t".join([str(x) for x in ["Chromosome", "Position(bp)", "Rate(cM/Mb)"]])
            + "\n"
        )
        for line in outLines:
            fOUT.write(
                "\t".join([str(x) for x in [chromName, line[1], line[3]]]) + "\n"
            )
            joinedRates.append([float(x) for x in [line[1], line[3]]])
        fOUT.write("\t".join([str(x) for x in [chromName, max_pos, "0.000000"]]) + "\n")
    return np.array(joinedRates), new_chrom_average


def get_deviation(oldRates, weights, joinedRates, winLen, jobN, new_average):
    """
    calculate the deviation between oldRates and joinedRates
    using a weighted average
    """
    D = []
    max_pos = int(max(np.take(oldRates, 0, 1)[-1], np.take(joinedRates, 0, 1)[-1]))
    out_wins = []
    win = [0, winLen]
    next_start = [0, 0]
    if jobN == 0:
        print(
            "\nCalculating the deviation in rates between the "
            "original map and back-lifted map "
            "({}bp windows)...".format(winLen)
        )
    while win[1] < max_pos:
        out_wins.append(win[0])
        flag = [0, 0]
        avg_rates = []
        rateFiles = [oldRates, joinedRates]
        for i in range(len(rateFiles)):
            newRates = []
            weights = []
            for j in range(next_start[i], rateFiles[i].shape[0] - 1):
                if rateFiles[i][j][0] > win[1]:
                    break
                newRates.append(rateFiles[i][j][1])
                BO = base_overlap(
                    win, [int(rateFiles[i][j][0]), int(rateFiles[i][j + 1][0])]
                )
                if BO != 0 and flag[i] == 0:
                    next_start[i] = j
                    flag[i] = 1
                weights.append(BO)
            newRates = np.array(newRates)
            weights = np.array(weights)
            if np.sum(weights) == 0.0:
                weightedAve = 0.0
            else:
                weightedAve = np.average(newRates, weights=weights)
            avg_rates.append(weightedAve)
        D.append(avg_rates[1] - avg_rates[0])
        win[0] += winLen
        win[1] += winLen
    # calc last window
    out_wins.append(win[0])
    win[1] = max_pos
    avg_rates = []
    for i in range(len(rateFiles)):
        newRates = []
        weights = []
        for j in range(next_start[i], rateFiles[i].shape[0] - 1):
            newRates.append(rateFiles[i][j][1])
            weights.append(
                base_overlap(
                    win, [int(rateFiles[i][j][0]), int(rateFiles[i][j + 1][0])]
                )
            )
        newRates = np.array(newRates)
        weights = np.array(weights)
        if np.sum(weights) == 0.0:
            weightedAve = 0.0
        else:
            weightedAve = np.average(newRates, weights=weights)
        avg_rates.append(weightedAve)
    D.append(avg_rates[1] - avg_rates[0])

    return np.array(out_wins), np.array(D)


def create_readme(postmapDir, cli_cmd, gmID, chainFile, newAssemb):
    """
    creates a readme file to include with new build
    explaining that this script was used to create it
    and listing the CLI command used
    """
    cli = " ".join([str(x) for x in cli_cmd]) + "\n"
    readme = """\tThis directory contains a build of the {newAssemb} genetic map.\n
    The map was generated automatically by lifting over {gmID} to {newAssemb} \
            using `liftOver_assemblies.py` in stdpopsim.
    The chain file `{chainFile}` was downloaded from the UCSC Genome Browser \
            and the UCSC liftOver tool was used to perform the lift.\n
    The following command was used:\n
    {cli_cmd}
    {date}""".format(
        newAssemb=newAssemb,
        gmID=gmID,
        chainFile=os.path.basename(chainFile),
        cli_cmd=cli,
        date=date.today(),
    )

    readmeFile = os.path.join(
        postmapDir, "README_%s_map" % (os.path.basename(postmapDir))
    )
    with open(readmeFile, "w") as fOUT:
        fOUT.write(readme)


def create_readme_dm3TOdm6(postmapDir, cli_cmd, gmID, chainFile, newAssemb):
    """
    creates a readme file to include with new build
    explaining that this script was used to create it
    and listing the CLI command used
    """
    cli = " ".join([str(x) for x in cli_cmd]) + "\n"
    readme = """\tThis directory contains a build of the {newAssemb} genetic map.\n
    The map was generated automatically by lifting over {gmID}
    (ftp://ftp.flybase.org/flybase/associated_files/Comeron.2012.10.15.xlsx)
    to {newAssemb} using `liftOver_assemblies_dm3TOdm6.py` in stdpopsim.
    The chain file `{chainFile}` was downloaded from the UCSC Genome Browser,
    and the UCSC liftOver tool was used to perform the lift over.\n
    The following command was used:
    {cli_cmd}
    {date}""".format(
        newAssemb=newAssemb,
        gmID=gmID,
        chainFile=os.path.basename(chainFile),
        cli_cmd=cli,
        date=date.today(),
    )

    readmeFile = os.path.join(postmapDir, "README.txt")
    with open(readmeFile, "w") as fOUT:
        fOUT.write(readme)


def step_plot(
    validationDir, oldRates, joinedRates, chromName, old_average, new_average
):
    """
    plots the rates along the chromosomes for eyeball test
    """
    pdf_file = os.path.join(validationDir, chromName + ".pdf")
    plt.figure(figsize=(8, 4))
    plt.hlines(
        old_average,
        xmin=0.0,
        xmax=np.max(np.take(oldRates, 0, 1)),
        color="dimgrey",
        linestyle="dashed",
        label="original avg",
    )
    plt.hlines(
        new_average,
        xmin=0.0,
        xmax=np.max(np.take(joinedRates, 0, 1)),
        color="red",
        linestyle="dashed",
        label="lifted avg",
    )
    plt.step(
        np.take(oldRates, 0, 1),
        np.take(oldRates, 1, 1),
        where="post",
        linewidth=0.6,
        color="dimgrey",
        label="original",
    )
    plt.step(
        np.take(joinedRates, 0, 1),
        np.take(joinedRates, 1, 1),
        where="post",
        linewidth=0.6,
        color="red",
        alpha=0.65,
        label="lifted",
    )
    plt.legend(title=chromName)
    plt.xlabel("Chromosome position")
    plt.ylabel("Rate (cM/Mb)")
    plt.savefig(pdf_file)


def deviation_plot(validationDir, win_pos, D, chromName, new_average):
    """
    plots the deviation between rates for the original
    and the re-lifted original
    """
    pdf_file = os.path.join(validationDir, chromName + "_validation.pdf")
    plt.figure(figsize=(10, 4))
    gs = gridspec.GridSpec(1, 3, width_ratios=[4, 1, 1])
    ax0 = plt.subplot(gs[0])
    ax0.hlines(
        new_average,
        xmin=0.0,
        xmax=np.max(win_pos),
        color="red",
        linestyle="dashed",
        label="chrom avg",
    )
    ax0.hlines(
        new_average * -1.0,
        xmin=0.0,
        xmax=np.max(win_pos),
        color="red",
        linestyle="dashed",
    )
    ax0.step(
        win_pos, D, where="post", linewidth=0.6, color="dimgrey", label="deviation"
    )
    ax0.legend(title=chromName)
    ax0.set_xlabel("Chromosome position")
    ax0.set_ylabel("Rate deviation (cM/Mb)")
    ax1 = plt.subplot(gs[1])
    ax1.violinplot(D, showmeans=True, widths=0.95)
    ax2 = plt.subplot(gs[2])
    ax2.boxplot(D, notch=True, showfliers=False)
    plt.savefig(pdf_file)


def allDev_plot(dev, validationDir, genome, dm3TOdm6):
    """
    plots the deviation between rates for the
    original and the re-lifted original
    for all chromosomes as a single boxplot
    """
    rates, labels = [], []
    if dm3TOdm6:
        for chrom in genome:
            for d in dev:
                if d[0][3:] == chrom:
                    rates.append(d[1])
                    labels.append(d[0])
    else:
        for chrom in genome.chromosomes:
            for d in dev:
                if d[0] == chrom.id:
                    rates.append(d[1])
                    labels.append(d[0])
    pdf_file = os.path.join(validationDir, "all_chromosome_deviations.pdf")
    plt.figure(figsize=(10, 4))
    plt.xlabel("Chromosome")
    plt.ylabel("Rate deviation (cM/Mb)")
    plt.xticks(rotation=45)
    plt.boxplot(rates, labels=labels, notch=True, showfliers=False)
    plt.savefig(pdf_file)


def worker(task_q, result_q, params):
    """
    does the work
    """
    (
        intermediatesDir,
        validationDir,
        inFiles,
        outFiles,
        chainFile,
        winLen,
        gapThresh,
        validation,
        adjAvg,
    ) = params
    while True:
        try:
            index_range, nth_job = task_q.get()
            tmp = []
            for i in index_range:
                # convert hapmap to bed
                (
                    tmp1,
                    chromName,
                    oldRates,
                    weights,
                    oldRateDict,
                    old_average,
                ) = hapmap_to_bed(inFiles[i], outFiles[i])
                if chromName is None:
                    break
                # liftOver files
                tmp2, tmp2unmapped = lift_over(tmp1, outFiles[i], chainFile)
                # convert bed back to hapmap
                joinedRates, new_average = bed_to_hapmap(
                    tmp2,
                    outFiles[i],
                    chromName,
                    winLen,
                    gapThresh,
                    i,
                    oldRateDict,
                    adjAvg,
                )
                # move intermediates to new directory
                cmd = [
                    str(x) for x in ["mv", tmp1, tmp2, tmp2unmapped, intermediatesDir]
                ]
                sp.run(cmd)
                # plot the rates
                if validation:
                    (
                        tmp1,
                        chromName,
                        oldRates,
                        weights,
                        oldRateDict,
                        old_average,
                    ) = hapmap_to_bed(validation[i], outFiles[i])
                    win_pos, D = get_deviation(
                        oldRates, weights, joinedRates, winLen, i, new_average
                    )
                    deviation_plot(validationDir, win_pos, D, chromName, new_average)
                    tmp.append([chromName, D])
                else:
                    step_plot(
                        validationDir,
                        oldRates,
                        joinedRates,
                        chromName,
                        old_average,
                        new_average,
                    )
            if validation:
                result_q.put(tmp)
        finally:
            task_q.task_done()


def remap(
    intermediatesDir,
    validationDir,
    inFiles,
    outFiles,
    chainFile,
    winLen,
    gapThresh,
    adjAvg,
    nProc,
    validation=False,
    genome=False,
    dm3TOdm6=False,
):
    """
    run all functions to translate and liftOver file types
    """
    # multithread remapping
    task_q = mp.JoinableQueue()
    result_q = mp.Queue()
    params = (
        intermediatesDir,
        validationDir,
        inFiles,
        outFiles,
        chainFile,
        winLen,
        gapThresh,
        validation,
        adjAvg,
    )
    create_proc(nProc, task_q, result_q, params)
    assign_task(range(len(inFiles)), task_q, nProc)
    try:
        task_q.join()
    except KeyboardInterrupt:
        print("KeyboardInterrupt")
        sys.exit(0)
    else:
        if validation:
            all_deviations = []
            while nProc:
                x = result_q.get()
                for y in x:
                    all_deviations.append(y)
                nProc -= 1
            allDev_plot(all_deviations, validationDir, genome, dm3TOdm6)
        else:
            pass
        print("\n")


def split_chroms(chroms, mapFile, preliftDir):
    """
    split original Comeron file by chromosome and
    convert to Hapmap format
    """
    for chrom in chroms:
        splitFile = os.path.join(preliftDir, "{}.txt".format(chrom))
        with open(mapFile, "r") as fIN, open(splitFile, "w") as fOUT:
            fOUT.write(
                "\t".join(
                    [str(x) for x in ["Chromosome", "Position(bp)", "Rate(cM/Mb)"]]
                )
                + "\n"
            )
            fIN.readline()
            for line in fIN:
                ar = line.split()
                if ar[0] == chrom:
                    fOUT.write(
                        "\t".join(
                            [
                                str(x)
                                for x in [
                                    "chr{}".format(ar[0]),
                                    int(float(ar[2])) - 1,
                                    float(ar[4]),
                                ]
                            ]
                        )
                        + "\n"
                    )
            # add last line of HapMap file
            fOUT.write(
                "\t".join(
                    [str(x) for x in ["chr{}".format(chrom), chroms[chrom], 0.00]]
                )
                + "\n"
            )
