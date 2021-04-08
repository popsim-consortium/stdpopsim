"""
Code for liftOver between assemblies
tailored specifically for
the Comeron2012 map
(dm3 to dm6)
"""

import sys
import os
import argparse
import multiprocessing as mp
import pandas as pd
import tarfile
import urllib.request
import shutil
import liftOver


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--winLen",
        dest="winLen",
        help="calculate a weighted average over intervals of length winLen",
        type=int,
        default=100000,
    )
    parser.add_argument(
        "--gapThresh",
        dest="gapThresh",
        help="gaps exceeding this length will have rate zero "
        "instead of either chromosome/window average rate",
        type=int,
        default=1000000,
    )
    parser.add_argument(
        "--useChromosomeAvg",
        help="replaces the gap window with the chromosome average rate",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--useAdjacentAvg",
        help="replaces the gap window with the average rate "
        "of the two adjacent windows",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--retainIntermediates",
        help="retains a directory of the intermediate liftOver files",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--nProc",
        dest="nProc",
        help="number of cores to use (default uses all available)",
        type=int,
        default=0,
    )
    args = parser.parse_args()

    # Set number of cores to use
    if args.nProc == 0:
        nProc = mp.cpu_count()
    else:
        nProc = args.nProc

    # Save cli command
    cli_cmd = sys.argv

    # Create temporary directories
    orig_label = "dm3"
    newAssemb = "dm6"
    cwd = os.getcwd()
    tmpDir = os.path.join(cwd, "liftOver_tmp_%s" % (orig_label))
    intermediatesDir = os.path.join(cwd, "liftOver_intermediates_%s" % (orig_label))
    intermediatesDir_val = os.path.join(cwd, "liftOver_validation_intermediates")
    validationDir = os.path.join(cwd, "liftOver_validation")
    preliftDir = os.path.join(tmpDir, "prelift")
    postliftDir = os.path.join(tmpDir, "comeron2012v2_maps")
    postliftDir_val = os.path.join(
        tmpDir, newAssemb + "_liftedOverTo_{assemb}".format(assemb=orig_label)
    )
    for d in [
        tmpDir,
        intermediatesDir,
        validationDir,
        preliftDir,
        postliftDir,
        intermediatesDir_val,
        postliftDir_val,
    ]:
        if not os.path.exists(d):
            os.mkdir(d)

    # get liftover files
    chainUrl = (
        "http://hgdownload.soe.ucsc.edu/goldenPath/"
        "dm3/liftOver/dm3ToDm6.over.chain.gz"
    )
    validationUrl = (
        "http://hgdownload.soe.ucsc.edu/goldenPath/"
        "dm6/liftOver/dm6ToDm3.over.chain.gz"
    )
    chainFile = os.path.join(tmpDir, "dm3ToDm6.over.chain.gz")
    validationChain = os.path.join(tmpDir, "dm6ToDm3.over.chain.gz")
    urllib.request.urlretrieve(chainUrl, filename=chainFile)
    urllib.request.urlretrieve(validationUrl, filename=validationChain)
    mapUrl = "ftp://ftp.flybase.org/flybase/associated_files/Comeron.2012.10.15.xlsx"
    mapFile = os.path.join(tmpDir, "Comeron.2012.10.15.xlsx")
    urllib.request.urlretrieve(mapUrl, filename=mapFile)

    # convert xlsx map to plain text
    mapTxtFile = os.path.join(tmpDir, "Comeron.2012.10.15.txt")
    with open(mapTxtFile, "w") as fOUT:
        pd.read_excel(mapFile).to_string(fOUT, index=False)

    dm3_chrom_sizes = {
        "X": 22422827,
        "2L": 23011544,
        "2R": 21146708,
        "3L": 24543557,
        "3R": 27905053,
    }

    # Split original map into seperate chromosomes
    liftOver.split_chroms(dm3_chrom_sizes, mapTxtFile, preliftDir)

    # Assign all chromosomes to files
    inFiles, outFiles, valFiles = [], [], []
    for chrom in dm3_chrom_sizes:
        inFiles.append(os.path.join(preliftDir, "{}.txt".format(chrom)))
        outFiles.append(
            os.path.join(
                postliftDir,
                "genetic_map_comeron2012v2_{assemb}_chr{id}.txt".format(
                    assemb=newAssemb, id=chrom
                ),
            )
        )
        valFiles.append(
            os.path.join(
                postliftDir_val,
                "genetic_map_{assemb}_chr{id}.txt".format(assemb=orig_label, id=chrom),
            )
        )

    # Begin remapping
    print("\nMapping the rates to the new assembly...")
    liftOver.remap(
        intermediatesDir,
        validationDir,
        inFiles,
        outFiles,
        chainFile,
        args.winLen,
        args.gapThresh,
        args.useAdjacentAvg,
        nProc,
    )

    # Write new README
    liftOver.create_readme_dm3TOdm6(
        postliftDir, cli_cmd, orig_label, chainFile, newAssemb
    )

    # Reverse liftOver for validation
    print(
        "Remapping the liftedOver rates back to the "
        "original assembly for validation..."
    )
    liftOver.remap(
        intermediatesDir_val,
        validationDir,
        outFiles,
        valFiles,
        validationChain,
        args.winLen,
        args.gapThresh,
        args.useAdjacentAvg,
        nProc,
        validation=inFiles,
        genome=dm3_chrom_sizes,
        dm3TOdm6=True,
    )

    # Tar new files and remove temporary directory
    print("Creating new tarball...")
    outTarFile = os.path.join(cwd, "comeron2012v2_maps.tar.gz")
    with tarfile.open(outTarFile, "w:gz") as tb:
        tb.add(postliftDir, arcname=os.path.basename(postliftDir))
    shutil.rmtree(tmpDir)
    if not args.retainIntermediates:
        shutil.rmtree(intermediatesDir)
    print("\nliftOver_comeron2012.py complete!")


if __name__ == "__main__":
    main()
