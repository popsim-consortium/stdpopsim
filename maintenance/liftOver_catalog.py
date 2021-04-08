"""
Code for liftOver between assemblies
from the species catalog
"""

import sys
import os
import argparse
import multiprocessing as mp
import tarfile
import urllib.request
import shutil
import stdpopsim
import liftOver


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", dest="species", help="stdpopsim species name")
    parser.add_argument("--map", dest="map", help="stdpopsim map name")
    parser.add_argument(
        "--chainFile", dest="chainFile", help="liftOver chain file from UCSC"
    )
    parser.add_argument(
        "--validationChain",
        dest="validationChain",
        help="liftOver chain file from UCSC",
    )
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
        help="gaps exceeding this length will have rate zero",
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
        help="replaces the gap window with the average rate of the "
        "two adjacent windows",
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

    # Grab metadata from stdpopsim
    species = stdpopsim.get_species(args.species)
    genome = species.genome
    gm = species.get_genetic_map(args.map)

    # Create temporary directories
    cwd = os.getcwd()
    tmpDir = os.path.join(cwd, "liftOver_tmp_%s" % (args.map))
    intermediatesDir = os.path.join(cwd, "liftOver_intermediates_%s" % (args.map))
    intermediatesDir_val = os.path.join(cwd, "liftOver_validation_intermediates")
    validationDir = os.path.join(cwd, "liftOver_validation")
    preliftDir = os.path.join(tmpDir, "prelift")
    newAssemb = os.path.basename(args.chainFile).split("To")[-1].split(".over")[0]
    postliftDir = os.path.join(
        tmpDir, gm.id + "_liftedOverTo_{assemb}".format(assemb=newAssemb)
    )
    postliftDir_val = os.path.join(
        tmpDir, newAssemb + "_liftedOverTo_{assemb}".format(assemb=gm.id)
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

    # Download map and extract all files from tarball
    tarFile = os.path.join(preliftDir, "tb.tar")
    print("\nDownloading tarball...")
    urllib.request.urlretrieve(gm.url, filename=tarFile)
    print("Extracting tarball...")
    tb = tarfile.open(tarFile)
    tb.extractall(preliftDir)
    tb.close

    # Assign all chromosomes to files
    newAssemb = os.path.basename(args.chainFile).split("To")[-1].split(".over")[0]
    inFiles, outFiles, valFiles = [], [], []
    for chrom in genome.chromosomes:
        inFiles.append(os.path.join(preliftDir, gm.file_pattern.format(id=chrom.id)))
        outFiles.append(
            os.path.join(
                postliftDir,
                "genetic_map_{assemb}_{id}.txt".format(assemb=newAssemb, id=chrom.id),
            )
        )
        valFiles.append(
            os.path.join(
                postliftDir_val,
                "genetic_map_{assemb}_{id}.txt".format(assemb=gm.id, id=chrom.id),
            )
        )

    # Begin remapping
    print("\nMapping the rates to the new assembly...")
    liftOver.remap(
        intermediatesDir,
        validationDir,
        inFiles,
        outFiles,
        args.chainFile,
        args.winLen,
        args.gapThresh,
        args.useAdjacentAvg,
        nProc,
    )

    # Write new README
    liftOver.create_readme(postliftDir, cli_cmd, gm.id, args.chainFile, newAssemb)

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
        args.validationChain,
        args.winLen,
        args.gapThresh,
        args.useAdjacentAvg,
        nProc,
        validation=inFiles,
        genome=genome,
    )

    # Tar new files and remove temporary directory
    print("Creating new tarball...")
    outTarFile = os.path.join(
        cwd, gm.id + "_liftedOverTo_{assemb}.tar.gz".format(assemb=newAssemb)
    )
    with tarfile.open(outTarFile, "w:gz") as tb:
        tb.add(postliftDir, arcname=os.path.basename(postliftDir))
    shutil.rmtree(tmpDir)
    if not args.retainIntermediates:
        shutil.rmtree(intermediatesDir)
    print("\nliftOver_catalog.py complete!")


if __name__ == "__main__":
    main()
