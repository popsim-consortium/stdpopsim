import allel
import operator
import pathlib
import numpy as np
import logging
import os
import glob
import tarfile
import stdpopsim

logger = logging.getLogger(__name__)
# make root directory for zarr annotations
annot_path = "annotations"
os.makedirs(annot_path, exist_ok=True)


def merged(intervals, *, closed: bool):
    """
    Merge overlapping and adjacent intervals.

    :param intervals: An iterable of (start, end) coordinates.
    :param bool closed: If True, [start, end] coordinates are closed,
        so [1, 2] and [3, 4] are adjacent intervals and will be merged.
        If False, [start, end) coordinates are half-open,
        so [1, 2) and [3, 4) are not adjacent and will not be merged.
    """

    def iter_merged(intervals, *, closed: bool):
        """
        Generate tuples of (start, end) coordinates for merged intervals.
        """
        intervals = sorted(intervals, key=operator.itemgetter(0))
        if len(intervals) == 0:
            return
        start, end = intervals[0]
        for a, b in intervals[1:]:
            assert a <= b
            if a > end + closed:
                # No intersection with the current interval.
                yield start, end
                start, end = a, b
            else:
                # Intersects, or is contiguous with, the current interval.
                end = max(end, b)
        yield start, end

    return list(iter_merged(intervals, closed=closed))


def gff_recarray_to_stdpopsim_intervals(gff):
    """
    Merge overlapping intervals and convert coordinates. GFF intervals are
    1-based [i,j], but stdpopsim intervals are 0-based [i-1,j).
    """
    intervals = np.array(merged(zip(gff.start, gff.end), closed=True))
    intervals[:, 0] = intervals[:, 0] - 1
    return intervals


def get_gff_recarray(url, sha256):
    local_path = pathlib.Path(url).name

    if not pathlib.Path(local_path).exists():
        logger.info(f"downloading {url}")
        stdpopsim.utils.download(url, local_path)

    logger.info("checking sha256")
    local_sha256 = stdpopsim.utils.sha256(local_path)
    if local_sha256 != sha256:
        logger.info(
            f"{local_path}: sha256: expected {sha256}, but found {local_sha256}. "
            "Delete the file to download it again."
        )
        exit(1)

    logger.info(f"loading {local_path} into numpy recarray")
    gff = allel.gff3_to_recarray(local_path)
    return gff


def make_tarfile(output_filename, source_dir, dest):
    if os.path.exists(output_filename):
        os.remove(output_filename)
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=dest)
    tar.close()


def download_process_annotations():
    """
    loop through all species and download annotation.
    from those annotations suck out what we want and
    store them in appropriately named files for upload
    """
    for spc in stdpopsim.all_species():
        if spc.annotations:
            for an in spc.annotations:
                CHROM_IDS = [chrom.id for chrom in spc.genome.chromosomes]
                logger.info(f"Downloading GFF file {an.id}")
                gff = get_gff_recarray(an.url, an.gff_sha256)
                logger.info(f"extracting annotations {an.id}")
                exons = gff[
                    np.where(
                        np.logical_and(
                            gff.source == an.annotation_source,
                            gff.type == an.annotation_type,
                        )
                    )
                ]
                logger.info(f"merging overlapping regions {an.id}")
                # create zarr store and zarr root
                spc_name_path = os.path.join(annot_path, spc.id)
                os.makedirs(spc_name_path, exist_ok=True)
                for chrom_id in CHROM_IDS:
                    chrom_exons = exons[np.where(exons.seqid == chrom_id)]
                    if len(chrom_exons) == 0:
                        continue
                    intervals = gff_recarray_to_stdpopsim_intervals(chrom_exons)
                    # double check that the intervals can be used in stdpopsim
                    stdpopsim.utils.check_intervals_validity(intervals)
                    out_file = os.path.join(
                        spc_name_path, an.file_pattern.format(id=chrom_id)
                    )
                    np.savetxt(out_file, intervals, fmt="%d")
                tf = spc_name_path + f"/{an.id}.tar.gz"
                make_tarfile(tf, spc_name_path, "")
                logger.info("made tarball at " + spc_name_path)
                for f in glob.glob(spc_name_path + "/*.txt"):
                    logger.info("removing " + f)
                    os.remove(f)
