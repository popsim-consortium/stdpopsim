"""
Miscellaneous utilities.
"""
import re
import os
import hashlib
import urllib.request
import shutil
import tarfile
import contextlib
import numpy as np


def is_valid_demographic_model_id(model_id):
    """
    Returns True if the specified string is a valid demographic model ID. This must
    be a string with the following pattern:

    {CamelCaseName}_{num populations}{First letter author name}{2 digit year}.
    """
    regex = re.compile(r"[A-Z][A-Za-z0-9]*_[1-9]\d*[A-Z]\d\d")
    return regex.fullmatch(model_id) is not None


def is_valid_dfe_id(model_id):
    """
    Returns True if the specified string is a valid dfe model ID. This must
    be a string with the following pattern:

    {CamelCaseName}_{First author initial}{2 digit year}.
    """
    regex = re.compile(r"[A-Z][A-Za-z]*_[A-Z]*\d\d")
    return regex.fullmatch(model_id) is not None


def is_valid_genetic_map_id(gmap_id):
    """
    Returns True if the specified string is a valid genetic map ID. This must
    be a string with the following pattern:

    {CamelCaseName}_{assembly ID}
    """
    regex = re.compile(r"[A-Z][A-Za-z0-9]*_\w+")
    return regex.fullmatch(gmap_id) is not None


def is_valid_species_id(species_id):
    """
    Returns True if the specified string is a valid species ID. This must
    be a 6 letter CamelCase identifier.
    """
    regex = re.compile(r"[A-Z][a-z]{,2}[A-Z][a-z]{,2}")
    return regex.fullmatch(species_id) is not None


def is_valid_species_name(name):
    """
    Returns True if the specified string is a valid species name. This
    must be two or more words with first letter capitalised.
    """
    # FIXME this only supports two words for now. See #329
    regex = re.compile(r"[A-Z][a-z]+ [a-z]+")
    return regex.fullmatch(name) is not None


def is_valid_species_common_name(common_name):
    """
    Returns True if the specified string is a valid species common name. This
    must start with a capital letter.
    """
    # FIXME any sensible restrictions we can make on common names? See #330.
    regex = re.compile(r"[A-Z].*")
    return regex.fullmatch(common_name) is not None


def download(url, filename):
    """
    Download url to the specified local file.
    """
    # TODO: what is a sensible timeout here?
    with urllib.request.urlopen(url, timeout=30) as f_in:
        with open(filename, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def sha256(filename):
    """
    Return the SHA256 hex digest for the specified file.
    """
    m = hashlib.sha256()
    BUFLEN = 4096 * m.block_size  # 256 Kib
    with open(filename, "rb") as f:
        while True:
            buf = f.read(BUFLEN)
            if len(buf) == 0:
                break
            m.update(buf)
    return m.hexdigest()


@contextlib.contextmanager
def cd(path):
    """
    Convenience function to change the current working directory in a context manager.
    """
    old_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old_dir)


def untar(filename, path):
    """
    Extract the optionally-gzipped tar file to the specifed path.
    """
    with tarfile.open(filename, "r") as tf:
        for info in tf.getmembers():
            # Due to security concerns, we only extract tarballs containing a
            # very restrictive set of file types. See the warning here:
            # https://docs.python.org/3/library/tarfile.html#tarfile.TarFile.extractall
            if not (info.isfile() or info.isdir()):
                raise ValueError(f"Tarball format error: member {info.name} not a file")
            if info.name.startswith("/") or info.name.startswith(".."):
                raise ValueError(f"Refusing to extract {info.name} outside of {path}")
        with cd(path):
            tf.extractall()


def read_bed(mask_fpath, chrom):
    """
    Returns intervals to keep based on a bed file specified by the mask_fpath.
    The mask must be in bed format (columns specify chrom, left, right) and
    additional columns are ignored. Intervals must be non-overlapping.

    Note that the chromosome name must match exactly (i.e. "22" is not equivalent
    to "chr22").
    """
    lines = np.loadtxt(
        mask_fpath,
        dtype={"names": ("chrom", "left", "right"), "formats": (object, int, int)},
        delimiter="\t",
        usecols=(0, 1, 2),
    )
    in_chrom = lines["chrom"] == f"{chrom}"
    lines = lines.compress(in_chrom)
    intervals = np.array([lines["left"], lines["right"]]).T
    return intervals


def mask_tree_sequence(ts, mask_intervals, exclude):
    """
    Return a masked tree sequence, based on the mask intervals and whether
    they are kept or excluded.
    """
    if exclude is False:
        ts = ts.keep_intervals(mask_intervals)
    else:
        ts = ts.delete_intervals(mask_intervals)
    return ts


def append_common_synonyms(genome):
    """
    For common chromosome IDs, add their obvious synonyms if they do not exist already.
    The input is a stdpopsim.Genome object, and synonyms are appended in place.
    """

    def add_if_unique(chrom, syn):
        if syn not in chrom.synonyms:
            chrom.synonyms.append(syn)

    for chrom in genome.chromosomes:
        # "1" -> "chr1" and "2L" -> "chr2L"
        if chrom.id[0].isdigit():
            add_if_unique(chrom, "chr" + chrom.id)
        # sex chroms
        if chrom.id in ["X", "Y", "W", "Z"]:
            add_if_unique(chrom, "chr" + chrom.id)
        # Mt
        if chrom.id == "MT":
            add_if_unique(chrom, "chr" + "M")


def _check_intervals_array_shape(intervals):
    if (len(intervals.shape) != 2) or (intervals.shape[1] < 2):
        raise ValueError(
            "Intervals must be 2D numpy arrays with 2 columns " "[left, right)."
        )


def _check_intervals_validity(intervals, start=0, end=np.inf):
    # is numpy?
    if not isinstance(intervals, np.ndarray):
        raise ValueError("Intervals must be a numpy array, did you try to pass a list?")
    # are the intervals ints?
    if (intervals.dtype != "int32") and (intervals.dtype != "int64"):
        raise ValueError(
            "Intervals must be arrays of type integer" "(i.e., with dtype='int')."
        )
    _check_intervals_array_shape(intervals)
    if np.any(intervals[:, 0] >= intervals[:, 1]):
        raise ValueError(
            "Left positions cannot be greater than or equal to right positions."
        )
    if np.any(intervals[:, :2] > end) or np.any(intervals[:, :2] < start):
        raise ValueError(f"Intervals must be within [{start}, {end}]")
    if not np.all(intervals[1:, 0] >= intervals[:-1, 1]):
        raise ValueError("Intervals must be non-overlapping.")


def mask_intervals(intervals, mask):
    """
    Removes the intervals of ``mask`` from those in ``intervals``.
    Both ``intervals`` and ``mask`` should be sets of intervals
    (i.e., sorted numpy arrays with two columns); then the result will
    be the the intervals in ``intervals`` except with those in
    ``mask`` removed.
    """
    _check_intervals_validity(intervals)
    _check_intervals_validity(mask)
    out = []
    last_mask_right = -1 * np.inf
    j = 0
    if j < len(mask):
        next_mask_left, next_mask_right = mask[j]
    else:
        # the mask is of zero length, return for efficiency
        return intervals
    for inter in intervals:
        left, right = inter
        while left < right:
            while left >= next_mask_left:
                last_mask_right = next_mask_right
                j += 1
                if j < len(mask):
                    next_mask_left, next_mask_right = mask[j]
                else:
                    next_mask_left, next_mask_right = np.inf, np.inf
            next_left = max(left, last_mask_right)
            next_right = min(right, next_mask_left)
            if next_left < next_right:
                out.append([next_left, next_right])
            left = max(next_left, next_right)
    out = np.array(out, dtype=intervals.dtype).reshape((len(out), 2))
    return out
