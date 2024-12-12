"""
Tests for the annotations management.
"""

from unittest import mock
import gzip
import tempfile
import os.path
import shutil
import pathlib
import pytest
import numpy as np
import logging
import stdpopsim
from stdpopsim import utils
import tests
import maintenance as maint

# Here we download all the annotations in one go and store
# then in the local cache directory, _test_cache. Tests are then run
# with the download URLs redirected to local files, which makes them
# much faster and takes network errors out of the equation. The
# tarball cache is also useful for developers as it means that these
# files are only downloaded once.

saved_urls = {}


def setup_module():
    destination = pathlib.Path("_test_cache/tarballs")
    for an in stdpopsim.all_annotations():
        key = an.id
        local_file = destination / (key + ".tar.gz")
        logging.info(f"key {key} local_file {local_file}")
        if not local_file.exists():
            cache_dir = local_file.parent
            cache_dir.mkdir(exist_ok=True, parents=True)
            print("Downloading", an.intervals_url)
            utils.download(an.intervals_url, local_file)
        # This assertion could fail if we update a file on AWS,
        # or a developer creates a new annotation with the wrong checksum
        # (in the latter case, this should at least be caught during CI tests).
        assert utils.sha256(local_file) == an.intervals_sha256, (
            f"SHA256 for {local_file} doesn't match the SHA256 for "
            f"{an.id}. If you didn't add this SHA256 yourself, "
            f"try deleting {local_file} and restarting the tests."
        )
        saved_urls[key] = an.intervals_url
        an.intervals_url = local_file.resolve().as_uri()
        an._cache.url = an.intervals_url


def teardown_module():
    for an in stdpopsim.all_annotations():
        an.intervals_url = saved_urls[an.id]
        an._cache.url = an.intervals_url


class AnnotationTestClass(stdpopsim.Annotation):
    """
    A Annotation that we can instantiate to get an annotation for testing.
    """

    def __init__(self):
        genome = stdpopsim.Genome(chromosomes=[])
        _species = stdpopsim.Species(
            id="TesSpe",
            ensembl_id="test_species",
            name="Test species",
            common_name="Testy McTestface",
            genome=genome,
        )
        super().__init__(
            species=_species,
            id="test_annotation",
            url="http://example.com/annotation.gff.gz",
            intervals_url="http://example.com/annotation.zip",
            intervals_sha256="1234",  # this shouldn't be checked anywhere
            gff_sha256="6789",
            description="test annotation",
            file_pattern="yolo_{id}.txt",
            annotation_source="your mom",
            annotation_type="test",
        )


# TODO: use this function for testing, or remove it.
def get_annotation_file():
    """
    Returns a GFF file as a bytes object.
    """
    b = "21\t.\tbiological_region\t6111184\t6111187\t0.999\t+\t.\tlogic_name=eponine\n"
    b += b
    b += b

    with tempfile.TemporaryDirectory() as map_dir:
        with open(os.path.join(map_dir, "test.gff"), "w") as f:
            print("##gff-version 3", file=f)
            print("##sequence-region   21 1 46709983", file=f)
            print(b, file=f)

        # For the zipfile to be in the right format, we must be in the right directory.
        with open(os.path.join(map_dir, "test.gff"), "rb") as f_in:
            # Now zip this
            with gzip.open("test.gff.gz", mode="wb") as gz_file:
                shutil.copyfileobj(f_in, gz_file)
        # Read back the tarball
        an = gzip.open("test.gff.gz", mode="rb")
    return an


class TestAnnotation(tests.CacheWritingTest):
    """
    Tests for the basic functionality of the Annotation class.
    """

    def test_cache_dirs(self):
        an = AnnotationTestClass()
        cache_dir = stdpopsim.get_cache_dir() / "annotations" / an.species.id
        assert an.cache_path.parent == cache_dir

    def test_str(self):
        an = AnnotationTestClass()
        assert len(str(an)) > 0


class TestAnnotationDownload(tests.CacheWritingTest):
    """
    Tests downloading code for the annotations.
    """

    def test_correct_url(self):
        an = AnnotationTestClass()
        with mock.patch("stdpopsim.utils.download", autospec=True) as mocked_get:
            # The destination file will be missing.
            with pytest.raises(FileNotFoundError):
                an.download()
        mocked_get.assert_called_once_with(an.intervals_url, filename=mock.ANY)

    def test_incorrect_url(self):
        an = AnnotationTestClass()
        an.intervals_url = "http://asdfwersdf.com/foozip"
        with pytest.raises(OSError):
            an.download()

    def test_download_over_cache(self):
        # TODO: The HomSap annotations are huge. Once we include a smaller
        # annotation set, we should instead use that, so tests are faster.
        species = stdpopsim.get_species("HomSap")
        an = species.get_annotations("ensembl_havana_104_CDS")
        an.download()
        assert an.is_cached()
        an.download()
        assert an.is_cached()


class TestGetChromosomeAnnotations(tests.CacheReadingTest):
    """
    Tests if we get chromosome level annotations
    using the Ensembl_GRCh38 human GFF.
    """

    @classmethod
    def setup_class(cls):
        species = stdpopsim.get_species("HomSap")
        cls.an = species.get_annotations("ensembl_havana_104_exons")

    def test_known_chromosome(self):
        cm = self.an.get_chromosome_annotations("21")
        assert isinstance(cm, np.ndarray)

    def test_known_chromosome_prefix(self):
        cm = self.an.get_chromosome_annotations("chr21")
        assert isinstance(cm, np.ndarray)

    def test_unknown_chromosome(self):
        for bad_chrom in ["", "ABD", None]:
            with pytest.raises(ValueError):
                self.an.get_chromosome_annotations(bad_chrom)


class TestGetChromosomeAnnotationsDroMel(tests.CacheReadingTest):
    """
    Tests if we get chromosome level annotations
    using the FlyBase_BDGP6.32.51_CDS GFF.
    """

    @classmethod
    def setup_class(cls):
        species = stdpopsim.get_species("DroMel")
        cls.an = species.get_annotations("FlyBase_BDGP6.32.51_CDS")

    def test_known_chromosome(self):
        cm = self.an.get_chromosome_annotations("2L")
        assert isinstance(cm, np.ndarray)

    def test_known_chromosome_prefix(self):
        cm = self.an.get_chromosome_annotations("chr2L")
        assert isinstance(cm, np.ndarray)

    def test_unknown_chromosome(self):
        for bad_chrom in ["", "ABD", None]:
            with pytest.raises(ValueError):
                self.an.get_chromosome_annotations(bad_chrom)


class TestGetChromosomeAnnotationsAraTha(tests.CacheReadingTest):
    """
    Tests if we get chromosome level annotations
    using the ara GFF.
    """

    @classmethod
    def setup_class(cls):
        species = stdpopsim.get_species("AraTha")
        cls.an = species.get_annotations("araport_11_exons")

    def test_known_chromosome(self):
        cm = self.an.get_chromosome_annotations("2")
        assert isinstance(cm, np.ndarray)

    def test_known_chromosome_prefix(self):
        cm = self.an.get_chromosome_annotations("chr2")
        assert isinstance(cm, np.ndarray)

    def test_unknown_chromosome(self):
        for bad_chrom in ["", "ABDEFG", None]:
            with pytest.raises(ValueError):
                self.an.get_chromosome_annotations(bad_chrom)


class TestIntervalMerge:
    """
    Tests for the Annotation stuff.
    """

    def test_merged_static(self):
        for closed in (True, False):
            assert maint.merged([], closed=closed) == []
            assert maint.merged([(10, 20), (15, 30)], closed=closed) == [(10, 30)]
            assert maint.merged([(10, 20), (20, 30)], closed=closed) == [(10, 30)]
            assert maint.merged([(10, 20), (22, 30)], closed=closed) == [
                (10, 20),
                (22, 30),
            ]
            assert maint.merged([(10, 20), (12, 16)], closed=closed) == [(10, 20)]
            assert maint.merged([(12, 16), (10, 20), (13, 15)], closed=closed) == [
                (10, 20)
            ]

        assert maint.merged([(10, 20), (21, 30)], closed=True) == [(10, 30)]
        assert maint.merged([(10, 20), (21, 30)], closed=False) == [(10, 20), (21, 30)]

    def test_merged_random(self):
        rng = np.random.default_rng(1234)
        for closed in (True, False):
            # Check merging is idempotent.
            for _ in range(100):
                starts = rng.integers(1, 1000, size=100)
                ends = starts + rng.integers(1, 100, size=len(starts))
                merged_intervals = maint.merged(zip(starts, ends), closed=closed)
                assert merged_intervals == maint.merged(merged_intervals, closed=closed)
