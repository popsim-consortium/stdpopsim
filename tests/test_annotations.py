"""
Tests for the annotations management.
"""
import unittest
from unittest import mock
import gzip
import tempfile
import os.path
import shutil
import pathlib

import pandas

import stdpopsim
from stdpopsim import utils
import tests


# Here we download all the annotations in one go and store
# then in the local cache directory, _test_cache. Tests are then run
# with the download URLs redirected to local files, which makes them
# much faster and takes network errors out of the equation. The
# tarball cache is also useful for developers as it means that these
# files are only downloaded once.

saved_urls = {}


def setUpModule():
    destination = pathlib.Path("_test_cache/zipfiles/")
    for an in stdpopsim.all_annotations():
        key = an._cache.cache_path
        local_file = destination / key
        if not local_file.exists():
            cache_dir = local_file.parent
            cache_dir.mkdir(exist_ok=True, parents=True)
            print("Downloading", an.zarr_url)
            utils.download(an.zarr_url, local_file)
        # This assertion could fail if we update a file on AWS,
        # or a developer creates a new annotation with the wrong checksum
        # (in the latter case, this should at least be caught during CI tests).
        assert utils.sha256(local_file) == an.zarr_sha256, (
                f"SHA256 for {local_file} doesn't match the SHA256 for "
                f"{an.id}. If you didn't add this SHA256 yourself, "
                f"try deleting {local_file} and restarting the tests."
            )
        saved_urls[key] = an.zarr_url
        an.zarr_url = local_file.resolve().as_uri()
        an._cache.url = an.zarr_url


def tearDownModule():
    for an in stdpopsim.all_annotations():
        an.zarr_url = saved_urls[an._cache.cache_path]
        an._cache.url = an.zarr_url


class AnnotationTestClass(stdpopsim.Annotation):
    """
    A Annotation that we can instantiate to get an annotation for testing.
    """

    def __init__(self):
        genome = stdpopsim.Genome(chromosomes=[])
        _species = stdpopsim.Species(
            id="TesSpe", name="Test species", common_name="Testy McTestface",
            genome=genome)
        super().__init__(
            species=_species,
            id="test_annotation",
            url="http://example.com/annotation.gff.gz",
            zarr_url="http://example.com/annotation.zip",
            zarr_sha256="1234",  # this shouldn't be checked anywhere
            description="test annotation",
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
        with open(os.path.join(map_dir, "test.gff"), 'rb') as f_in:
            # Now zip this
            with gzip.open("test.gff.gz", mode="wb") as gz_file:
                shutil.copyfileobj(f_in, gz_file)
        # Read back the tarball
        an = gzip.open("test.gff.gz", mode='rb')
    return an


class TestAnnotation(tests.CacheWritingTest):
    """
    Tests for the basic functionality of the Annotation class.
    """

    def test_cache_dirs(self):
        an = AnnotationTestClass()
        cache_dir = stdpopsim.get_cache_dir() / "annotations" / an.species.id
        self.assertEqual(an.cache_path.parent, cache_dir)

    def test_str(self):
        an = AnnotationTestClass()
        self.assertGreater(len(str(an)), 0)


class TestAnnotationDownload(tests.CacheWritingTest):
    """
    Tests downloading code for the annotations.
    """

    def test_correct_url(self):
        an = AnnotationTestClass()
        with mock.patch("stdpopsim.utils.download", autospec=True) as mocked_get:
            # The destination file will be missing.
            with self.assertRaises(FileNotFoundError):
                an.download()
        mocked_get.assert_called_once_with(an.zarr_url, filename=unittest.mock.ANY)

    def test_incorrect_url(self):
        an = AnnotationTestClass()
        an.zarr_url = 'http://asdfwersdf.com/foozip'
        with self.assertRaises(OSError):
            an.download()

    def test_download_over_cache(self):
        # TODO: The HomSap annotations are huge. Once we include a smaller
        # annotation set, we should instead use that, so tests are faster.
        species = stdpopsim.get_species("HomSap")
        an = species.get_annotations("Ensembl_GRCh38_gff3")
        an.download()
        self.assertTrue(an.is_cached())
        an.download()
        self.assertTrue(an.is_cached())


class TestGetChromosomeAnnotations(tests.CacheReadingTest):
    """
    Tests if we get chromosome level annotations
    using the Ensembl_GRCh38 human GFF.
    """
    # TODO: The HomSap annotations are huge. Once we include a smaller
    # annotation set, we should instead use that, so tests are faster.
    species = stdpopsim.get_species("HomSap")
    an = species.get_annotations("Ensembl_GRCh38_gff3")

    def test_known_chromosome(self):
        cm = self.an.get_chromosome_annotations("21")
        self.assertIsInstance(cm, pandas.DataFrame)

    def test_known_chromosome_prefix(self):
        cm = self.an.get_chromosome_annotations("chr21")
        self.assertIsInstance(cm, pandas.DataFrame)

    def test_unknown_chromosome(self):
        for bad_chrom in ["", "ABD", None]:
            with self.assertRaises(ValueError):
                self.an.get_chromosome_annotations(bad_chrom)

    def test_get_genes(self):
        g = self.an.get_genes_from_chromosome("21")
        self.assertIsInstance(g, pandas.DataFrame)

    def test_get_genes_full(self):
        g = self.an.get_genes_from_chromosome("21", full_table=True)
        self.assertIsInstance(g, pandas.DataFrame)

    def test_bad_annot_type(self):
        bad_annot = "foo"
        with self.assertRaises(ValueError):
            self.an.get_annotation_type_from_chromomosome(bad_annot, "21")
