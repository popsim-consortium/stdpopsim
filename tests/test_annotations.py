"""
Tests for the annotations management.
"""
import unittest
from unittest import mock
import gzip
import tempfile
import os.path
import shutil
import urllib.request
import pathlib

import pandas
import stdpopsim
from stdpopsim import annotations
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
        key = an.file_name
        local_file = destination / key
        if not local_file.exists():
            cache_dir = local_file.parent
            cache_dir.mkdir(exist_ok=True, parents=True)
            print("Downloading", an.zarr_url)
            urllib.request.urlretrieve(an.zarr_url, local_file)
        saved_urls[key] = an.zarr_url
        an.zarr_url = local_file.resolve().as_uri()


def tearDownModule():
    for an in stdpopsim.all_annotations():
        an.zarr_url = saved_urls[an.file_name]


class AnnotationTestClass(annotations.Annotation):
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
            file_name="annotation.gff.gz")


# number of chromosomes etc.
def get_annotation_file(custom_file_f=None, filter=None):
    """
    Returns a GFF file as a bytes object.

    :param func custom_file_f: A function that accepts a single parameter
        (a folder name), which may be used to create additional files under
        the given folder. All files in the folder will be included in the
        returned tarball.
    :param func filter: A function which is passed as the ``filter`` argument
        to ``TarFile.add()``. This function can be used to change the info
        field for each file in the returned tarball. See tarfile documentation
        for more details.
    """
    b = "21\t.\tbiological_region\t6111184\t6111187\t0.999\t+\t.\tlogic_name=eponine\n"
    b += b
    b += b

    with tempfile.TemporaryDirectory() as map_dir:
        with open(os.path.join(map_dir, "test.gff"), "w") as f:
            print("##gff-version 3", file=f)
            print("##sequence-region   21 1 46709983", file=f)
            print(b, file=f)

        if custom_file_f is not None:
            # Do arbitrary things under map_dir.
            custom_file_f(map_dir)

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
        gm = AnnotationTestClass()
        cache_dir = stdpopsim.get_cache_dir() / "annotations"
        self.assertEqual(gm.annot_cache_dir, cache_dir)
        self.assertEqual(gm.species_cache_dir, gm.annot_cache_dir / gm.species.id)

    def test_str(self):
        gm = AnnotationTestClass()
        self.assertGreater(len(str(gm)), 0)

    def test_is_cached(self):
        gm = AnnotationTestClass()
        os.makedirs(gm.species_cache_dir, exist_ok=True)
        self.assertTrue(gm.is_cached())
        shutil.rmtree(gm.annot_cache_dir)
        self.assertFalse(gm.is_cached())


class TestAnnotationDownload(tests.CacheWritingTest):
    """
    Tests downloading code for the annotations.
    """

    def test_correct_url(self):
        gm = AnnotationTestClass()
        with mock.patch("urllib.request.urlretrieve", autospec=True) as mocked_get:
            gm.download()
        mocked_get.assert_called_once_with(gm.zarr_url, filename=unittest.mock.ANY)

    def test_incorrect_url(self):
        gm = AnnotationTestClass()
        gm.zarr_url = 'http://asdfwersdf.com/foozip'
        with self.assertRaises(urllib.error.URLError):
            gm.download()

    def test_download_over_cache(self):
        for gm in stdpopsim.all_annotations():
            gm.download()
            self.assertTrue(gm.is_cached())
            gm.download()
            self.assertTrue(gm.is_cached())

    def test_multiple_threads_downloading(self):
        gm = next(stdpopsim.all_annotations())
        gm.download()
        saved = gm.is_cached
        try:
            # Trick the download code into thinking there's several happening
            # concurrently
            gm.is_cached = lambda: False
            with self.assertWarns(Warning):
                gm.download()
        finally:
            gm.is_cached = saved


class TestGetChromosomeAnnotations(tests.CacheReadingTest):
    """
    Tests if we get chromosome level annotations
    using the Ensembl_GRCh38 human GFF.
    """
    species = stdpopsim.get_species("HomSap")
    an = species.get_annotations("Ensembl_GRCh38_gff3")

    def test_known_chromosome(self):
        print(self.cache_dir)
        cm = self.an.get_chromosome_annotations("21")
        self.assertIsInstance(cm, pandas.DataFrame)

    def test_known_chromosome_prefix(self):
        print(self.cache_dir)
        cm = self.an.get_chromosome_annotations("chr21")
        self.assertIsInstance(cm, pandas.DataFrame)

    def test_unknown_chromosome(self):
        for bad_chrom in ["", "ABD", None]:
            with self.assertRaises(ValueError):
                self.an.get_chromosome_annotations(bad_chrom)

    def test_unknown_annot_file(self):
        real_fn = self.an.file_name
        self.an.file_name = "foo"
        bad_chrom = "chrF"
        with self.assertRaises(ValueError):
            self.an.get_chromosome_annotations(bad_chrom)
        self.an.file_name = real_fn

    def test_get_genes(self):
        g = self.an.get_genes_from_chromosome("chr21")
        self.assertIsInstance(g, pandas.DataFrame)

    def test_get_genes_full(self):
        g = self.an.get_genes_from_chromosome("chr21", full_table=True)
        self.assertIsInstance(g, pandas.DataFrame)

    def test_bad_annot_type(self):
        bad_annot = "foo"
        with self.assertRaises(ValueError):
            self.an.get_annotation_type_from_chromomosome(bad_annot, '21')
