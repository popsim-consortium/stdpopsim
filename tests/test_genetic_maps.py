"""
Tests for the genetic maps management.
"""
import unittest
from unittest import mock
import tarfile
import tempfile
import os.path
import shutil
import urllib.request
import pathlib

import msprime

import stdpopsim
from stdpopsim import genetic_maps
import tests


# Here we download all the genetic map tarballs in one go and store
# then in the local cache directory, _test_cache. Tests are then run
# with the download URLs redirected to local files, which makes them
# much faster and takes network errors out of the equation. The
# tarball cache is also useful for developers as it means that these
# files are only downloaded once.

saved_urls = {}


def setUpModule():
    destination = pathlib.Path("_test_cache/tarballs")
    for genetic_map in stdpopsim.all_genetic_maps():
        key = genetic_map.name
        local_file = destination / (key + ".tar.gz")
        if not local_file.exists():
            cache_dir = local_file.parent
            cache_dir.mkdir(exist_ok=True, parents=True)
            # print("Downloading", genetic_map.url)
            urllib.request.urlretrieve(genetic_map.url, local_file)
        saved_urls[key] = genetic_map.url
        genetic_map.url = local_file.resolve().as_uri()


def tearDownModule():
    for genetic_map in stdpopsim.all_genetic_maps():
        genetic_map.url = saved_urls[genetic_map.name]


class GeneticMapTestClass(genetic_maps.GeneticMap):
    """
    A genetic map that we can instantiate to get a genetic map for testing.
    """

    def __init__(self):
        genome = stdpopsim.Genome(chromosomes=[])
        _species = stdpopsim.Species(
            id="TesSpe", name="Test species", common_name="Testy McTestface",
            genome=genome)
        super().__init__(
            species=_species,
            name="test_map",
            url="http://example.com/genetic_map.tar.gz",
            file_pattern="prefix_{name}.txt")


# TODO add some parameters here to check different compression options,
# number of chromosomes etc.
def get_genetic_map_tarball(custom_file_f=None, filter=None):
    """
    Returns a genetic map in hapmap format in a tarball as a bytes object.

    :param func custom_file_f: A function that accepts a single parameter
        (a folder name), which may be used to create additional files under
        the given folder. All files in the folder will be included in the
        returned tarball.
    :param func filter: A function which is passed as the ``filter`` argument
        to ``TarFile.add()``. This function can be used to change the info
        field for each file in the returned tarball. See tarfile documentation
        for more details.
    """
    with tempfile.TemporaryDirectory() as map_dir:
        for j in range(1, 10):
            # TODO Have a way to put in different maps??
            with open(os.path.join(map_dir, "prefix_chr{}.txt".format(j)), "w") as f:
                print("Chromosome  Position(bp)    Rate(cM/Mb)     Map(cM)", file=f)
                print("chr1        55550   2.981822        0.000000", file=f)
                print("chr1        82571   2.082414        0.080572", file=f)
                print("chr1        88169   0               0.092229", file=f)

        if custom_file_f is not None:
            # Do arbitrary things under map_dir.
            custom_file_f(map_dir)

        # For the tarfile to be in the right format, we must be in the right directory.
        with genetic_maps.cd(map_dir):
            # Now tar up this map_directory
            with tempfile.TemporaryFile('wb+') as tmp_file:
                with tarfile.open(fileobj=tmp_file, mode="w:gz") as tar_file:
                    for filename in os.listdir("."):
                        tar_file.add(filename, filter=filter)
                # Read back the tarball
                tmp_file.seek(0)
                tarball = tmp_file.read()
    return tarball


class TestGeneticMapTarball(unittest.TestCase):
    """
    Tests that we correctly encode a genetic map in the tarball test function.
    """
    def get_maps(self, tarball):
        maps = {}
        with tempfile.TemporaryFile('wb+') as f:
            f.write(tarball)
            f.seek(0)
            with tarfile.open(fileobj=f, mode='r') as tar_file:
                with tempfile.TemporaryDirectory() as extract_dir:
                    with genetic_maps.cd(extract_dir):
                        tar_file.extractall()
                        for fn in os.listdir(extract_dir):
                            maps[fn] = msprime.RecombinationMap.read_hapmap(fn)
        return maps

    def test_no_args(self):
        tarball = get_genetic_map_tarball()
        maps = self.get_maps(tarball)
        self.assertGreater(len(maps), 0)


class TestGeneticMap(tests.CacheWritingTest):
    """
    Tests for the basic functionality of the genetic map class.
    """

    def test_cache_dirs(self):
        gm = GeneticMapTestClass()
        cache_dir = stdpopsim.get_cache_dir() / "genetic_maps"
        self.assertEqual(gm.cache_dir, cache_dir)
        self.assertEqual(gm.species_cache_dir, gm.cache_dir / gm.species.id)
        self.assertEqual(gm.map_cache_dir, gm.species_cache_dir / gm.name)

    def test_str(self):
        gm = GeneticMapTestClass()
        self.assertGreater(len(str(gm)), 0)

    def test_is_cached(self):
        gm = GeneticMapTestClass()
        os.makedirs(gm.map_cache_dir, exist_ok=True)
        self.assertTrue(gm.is_cached())
        shutil.rmtree(gm.map_cache_dir)
        self.assertFalse(gm.is_cached())


class TestGeneticMapDownload(tests.CacheWritingTest):
    """
    Tests downloading code for the genetic maps.
    """

    def test_correct_url(self):
        gm = GeneticMapTestClass()
        with mock.patch("urllib.request.urlretrieve") as mocked_get:
            # The destination file will be missing.
            with self.assertRaises(FileNotFoundError):
                gm.download()
        mocked_get.assert_called_once_with(gm.url, filename=unittest.mock.ANY)

    def test_download_over_cache(self):
        for gm in stdpopsim.all_genetic_maps():
            gm.download()
            self.assertTrue(gm.is_cached())
            gm.download()
            self.assertTrue(gm.is_cached())

    def test_multiple_threads_downloading(self):
        gm = next(stdpopsim.all_genetic_maps())
        gm.download()
        saved = gm.is_cached
        try:
            # Trick the download code into thinking there's several happening
            # concurrently
            gm.is_cached = lambda: False
            with self.assertWarns(UserWarning):
                gm.download()
        finally:
            gm.is_cached = saved


class TestGeneticMapDownloadSecurity(tests.CacheWritingTest):
    """
    Security related tests for the genetic map downloading code.

    For an extended discussion of tar-related security issues see:
    https://bugs.python.org/issue21109
    """

    def test_tar_path_traversal_attack(self):
        # Test for vulnerability to path-traversal attacks.
        def mock_retrieve_factory():
            for dest in ("../nonexistant", "/nonexistant"):
                def retrieve(url, filename):
                    tarball = get_genetic_map_tarball(
                            custom_file_f=lambda map_dir: os.symlink(
                                dest, os.path.join(map_dir, "my-link")))
                    with open(filename, "wb") as f:
                        f.write(tarball)
                yield retrieve, f"path-traversal: {dest}"
        self.assert_bad_tar(mock_retrieve_factory)

    def test_bad_tar_members(self):
        # Pretend we downloaded a tarball containing a FIFO or device file.
        # There is no reasonable use for these types of files in a genetic
        # map tarball, and their presence likely indicates a maliciously
        # crafted tarball.
        # Creating a character or block special device file requires root
        # privileges, so we instead modify the ``type`` field of each file
        # in the tar.
        def mock_retrieve_factory():
            for file_type, assert_msg in zip(
                    (tarfile.FIFOTYPE, tarfile.CHRTYPE, tarfile.BLKTYPE),
                    ("FIFO", "char device", "block device")):
                def retrieve(url, filename):
                    def filt(info):
                        info.type = file_type
                        return info
                    tarball = get_genetic_map_tarball(filter=filt)
                    with open(filename, "wb") as f:
                        f.write(tarball)
                yield retrieve, assert_msg

        self.assert_bad_tar(mock_retrieve_factory)

    def assert_bad_tar(self, mock_retrieve_factory):
        gm = GeneticMapTestClass()
        for retrieve, assert_msg in mock_retrieve_factory():
            with mock.patch("urllib.request.urlretrieve", new=retrieve):
                with self.assertRaises(ValueError, msg=assert_msg):
                    gm.download()
            self.assertFalse(gm.is_cached())


class TestAllGeneticMaps(tests.CacheReadingTest):
    """
    Tests if the all_genetic_maps() function works correctly.
    """
    def test_non_empty(self):
        self.assertGreater(len(list(stdpopsim.all_genetic_maps())), 0)

    def test_types(self):
        for gm in stdpopsim.all_genetic_maps():
            self.assertIsInstance(gm, genetic_maps.GeneticMap)


class TestGetChromosomeMap(tests.CacheReadingTest):
    """
    Tests if we get chromosome maps using the HapmapII_GRCh37 human map.
    """
    species = stdpopsim.get_species("HomSap")
    genetic_map = species.get_genetic_map("HapmapII_GRCh37")

    def test_warning_from_no_mapped_chromosome(self):
        chrom = self.species.genome.get_chromosome("chrY")
        with self.assertWarns(Warning):
            cm = self.genetic_map.get_chromosome_map(chrom.id)
            self.assertIsInstance(cm, msprime.RecombinationMap)
            self.assertEqual(chrom.length, cm.get_sequence_length())

    def test_known_chromosome(self):
        chrom = self.species.genome.get_chromosome("chr22")
        cm = self.genetic_map.get_chromosome_map(chrom.id)
        self.assertIsInstance(cm, msprime.RecombinationMap)

    def test_unknown_chromosome(self):
        for bad_chrom in ["", "ABD", None]:
            with self.assertRaises(ValueError):
                self.genetic_map.get_chromosome_map(bad_chrom)
