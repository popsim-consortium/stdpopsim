"""
Tests for the genetic maps management.
"""
import unittest
from unittest import mock
import tarfile
import tempfile
import os.path
import shutil

import appdirs
import msprime

from stdpopsim import genetic_maps

_cachedir = None


def setUpModule():
    global _cachedir
    _cachedir = tempfile.TemporaryDirectory()


def tearDownModule():
    global _cachedir
    del _cachedir


class GeneticMapTestClass(genetic_maps.GeneticMap):
    species = "test_species"
    name = "test_map"
    url = "http://example.com/genetic_map.tar.gz"
    file_pattern = "prefix_{name}.txt"


class GeneticMapCacheTestClass(GeneticMapTestClass):
    """
    Ensure we're not touching the global cache directory.
    """
    # def __init__(self):
    #     self.cache_dir = _cachedir.name
    #     self.species_cache_dir = os.path.join(self.cache_dir, self.species)
    #     self.map_cache_dir = os.path.join(self.species_cache_dir, self.name)


# TODO add some parameters here to check different compression options,
# number of chromosomes etc.
def get_genetic_map_tarball():
    """
    Returns a genetic map in hapmap format in a tarball as a bytes object.
    """
    with tempfile.TemporaryDirectory() as map_dir:
        for j in range(1, 10):
            # TODO Have a way to put in different maps??
            with open(os.path.join(map_dir, "prefix_chr{}.txt".format(j)), "w") as f:
                print("Chromosome  Position(bp)    Rate(cM/Mb)     Map(cM)", file=f)
                print("chr1        55550   2.981822        0.000000", file=f)
                print("chr1        82571   2.082414        0.080572", file=f)
                print("chr1        88169   0               0.092229", file=f)
        # For the tarfile to be in the right format, we must be in the right directory.
        with genetic_maps.cd(map_dir):
            # Now tar up this map_directory
            with tempfile.TemporaryFile('wb+') as tmp_file:
                with tarfile.open(fileobj=tmp_file, mode="w:gz") as tar_file:
                    for filename in os.listdir("."):
                        tar_file.add(filename)
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


class TestGeneticMap(unittest.TestCase):
    """
    Tests for the basic functionality of the genetic map class.
    """

    def test_cache_dirs(self):
        gm = GeneticMapTestClass()
        cache_dir = appdirs.user_cache_dir("stdpopsim", "popgensims")
        self.assertEqual(gm.cache_dir, cache_dir)
        self.assertEqual(gm.species_cache_dir, os.path.join(gm.cache_dir, gm.species))
        self.assertEqual(gm.map_cache_dir, os.path.join(gm.species_cache_dir, gm.name))

    def test_str(self):
        gm = GeneticMapTestClass()
        self.assertGreater(len(str(gm)), 0)

    def test_is_cached(self):
        gm = GeneticMapTestClass()
        os.makedirs(gm.map_cache_dir, exist_ok=True)
        self.assertTrue(gm.is_cached())
        shutil.rmtree(gm.map_cache_dir)
        self.assertFalse(gm.is_cached())


class TestGeneticMapDownload(unittest.TestCase):

    def test_correct_url(self):
        gm = GeneticMapCacheTestClass()
        with mock.patch("requests.get") as mocked_get:
            # We're trying to untar an empty file, which fails.
            with self.assertRaises(tarfile.ReadError):
                gm.download()
        mocked_get.assert_called_once_with(gm.url, stream=True)

    def test_file_download(self):
        gm = GeneticMapCacheTestClass()
        tarball = get_genetic_map_tarball()

        def iter_content(chunk_size):
            yield tarball
        with mock.patch("requests.get") as mocked_get:
            mocked_get.return_value.iter_content = iter_content
            gm.download()
        self.assertTrue(gm.is_cached())
        recomb_map = gm.get_chromosome_map("chr1")
        # TODO add some real tests for a recomb_map.
        self.assertIsNotNone(recomb_map)
