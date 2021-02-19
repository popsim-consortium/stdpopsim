"""
Tests for the genetic maps management.
"""
import unittest
from unittest import mock
import tarfile
import tempfile
import os.path
import pathlib

import msprime

import stdpopsim
from stdpopsim import utils
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
        key = genetic_map.id
        local_file = destination / (key + ".tar.gz")
        if not local_file.exists():
            cache_dir = local_file.parent
            cache_dir.mkdir(exist_ok=True, parents=True)
            print("Downloading", genetic_map.url)
            utils.download(genetic_map.url, local_file)
        # This assertion could fail if we update a file on AWS,
        # or a developer creates a new genetic map with the wrong checksum
        # (in the latter case, this should at least be caught during CI tests).
        assert utils.sha256(local_file) == genetic_map.sha256, (
            f"SHA256 for {local_file} doesn't match the SHA256 for "
            f"{genetic_map.id}. If you didn't add this SHA256 yourself, "
            f"try deleting {local_file} and restarting the tests."
        )
        saved_urls[key] = genetic_map.url
        genetic_map.url = local_file.resolve().as_uri()
        genetic_map._cache.url = genetic_map.url


def tearDownModule():
    for genetic_map in stdpopsim.all_genetic_maps():
        genetic_map.url = saved_urls[genetic_map.id]
        genetic_map._cache.url = genetic_map.url


class GeneticMapTestClass(stdpopsim.GeneticMap):
    """
    A genetic map that we can instantiate to get a genetic map for testing.
    """

    def __init__(self):
        genome = stdpopsim.Genome(chromosomes=[])
        _species = stdpopsim.Species(
            id="TesSpe",
            name="Test species",
            common_name="Testy McTestface",
            genome=genome,
        )
        super().__init__(
            species=_species,
            id="test_map",
            url="http://example.com/genetic_map.tar.gz",
            sha256="1234",  # url doesn't exist, so this will never be checked
            file_pattern="prefix_{name}.txt",
        )


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
        with utils.cd(map_dir):
            # Now tar up this map_directory
            with tempfile.TemporaryFile("wb+") as tmp_file:
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
        with tempfile.TemporaryFile("wb+") as f:
            f.write(tarball)
            f.seek(0)
            with tarfile.open(fileobj=f, mode="r") as tar_file:
                with tempfile.TemporaryDirectory() as extract_dir:
                    with utils.cd(extract_dir):
                        tar_file.extractall()
                        for fn in os.listdir(extract_dir):
                            maps[fn] = msprime.RateMap.read_hapmap(fn)
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
        cache_dir = stdpopsim.get_cache_dir() / "genetic_maps" / gm.species.id / gm.id
        self.assertEqual(gm.map_cache_dir, cache_dir)

    def test_str(self):
        gm = GeneticMapTestClass()
        self.assertGreater(len(str(gm)), 0)


class TestGeneticMapDownload(tests.CacheWritingTest):
    """
    Tests downloading code for the genetic maps.
    """

    def test_correct_url(self):
        gm = GeneticMapTestClass()
        with mock.patch("stdpopsim.utils.download", autospec=True) as mocked_get:
            # The destination file will be missing.
            with self.assertRaises(FileNotFoundError):
                gm.download()
        mocked_get.assert_called_once_with(gm.url, unittest.mock.ANY)

    def test_download_over_cache(self):
        species = stdpopsim.get_species("DroMel")
        gm = species.get_genetic_map("ComeronCrossover_dm6")
        gm.download()
        self.assertTrue(gm.is_cached())
        gm.download()
        self.assertTrue(gm.is_cached())

    def test_multiple_threads_downloading(self):
        species = stdpopsim.get_species("DroMel")
        gm = species.get_genetic_map("ComeronCrossover_dm6")
        gm.download()
        saved = gm._cache.is_cached
        try:
            # Trick the download code into thinking there's several happening
            # concurrently
            gm._cache.is_cached = lambda: False
            with self.assertWarns(UserWarning):
                gm.download()
        finally:
            gm._cache.is_cached = saved


class TestAllGeneticMaps(tests.CacheReadingTest):
    """
    Tests if the all_genetic_maps() function works correctly.
    """

    def test_non_empty(self):
        self.assertGreater(len(list(stdpopsim.all_genetic_maps())), 0)

    def test_types(self):
        for gm in stdpopsim.all_genetic_maps():
            self.assertIsInstance(gm, stdpopsim.GeneticMap)

    def test_ids(self):
        for gm in stdpopsim.all_genetic_maps():
            self.assertIsInstance(gm.id, str)
            self.assertTrue(utils.is_valid_genetic_map_id(gm.id))


class TestGetChromosomeMap(tests.CacheReadingTest):
    """
    Tests if we get chromosome maps using the HapMapII_GRCh37 human map.
    """

    def test_warning_from_no_mapped_chromosome(self):
        species = stdpopsim.get_species("HomSap")
        genetic_map = species.get_genetic_map("HapMapII_GRCh37")
        chrom = species.genome.get_chromosome("chrY")
        with self.assertWarns(Warning):
            cm = genetic_map.get_chromosome_map(chrom.id)
        self.assertIsInstance(cm, msprime.RateMap)
        self.assertEqual(chrom.length, cm.sequence_length)

    def test_known_chromosome(self):
        species = stdpopsim.get_species("CanFam")
        genetic_map = species.get_genetic_map("Campbell2016_CanFam3_1")
        chrom = species.genome.get_chromosome("1")
        cm = genetic_map.get_chromosome_map(chrom.id)
        self.assertIsInstance(cm, msprime.RateMap)
        self.assertEqual(chrom.length, cm.sequence_length)

    def test_warning_for_long_recomb_map(self):
        species = stdpopsim.get_species("HomSap")
        genetic_map = species.get_genetic_map("HapMapII_GRCh37")
        chrom = species.genome.get_chromosome("chr1")
        with self.assertWarns(Warning):
            cm = genetic_map.get_chromosome_map(chrom.id)
        self.assertIsInstance(cm, msprime.RateMap)
        self.assertLess(chrom.length, cm.sequence_length)

    def test_unknown_chromosome(self):
        species = stdpopsim.get_species("HomSap")
        genetic_map = species.get_genetic_map("HapMapII_GRCh37")
        for bad_chrom in ["", "ABD", None]:
            with self.assertRaises(ValueError):
                genetic_map.get_chromosome_map(bad_chrom)
