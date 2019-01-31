"""
Infrastructure for defining chromosome information for different species.
"""
import os.path
import tempfile
import tarfile
import logging
import contextlib
import shutil

import appdirs
import requests
import msprime


logger = logging.getLogger(__name__)


genetic_maps = {}


def get_genetic_map(species, name):
    """
    Returns the genetic map with the specified name for the specified species.
    Raises a ValueError if the map has not been registered.
    """
    key = "{}/{}".format(species, name)
    if key not in genetic_maps:
        raise ValueError("Unknown genetic map '{}'".format(key))
    return genetic_maps[key]


def register_genetic_map(genetic_map):
    """
    Registers the specified recombination map so that it can be loaded on demand.
    """
    key = "{}/{}".format(genetic_map.species, genetic_map.name)
    logger.debug("Registering genetic map '{}'".format(key))
    genetic_maps[key] = genetic_map


@contextlib.contextmanager
def cd(path):
    """
    Convenience function to change the current working directory in a context
    manager.
    """
    old_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old_dir)


class GeneticMap(object):
    """
    Class representing a genetic map for a species. Provides functionality for
    downloading and cacheing recombination maps from a remote URL.
    """
    species = None
    name = None
    url = None
    file_pattern = None

    def __init__(self):
        self.cache_dir = appdirs.user_cache_dir("stdpopsim", "popgensims")
        self.species_cache_dir = os.path.join(self.cache_dir, self.species)
        self.map_cache_dir = os.path.join(self.species_cache_dir, self.name)

    def is_cached(self):
        """
        Returns True if this map is cached locally.
        """
        return os.path.exists(self.map_cache_dir)

    def download(self):
        """
        Downloads this genetic map from the source URL and stores it in the
        cache directory. If the map directory already exists it is first
        removed.
        """
        if self.is_cached():
            logger.info("Clearing cache {}".format(self.map_cache_dir))
            shutil.rmtree(self.map_cache_dir)
        logger.debug("Making species cache directory {}".format(self.species_cache_dir))
        os.makedirs(self.species_cache_dir, exist_ok=True)

        logger.debug("Downloading genetic map '{}' from {}".format(self.name, self.url))
        r = requests.get(self.url, stream=True)
        with tempfile.TemporaryDirectory() as tempdir:
            download_file = os.path.join(tempdir, "downloaded")
            extract_dir = os.path.join(tempdir, "extracted")
            with open(download_file, 'wb') as f:
                for chunk in r.iter_content(chunk_size=1024):
                    f.write(chunk)
            logger.debug("Extracting genetic map")
            os.makedirs(extract_dir)
            with tarfile.open(download_file, 'r') as tf:
                with cd(extract_dir):
                    tf.extractall()
            # If this has all gone OK up to here we can now move the
            # extracted directory into the cache location. This should
            # minimise the chances of having malformed maps in the cache.
            logger.info("Storing map in {}".format(self.map_cache_dir))
            shutil.move(extract_dir, self.map_cache_dir)

    def get_chromosome_map(self, name):
        """
        Returns the genetic map for the chromosome with the specified name.
        """
        # TODO look up a list of known names to give a good error message.
        if not self.is_cached():
            self.download()
        map_file = os.path.join(self.map_cache_dir, self.file_pattern.format(name=name))
        return msprime.RecombinationMap.read_hapmap(map_file)


class Chromosome(object):

    def __init__(self, name, length, mean_recombination_rate, mean_mutation_rate):
        self.name = name
        self.length = length
        self.mean_recombination_rate = mean_recombination_rate
        self.mean_mutation_rate = mean_mutation_rate

    def _get_recombination_map(self, species, map_name, chr_name):
        genetic_map = get_genetic_map(species, map_name)
        return genetic_map.get_chromosome_map(chr_name)
