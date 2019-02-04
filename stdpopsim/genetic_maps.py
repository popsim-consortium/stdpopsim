"""
Infrastructure for managing genetic maps.
"""
import os.path
import tempfile
import tarfile
import logging
import contextlib
import shutil
import inspect

import appdirs
import requests
import msprime


logger = logging.getLogger(__name__)


registered_maps = {}


def get_genetic_map(species, name):
    """
    Returns the genetic map with the specified name for the specified species.
    Raises a ValueError if the map has not been registered.
    """
    key = "{}/{}".format(species, name)
    if key not in registered_maps:
        raise ValueError("Unknown genetic map '{}'".format(key))
    return registered_maps[key]


def register_genetic_map(genetic_map):
    """
    Registers the specified recombination map so that it can be loaded on demand.

    A key is generated from each genetic map based on its class and module names,
    giving the name of the map and species, respectively.
    """
    key = "{}/{}".format(genetic_map.species, genetic_map.name)
    logger.debug("Registering genetic map '{}'".format(key))
    registered_maps[key] = genetic_map


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


class classproperty(object):
    """
    Define a 'class property'. Used below for defining GeneticMap name and species.

    https://stackoverflow.com/questions/5189699/how-to-make-a-class-property
    """
    def __init__(self, f):
        self.f = f

    def __get__(self, obj, owner):
        return self.f(owner)


class GeneticMap(object):
    """
    Class representing a genetic map for a species. Provides functionality for
    downloading and cacheing recombination maps from a remote URL.

    Specific genetic maps are defined by subclassing this abstract superclass
    and registering the map.
    """

    url = None
    """
    The URL where this genetic map can be obtained.
    """

    file_pattern = None
    """
    The pattern used to map name individual chromosome to files, suitable for use
    with Python's :meth:`str.format` method.
    """

    def __init__(self):
        self.cache_dir = appdirs.user_cache_dir("stdpopsim", "popgensims")
        self.species_cache_dir = os.path.join(self.cache_dir, self.species)
        self.map_cache_dir = os.path.join(self.species_cache_dir, self.name)

    # We use a bit of trickery here to dynamically get the species name and
    # map ID from the subclass.
    @classproperty
    def name(cls):
        """
        The name of this GeneticMap. This is equal to the name of the class
        defining it.
        """
        return cls.__name__

    @classproperty
    def species(cls):
        """
        The species that this GeneticMap subclass is for. Equal to the name
        of the module in which the class is defined.
        """
        mod = inspect.getmodule(cls).__name__
        species = mod.split(".")[-1]
        return species

    def __str__(self):
        s = "GeneticMap:\n"
        s += "\tspecies   = {}\n".format(self.species)
        s += "\tname      = {}\n".format(self.name)
        s += "\turl       = {}\n".format(self.url)
        s += "\tcached    = {}\n".format(self.is_cached())
        s += "\tcache_dir = {}\n".format(self.map_cache_dir)
        return s

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

        logger.info("Downloading genetic map '{}' from {}".format(self.name, self.url))
        response = requests.get(self.url, stream=True)
        with tempfile.TemporaryDirectory() as tempdir:
            download_file = os.path.join(tempdir, "downloaded")
            extract_dir = os.path.join(tempdir, "extracted")
            with open(download_file, 'wb') as f:
                for chunk in response.iter_content(chunk_size=1024):
                    f.write(chunk)
            logger.debug("Extracting genetic map")
            os.makedirs(extract_dir)
            with tarfile.open(download_file, 'r') as tf:
                for info in tf.getmembers():
                    # TODO test for any prefixes on the name; we should just
                    # expand to a normal file. See  the warning here:
                    # https://docs.python.org/3.5/library/tarfile.html#tarfile.TarFile.extractall
                    if not info.isfile():
                        raise ValueError(
                            "Tarball format error: member {} not a file".format(
                                info.name))
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
