"""
Infrastructure for managing genetic maps.
"""
import pathlib
import tempfile
import tarfile
import logging
import contextlib
import warnings
import os
import urllib.request

import msprime

from . import cache

logger = logging.getLogger(__name__)


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


# TODO change this to use attrs
class GeneticMap:
    """
    Class representing a genetic map for a species. Provides functionality for
    downloading and cacheing recombination maps from a remote URL.

    .. todo: Document the attributes in this class

    :ivar url: The URL where the packed and compressed genetic map can be obtained.
    :vartype url: str
    :ivar file_pattern: The pattern used to map name individual chromosome to
        files, suitable for use with Python's :meth:`str.format` method.
    :vartype file_pattern: str
    """

    def __init__(
            self, species, id=None, url=None, file_pattern=None,
            description=None, long_description=None, citations=None):
        self.id = id
        self.species = species
        self.description = description
        self.long_description = long_description
        self.url = url
        self.file_pattern = file_pattern
        self.description = description
        self.citations = citations

    @property
    def cache_dir(self):
        return pathlib.Path(cache.get_cache_dir()) / "genetic_maps"

    @property
    def species_cache_dir(self):
        return self.cache_dir / self.species.id

    @property
    def map_cache_dir(self):
        return self.species_cache_dir / self.id

    def __str__(self):
        s = "GeneticMap:\n"
        s += "\tspecies   = {}\n".format(self.species.name)
        s += "\tid        = {}\n".format(self.id)
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
            logger.info(f"Clearing cache {self.map_cache_dir}")
            with tempfile.TemporaryDirectory(dir=self.species_cache_dir) as tempdir:
                # Atomically move to a temporary directory, which will be automatically
                # deleted on exit.
                dest = pathlib.Path(tempdir) / "will_be_deleted"
                os.rename(self.map_cache_dir, dest)
        logger.debug(f"Checking species cache directory {self.species_cache_dir}")
        os.makedirs(self.species_cache_dir, exist_ok=True)

        logger.info(f"Downloading genetic map '{self.id}' from {self.url}")
        # os.rename will not work on some Unixes if the source and dest are on
        # different file systems. Keep the tempdir in the same directory as
        # the destination to ensure it's on the same file system.
        with tempfile.TemporaryDirectory(dir=self.species_cache_dir) as tempdir:
            download_file = os.path.join(tempdir, "downloaded")
            extract_dir = os.path.join(tempdir, "extracted")
            urllib.request.urlretrieve(self.url, filename=download_file)
            logger.debug("Extracting genetic map")
            os.makedirs(extract_dir)
            with tarfile.open(download_file, 'r') as tf:
                for info in tf.getmembers():
                    # TODO test for any prefixes on the name; we should just
                    # expand to a normal file. See  the warning here:
                    # https://docs.python.org/3.5/library/tarfile.html#tarfile.TarFile.extractall
                    if not info.isfile():
                        raise ValueError(
                            f"Tarball format error: member {info.name} not a file")
                with cd(extract_dir):
                    tf.extractall()
            # If this has all gone OK up to here we can now move the
            # extracted directory into the cache location. This should
            # minimise the chances of having malformed maps in the cache.
            logger.info("Storing map in {}".format(self.map_cache_dir))
            # os.rename is atomic, and will raise an OSError if the directory
            # already exists. Therefore, if we see the map exists we assume
            # that some other thread has already dowloaded it and raise a
            # warning.
            try:
                os.rename(extract_dir, self.map_cache_dir)
            except OSError:
                warnings.warn(
                    "Error occured renaming map directory. Are several threads/processes"
                    "downloading this map at the same time?")

    def get_chromosome_map(self, id):
        """
        Returns the genetic map for the chromosome with the specified ``id``.

        :param str id: The chromosome identifier.
             A complete list of chromosome IDs for each species can be found in the
             "Genome" subsection for the species in the :ref:`sec_catalog`.
        :rtype: :class:`msprime.RecombinationMap`
        :return: A :class:`msprime.RecombinationMap` object.
        """
        chrom = self.species.genome.get_chromosome(id)
        if not self.is_cached():
            self.download()
        # We assume that if the map file does not exist this is a property of the
        # map itself and not a download error. If a failure occurs reading the map
        # this is propagated to the user, as this indicates a corrupted map which
        # needs to be redownloaded.
        map_file = os.path.join(self.map_cache_dir, self.file_pattern.format(id=id))
        if os.path.exists(map_file):
            ret = msprime.RecombinationMap.read_hapmap(map_file)
        else:
            warnings.warn(
                "Warning: recombination map not found for chromosome: '{}'"
                " on map: '{}', substituting a flat map with chromosome "
                "recombination rate {}".format(
                    id, self.id, chrom.recombination_rate))
            ret = msprime.RecombinationMap.uniform_map(
                    chrom.length, chrom.recombination_rate)
        return ret
