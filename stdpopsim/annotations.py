"""
Infrastructure for defining information about genome annotation.
"""

import logging
import os
import attr
import numpy as np
import stdpopsim

logger = logging.getLogger(__name__)


@attr.s(kw_only=True)
class Annotation:
    """
    Class representing an annotation track.

    :ivar str ~.id: String that uniquely identifies the annotation.
    :ivar species: The species to which this annotation applies.
    :vartype species: :class:`.Species`
    :ivar str url: The URL where the packed and compressed GFF3 can be found.
    :ivar str intervals_url: The URL of the intervals cache of the annotations.
    :ivar str intervals_sha256: The SHA256 checksum of the annotations cache.
    :ivar str ~.description: One line description of the annotation.
    :ivar citations: List of citations for the annotation.
    :vartype citations: list of :class:`.Citation`
    :ivar file_pattern: The pattern used to map individual chromosome id strings
        to files
        :ivar assembly_name: The name of the genome assembly.
    :vartype assembly_name: str    
    :ivar str assembly_name: The name of the assembly the annotation is based on
    """

    id = attr.ib()
    species = attr.ib()
    url = attr.ib()
    gff_sha256 = attr.ib()
    intervals_url = attr.ib()
    intervals_sha256 = attr.ib()
    description = attr.ib()
    citations = attr.ib(factory=list)
    file_pattern = attr.ib()
    annotation_source = attr.ib()
    annotation_type = attr.ib()
    assembly_name = attr.ib()

    def __attrs_post_init__(self):
        self._cache = stdpopsim.CachedData(
            namespace=f"annotations/{self.species.id}/{self.id}",
            url=self.intervals_url,
            sha256=self.intervals_sha256,
            extract=True,
        )
        # logging.info(f"annotation namespace = {self._cache.namespace}")

    @property
    def cache_path(self):
        return self._cache.cache_path

    def __str__(self):
        s = "GTF Annotation:\n"
        s += "\tspecies   = {}\n".format(self.species.name)
        s += "\tid        = {}\n".format(self.id)
        s += "\turl       = {}\n".format(self.url)
        s += "\tintervals url  = {}\n".format(self.intervals_url)
        s += "\tcached    = {}\n".format(self.is_cached())
        s += "\tcache_path = {}\n".format(self.cache_path)
        return s

    def is_cached(self):
        """
        Returns True if this annotation is cached locally.
        """
        return self._cache.is_valid()

    def download(self):
        """
        Downloads the intervals URL and stores it in the cache directory.
        """
        self._cache.download()

    def get_chromosome_annotations(self, id):
        """
        Returns the numpy interval array for the chromosome with the specified id.
        """
        chrom = self.species.genome.get_chromosome(id)
        if not self.is_cached():
            self.download()
        file_path = os.path.join(self.cache_path, self.file_pattern.format(id=chrom.id))
        ret = np.loadtxt(file_path, dtype="int32")
        if len(ret) == 0:
            raise ValueError(f"No annotations found for {id}")
        return ret
