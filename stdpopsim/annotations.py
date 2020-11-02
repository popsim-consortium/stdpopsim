"""
Infrastructure for defining information about genome annotation.
"""
import logging

import attr
import pandas
import zarr

import stdpopsim

logger = logging.getLogger(__name__)


def zarr_to_dataframe(path):
    """
    converts zarr annotation file to
    pandas dataframe for manipulations
    """
    z = zarr.open(path)
    df = pandas.DataFrame(
        {'seqid': z['seqid'],
         'source': z['source'],
         'type': z['type'],
         'start': z['start'],
         'end': z['end'],
         'score': z['score'],
         'strand': z['strand'],
         'phase': z['phase']}
    )
    return df


@attr.s(kw_only=True)
class Annotation:
    """
    Class representing a GFF3 annotation file.

    :ivar str ~.id: String that uniquely identifies the annotation.
    :ivar species: The species to which this annotation applies.
    :vartype species: :class:`.Species`
    :ivar str url: The URL where the packed and compressed GFF3 can be found.
    :ivar str zarr_url: The URL of the zarr cache of the GFF3.
    :ivar str zarr_sha256: The SHA256 checksum of the zarr cache.
    :ivar str ~.description: One line description of the annotation.
    :ivar citations: List of citations for the annotation.
    :vartype citations: list of :class:`.Citation`
    """
    id = attr.ib()
    species = attr.ib()
    url = attr.ib()
    zarr_url = attr.ib()
    zarr_sha256 = attr.ib()
    description = attr.ib()
    citations = attr.ib(factory=list)

    def __attrs_post_init__(self):
        self._cache = stdpopsim.CachedData(
            namespace=f"annotations/{self.species.id}",
            url=self.zarr_url,
            sha256=self.zarr_sha256,
            extract=False,
        )

    @property
    def cache_path(self):
        return self._cache.cache_path

    def __str__(self):
        s = "GTF Annotation:\n"
        s += "\tspecies   = {}\n".format(self.species.name)
        s += "\tid        = {}\n".format(self.id)
        s += "\turl       = {}\n".format(self.url)
        s += "\tzarr url  = {}\n".format(self.zarr_url)
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
        Downloads the zarr URL and stores it in the cache directory.
        """
        self._cache.download()

    def get_chromosome_annotations(self, id):
        """
        Returns the pandas dataframe for the chromosome with the specified id.
        """
        chrom = self.species.genome.get_chromosome(id)
        if not self.is_cached():
            self.download()
        bed = zarr_to_dataframe(str(self.cache_path))
        assert type(bed) == pandas.DataFrame
        ret = bed[bed.seqid == chrom.id]
        if len(ret) == 0:
            raise ValueError(f"No annotations found for {id}")
        return ret

    def get_annotation_type_from_chromomosome(self, a_type, chrom_id, full_table=False):
        """
        Returns all elements of type a_type from chromosome specified
        """
        annots = self.get_chromosome_annotations(chrom_id)
        subset = annots[annots.type == a_type]
        if subset.empty:
            raise ValueError(f"annotation type '{a_type}' not found on chrom"
                             f" {chrom_id}")
        if full_table:
            return subset
        else:
            return subset[['start', 'end']]

    def get_genes_from_chromosome(self, chrom_id, full_table=False):
        """
        Returns all elements of type gene from annotation
        """
        return self.get_annotation_type_from_chromomosome('gene', chrom_id,
                                                          full_table)
