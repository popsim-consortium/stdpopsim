"""
Infrastructure for defining information about genome annotation.
"""
import logging
import attr
import pandas
import pathlib
import warnings
import os
import urllib.request
from . import cache
import zarr
import tempfile

TEMP="this is a test"

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


@attr.s
class Annotation(object):
    """
    Class represnting a Annotation file
    assume GFF3/GTF or similar

    :ivar url: The URL where the packed and compressed GTF can be found
    :vartype url: str
    :ivar species_id: species id
    :vartype id: str
    :ivar species: a `stdpopsim.species` instance
    :ivar annotation_description: description of annotation file
    :vartype annotation_description: str
    """
    url = attr.ib(default=None)
    zarr_url = attr.ib(default=None)
    species = attr.ib(default=None)
    id = attr.ib(default=None)
    file_name = attr.ib(default=None)
    description = attr.ib(default=None)
    citations = attr.ib(factory=list)
    long_description = attr.ib(default=None)

    def __attrs_post_init__(self):
        self.file_name = os.path.basename(self.zarr_url)

    @property
    def annot_cache_dir(self):
        return pathlib.Path(cache.get_cache_dir()) / "annotations"

    @property
    def species_cache_dir(self):
        return self.annot_cache_dir / self.species.id

    def __str__(self):
        s = "GTF Annotation:\n"
        s += "\tspecies   = {}\n".format(self.species.name)
        s += "\tid        = {}\n".format(self.id)
        s += "\turl       = {}\n".format(self.url)
        s += "\tzarr url       = {}\n".format(self.zarr_url)
        s += "\tcached    = {}\n".format(self.is_cached())
        s += "\tcache_dir = {}\n".format(self.species_cache_dir)
        return s

    def is_cached(self):
        """
        Returns True if this annotation is cached locally.
        """
        return os.path.exists(self.species_cache_dir)

    def download(self):
        """
        Downloads the zarr from the source URL and stores it in the
        cache directory. If the annotation directory already exists it is first
        removed.
        """
        self.file_name = os.path.basename(self.zarr_url)
        if self.is_cached():
            logger.info(f"Clearing cache {self.species_cache_dir}")
            with tempfile.TemporaryDirectory() as tempdir:
                dest = pathlib.Path(tempdir) / "will_be_deleted"
                os.rename(self.annot_cache_dir, dest)
        logger.debug(f"Checking species cache directory {self.species_cache_dir}")
        os.makedirs(self.species_cache_dir, exist_ok=True)
        download_file = f'{self.species_cache_dir}/{self.file_name}'
        logger.info(f"Downloading Zarr file '{self.id}' from {self.zarr_url}")
        logger.info(f"download_file: {download_file}")
        logger.info(f"species_cache_dir: {self.species_cache_dir}")
        if os.path.exists(download_file):
            warnings.warn("multiple downloads?")
        try:
            urllib.request.urlretrieve(self.zarr_url, filename=download_file)
        except urllib.error.URLError:
            print(f"could not connect to {self.zarr_url}")
            raise
        logger.debug("Download Zarr complete")
        logger.info(f"Storing Zarr in {self.species_cache_dir}")

    def get_chromosome_annotations(self, id):
        """
        Returns the pandas dataframe for
        the chromosome with the specified id.
        """
        if not self.is_cached():
            self.download()
        annot_file = os.path.join(self.species_cache_dir, self.file_name)
        if id is None:
            raise ValueError("bad chrom id")
        chr_prefix = "chr"  # building this in for future generalization
        if id.startswith(chr_prefix):
            id = id[len(chr_prefix):]
        if os.path.exists(annot_file):
            bed = zarr_to_dataframe(annot_file)
            assert type(bed) == pandas.DataFrame
            ret = bed[bed.seqid == id]
            if len(ret) == 0:
                raise ValueError
        else:
            ret = None
            raise ValueError(
                "Warning: annotation file not found for chromosome: '{}'"
                " on annotation: '{}', no annotations will be used".format(
                    id, self.id))
        return ret

    def get_annotation_type_from_chromomosome(self, a_type, chrom_id, full_table=False):
        """
        Returns all elements of
        type a_type from chromosome
        specified
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
        Returns all elements of
        type gene from annotation
        """
        return self.get_annotation_type_from_chromomosome('gene', chrom_id,
                                                          full_table)
