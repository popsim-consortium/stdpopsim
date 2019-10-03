"""
Citation management for stdpopsim. Provides utilities for printing
citation information associated with different entities derived from the
literature that are used within a simulation.
"""
import attr


@attr.s
class Citation(object):
    """
    A reference to the literature that should be acknowledged by users of
    stdpopsim.

    :ivar doi: The DOI for the publication providing the definitive reference.
    :vartype doi: str
    :ivar author: Short author list, .e.g, "Author 1 et. al".
    :vartype author: str
    :ivar year: Year of publication as a 4 digit integer, e.g. 2008.
    :vartype year: int
    """
    doi = attr.ib(type=str, kw_only=True)
    author = attr.ib(type=str, kw_only=True)
    year = attr.ib(type=int, kw_only=True)

    def __str__(self):
        return f"{self.author}, {self.year}: {self.doi}"
