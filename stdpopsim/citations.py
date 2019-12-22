"""
Citation management for stdpopsim. Provides utilities for printing
citation information associated with different entities derived from the
literature that are used within a simulation.
"""
import sys
import collections
import urllib.request

import attr


class CiteReason(object):
    ENGINE = "simulation engine"
    DEM_MODEL = "demographic model"
    GEN_MAP = "genetic map"
    POP_SIZE = "population size"
    GEN_TIME = "generation time"
    MUT_RATE = "mutation rate"
    REC_RATE = "recombination rate"
    ASSEMBLY = "genome assembly"


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
    reasons = attr.ib(factory=set, kw_only=True)

    def __str__(self):
        return f"{self.author}, {self.year}: {self.doi}"

    def print(self, file=sys.stdout):
        print("[", end="", file=file)
        for i, reason in enumerate(self.reasons):
            end = "" if i == len(self.reasons)-1 else ", "
            print(reason, end=end, file=file)
        print("]", file=file)
        print("  "+str(self), file=file)

    def because(self, reasons):
        """Returns a new Citation with the given reasons."""
        if not isinstance(reasons, set):
            reasons = {reasons}
        return self.__class__(
                author=self.author, year=self.year, doi=self.doi,
                reasons=self.reasons | reasons)

    @staticmethod
    def merge(citations):
        """Returns a deduplicated list of Citation objects."""
        cset = collections.OrderedDict()
        for citation in citations:
            if citation.doi in cset:
                citation = cset[citation.doi].because(citation.reasons)
            cset[citation.doi] = citation
        return list(cset.values())

    def fetch_bibtex(self):
        """Retrieve the bibtex of a citation from Crossref."""
        req = urllib.request.Request(self.doi)
        req.add_header("Accept", "text/bibliography; style=bibtex")
        with urllib.request.urlopen(req) as con:
            return con.read().decode()
