"""
Citation management for stdpopsim. Provides utilities for printing
citation information associated with different entities derived from the
literature that are used within a simulation.
"""
import collections
import urllib.request

import attr


# TODO this should be an enum.Enum


class CiteReason:
    ENGINE = "simulation engine"
    DEM_MODEL = "demographic model"
    GEN_MAP = "genetic map"
    POP_SIZE = "population size"
    GEN_TIME = "generation time"
    MUT_RATE = "mutation rate"
    REC_RATE = "recombination rate"
    GENE_CONVERSION = "gene conversion parameters"
    ASSEMBLY = "genome assembly"
    ANNOTATION = "genome annotation"
    DFE = "distribution of fitness effects"
    STDPOPSIM = "stdpopsim"


@attr.s
class Citation:
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

    def displaystr(self):
        s = ["["]
        for i, reason in enumerate(self.reasons):
            if i != 0:
                s.append(", ")
            s.append(reason)
        s.append("]\n")
        s.append(str(self))
        return "".join(s)

    def because(self, reasons):
        """Returns a new Citation with the given reasons."""
        if not isinstance(reasons, set):
            reasons = {reasons}
        return self.__class__(
            author=self.author,
            year=self.year,
            doi=self.doi,
            reasons=self.reasons | reasons,
        )

    @staticmethod
    def merge(citations):
        """Returns a deduplicated list of Citation objects."""
        cset = collections.OrderedDict()
        for citation in citations:
            if citation.doi in cset:
                citation = cset[citation.doi].because(citation.reasons)
            cset[citation.doi] = citation
        return list(cset.values())

    def assert_valid(self):
        """
        Checks that this citation is valid by checking the types and values
        of the instance variables.
        """
        assert isinstance(self.author, str)
        assert len(self.author) > 0
        assert isinstance(self.doi, str)
        assert len(self.doi) > 0
        parsed = urllib.parse.urlparse(self.doi)
        assert parsed.scheme.startswith("http")
        assert isinstance(self.year, int)
        assert self.year > 0
        for reason in self.reasons:
            assert isinstance(reason, str)
            # TODO check that it's in the set of accepted reasons.

    def fetch_bibtex(self):
        """Retrieve the bibtex of a citation from Crossref."""
        req = urllib.request.Request(self.doi)
        req.add_header("Accept", "text/bibliography; style=bibtex")
        with urllib.request.urlopen(req) as con:
            return con.read().decode()


_stdpopsim_citation = Citation(
    doi="https://doi.org/10.7554/eLife.54967",
    year="2020",
    author="Adrion et al.",
    reasons={CiteReason.STDPOPSIM},
)
