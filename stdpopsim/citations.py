"""
Citation management for stdpopsim. Provides utilities for printing
citation information associated with different entities derived from the
literature that are used within a simulation.
"""


class CitableMixin(object):
    """
    A simple class that provides class attributes which define citation
    information about an particular entity that should be acknowledged
    by citation. This is intendted to be used as a simple mixin class;
    it does not store or change any state of inheriting objects.
    """

    doi = None
    """
    The DOI for the publication providing the definitive reference.
    """

    author = None
    """
    Short author list.
    """

    year = None
    """
    Year of publication as a 4 digit integer, e.g. 2008.
    """
