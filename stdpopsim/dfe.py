"""
Common infrastructure for specifying DFEs.
"""
# import numpy as np
import textwrap
import attr
import collections.abc
from . import ext


@attr.s(kw_only=True)
class DFE:
    """
        Class representing a "Distribution of Fitness Effects", i.e., a DFE.

        Instances of this class are constructed by DFE implementors, following the
        :ref:`developer documentation <sec_development_dfe_model>`. To instead
        obtain a pre-specified model as listed in the :ref:`sec_catalog`,
        see :class:`Species.get_dfe`.

    ``proportions`` and ``mutation_types`` must be lists of the same length,
    and ``proportions should be nonnegative numbers summing to 1.

        :ivar ~.mutation_types: a list of MutationTypes associated with the DFE.
        :vartype ~.mutation_types: list
        :ivar ~.proportions: a list of MutationTypes proportions
        :vartype ~.proportions: list
        :ivar ~.id: The unique identifier for this model. DFE IDs should be
            short and memorable, and conform to the stdpopsim
            :ref:`naming conventions <sec_development_naming_conventions>`
            for DFE models.
        :vartype ~.id: str
        :ivar ~.description: A short description of this model as it would be used in
            written text, e.g., "Lognormal DFE". This should
            describe the DFE itself and not contain author or year information.
        :vartype ~.description: str
        :ivar long_description: A concise, but detailed, summary of the DFE model.
        :vartype long_description: str
        :ivar citations: A list of :class:`Citation`, that describe the primary
            reference(s) for the DFE model.
        :vartype citations: list of :class:`Citation`
    """

    id = attr.ib()
    description = attr.ib()
    long_description = attr.ib()
    mutation_types = attr.ib(default=attr.Factory(list))
    proportions = attr.ib(default=attr.Factory(list))
    citations = attr.ib(default=attr.Factory(list))

    def __attrs_post_init__(self):
        self.citations = [] if self.citations is None else self.citations
        if self.proportions == [] and len(self.mutation_types) == 1:
            self.proportions = [1]

        if not (
            isinstance(self.proportions, collections.abc.Collection)
            and isinstance(self.mutation_types, collections.abc.Collection)
            and len(self.proportions) == len(self.mutation_types)
        ):
            raise ValueError(
                "proportions and mutation_types must be lists of the same length."
            )

        for p in self.proportions:
            if not isinstance(p, (float, int)) or p < 0:
                raise ValueError("proportions must be nonnegative numbers.")
        # sum_p = sum(self.proportions)
        #    if not np.isclose(sum_p, 1):
        #       raise ValueError("proportions must sum 1.0.")

        for m in self.mutation_types:
            if not isinstance(m, ext.MutationType):
                raise ValueError(
                    "mutation_types must be a list of MutationType objects."
                )

    def __str__(self):
        long_desc_lines = [
            line.strip()
            for line in textwrap.wrap(textwrap.dedent(self.long_description))
        ]
        long_desc = "\n║                     ".join(long_desc_lines)
        s = (
            "DFE:\n"
            f"║  id               = {self.id}\n"
            f"║  description      = {self.description}\n"
            f"║  long_description = {long_desc}\n"
            f"║  citations        = {[cite.doi for cite in self.citations]}\n"
        )
        return s


def neutral_DFE(convert_to_substitution=True):
    id = "neutral"
    description = "neutral DFE"
    long_description = "strictly neutral mutations"
    neutral = ext.MutationType(convert_to_substitution=convert_to_substitution)
    return DFE(
        id=id,
        description=description,
        long_description=long_description,
        mutation_types=[neutral],
        proportions=[1.0],
    )
