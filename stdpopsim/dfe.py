"""
Common infrastructure for specifying DFEs.
"""
# import numpy as np
import textwrap
import attr
import collections.abc
import numpy as np


def _copy_converter(x):
    if isinstance(x, list):
        x = x.copy()
    return x


@attr.s(kw_only=True)
class MutationType(object):
    """
    Class representing a "type" of mutation.  The design closely mirrors SLiM's
    MutationType class.

    The main thing that mutation types carry is a way of drawing a selection
    coefficient for each new mutation. This ``distribution_type`` should be one
    of (see the SLiM manual for more information on these):

    - ``f``: fixed, one parameter (the selection coefficient)
    - ``e``: exponential, one parameter (mean)
    - ``g``: gamma, two parameters (mean, shape)
    - ``n``: normal, two parameters (mean, SD)
    - ``w``: Weibull, two parameters (scale, shape)
    - ``u``: Uniform, two parameters (min, max)
    - ``lp``: positive logNormal, two parameters (mean and sd on log scale; see rlnorm)
    - ``ln``: negative logNormal, two parameters (mean and sd on log scale; see rlnorm)

    Type "lp" is always positive, and type "ln" is always negative: both use
    the same log-normal distribution, but "ln" is multiplied by -1.  For
    exponential and gamma, a negative mean can be provided, obtaining always
    negative values.

    Instead of a single dominance coefficient (which would be specified with
    `dominance_coeff`), a discretized relationship between dominance and
    selection coefficient can be implemented: if dominance_coeff_list is
    provided, mutations with selection coefficient s for which
    dominance_coeff_breaks[k-1] <= s <= dominance_coeff_breaks[k] will have
    dominance coefficient dominance_coeff[k]. In other words, the first entry
    of dominance_coeff_list applies to any mutations with selection coefficient
    below the first entry of dominance_coeff_breaks; the second entry of
    dominance_coeff_list applies to mutations with selection coefficient
    between the first and second entries of dominance_coeff_breaks, and so
    forth. The list of breaks must therefore be of length one less than the
    list of dominance coefficients.

    :ivar distribution_type: A one-letter abbreviation for the distribution of
        fitness effects that each new mutation of this type draws from (see below).
    :vartype distribution_type: str
    :ivar distribution_args: Arguments for the distribution type.
    :vartype distribution_type: list
    :ivar dominance_coeff: The dominance coefficient (negative = underdominance,
        0 = recessive, 0.5 = additive, 1.0 = completely dominant, > 1.0 = overdominant)
        Default: 0.5.
    :vartype dominance_coeff: float
    :ivar convert_to_substitution: Whether to retain any fixed mutations in the
        simulation: if not, we cannot ask about their frequency once fixed.
        (Either way, they will remain in the tree sequence).
    :vartype convert_to_substitution: bool
    :ivar dominance_coeff_list: Either None (the default) or a list of floats describing
        a list of dominance coefficients, to apply to different selection coefficients
        (see details). Cannot be specified along with dominance_coeff.
    :vartype dominance_coeff_list: list of floats
    :ivar dominance_coeff_breaks: Either None (the default) or a list of floats
        describing the intervals of selection coefficient over which each of the entries
        of dominance_coeff_list applies (see details). Must be of length one shorter than
        dominance_coeff_list.
    :vartype dominance_coeff_breaks: list of floats
    """

    dominance_coeff = attr.ib(default=None, type=float)
    distribution_type = attr.ib(default="f", type=str)
    distribution_args = attr.ib(
        factory=lambda: [0], type=list, converter=_copy_converter
    )
    convert_to_substitution = attr.ib(default=True, type=bool)
    dominance_coeff_list = attr.ib(default=None, type=list, converter=_copy_converter)
    dominance_coeff_breaks = attr.ib(default=None, type=list, converter=_copy_converter)

    def __attrs_post_init__(self):
        if self.dominance_coeff is None and self.dominance_coeff_list is None:
            self.dominance_coeff = 0.5

        if self.dominance_coeff is not None:
            if (self.dominance_coeff_list is not None) or (
                self.dominance_coeff_breaks is not None
            ):
                raise ValueError(
                    "Cannot specify both dominance_coeff and dominance_coeff_list."
                )
            if not isinstance(self.dominance_coeff, (float, int)):
                raise ValueError("dominance_coeff must be a number.")
            if not np.isfinite(self.dominance_coeff):
                raise ValueError(
                    f"Invalid dominance coefficient {self.dominance_coeff}."
                )

        if self.dominance_coeff_list is not None:
            # disallow the inefficient and annoying length-one case
            if len(self.dominance_coeff_list) < 2:
                raise ValueError("dominance_coeff_list must have at least 2 elements.")
            for h in self.dominance_coeff_list:
                if not isinstance(h, (float, int)):
                    raise ValueError("dominance_coeff_list must be a list of numbers.")
                if not np.isfinite(h):
                    raise ValueError(f"Invalid dominance coefficient {h}.")
            if self.dominance_coeff_breaks is None:
                raise ValueError(
                    "A list of dominance coefficients provided but no breaks."
                )
            if len(self.dominance_coeff_list) != len(self.dominance_coeff_breaks) + 1:
                raise ValueError(
                    "len(dominance_coeff_list) must be equal "
                    "to len(dominance_coeff_breaks) + 1"
                )
            lb = -1 * np.inf
            for b in self.dominance_coeff_breaks:
                if not isinstance(b, (float, int)):
                    raise ValueError(
                        "dominance_coeff_breaks must be a list of numbers."
                    )
                if not np.isfinite(b):
                    raise ValueError(f"Invalid dominance coefficient break {b}.")
                if b < lb:
                    raise ValueError("dominance_coeff_breaks must be nondecreasing.")
                lb = b

        if not isinstance(self.distribution_type, str):
            raise ValueError("distribution_type must be str.")

        if not isinstance(self.distribution_args, list):
            raise ValueError("distribution_args must be list.")

        for i in range(len(self.distribution_args)):
            if not isinstance(self.distribution_args[i], (float, int)):
                raise ValueError(f"distribution_args[{i}] is not a number.")
            if not np.isfinite(self.distribution_args[i]):
                raise ValueError(f"distribution_args[{i}] is an invalid parameter.")

        if not isinstance(self.convert_to_substitution, bool):
            raise ValueError("convert_to_substitution must be bool.")

        # To add a new distribution type: validate the
        # distribution_args here, and add unit tests.
        if self.distribution_type == "f":
            # Fixed-value selection coefficent.
            if len(self.distribution_args) != 1:
                raise ValueError(
                    "Fixed-value mutation types (distribution_type='f') "
                    "take a single selection-coefficient parameter."
                )
        elif self.distribution_type == "g":
            # Gamma-distributed selection coefficient with (mean, shape)
            # parameterisation. A negative value for the mean is permitted,
            # and indicates a reflection of the horizontal axis.
            # See Eidos documentation for rgamma().
            if len(self.distribution_args) != 2:
                raise ValueError(
                    "Gamma-distributed sel. coefs. (distribution_type='g') "
                    "use a (mean, shape) parameterisation."
                )
            if self.distribution_args[1] <= 0:
                raise ValueError("The shape parameter must be positive.")
        elif self.distribution_type == "e":
            # An exponentially-distributed fitness effect (mean).
            # See Eidos documentation for rexp().
            if len(self.distribution_args) != 1:
                raise ValueError(
                    "Exponentially-distributed sel. coefs. (distribution_type='e') "
                    "use a (mean) parameterisation."
                )
        elif self.distribution_type == "n":
            # An normally-distributed fitness effect (mean, standard deviation).
            # See Eidos documentation for rnorm().
            if len(self.distribution_args) != 2:
                raise ValueError(
                    "Normally-distributed sel. coefs. (distribution_type='n') "
                    "use a (mean, sd) parameterisation."
                )
            if self.distribution_args[1] < 0:
                raise ValueError("The sd parameter must be nonnegative.")
        elif self.distribution_type == "w":
            # A Weibull-distributed fitness effect (scale, shape).
            # See Eidos documentation for rweibull().
            if len(self.distribution_args) != 2:
                raise ValueError(
                    "Weibull-distributed sel. coef. (distribution_type='w') "
                    "use a (scale, shape) parameterisation."
                )
            if self.distribution_args[0] <= 0:
                raise ValueError("The scale parameter must be positive.")
            if self.distribution_args[1] <= 0:
                raise ValueError("The shape parameter must be positive.")
        elif self.distribution_type in ("lp", "ln"):
            # A lognormally-distributed fitness effect (meanlog, sdlog),
            # either positive or negative.
            # See Eidos documentation for rlnorm().
            if len(self.distribution_args) != 2:
                raise ValueError(
                    "Lognormally-distributed sel. coefs. (distribution_type='lp'/'ln') "
                    "use a (meanlog, sdlog) parameterisation, requiring sdlog > 0."
                )
            if self.distribution_args[1] < 0:
                raise ValueError("The sdlog parameter must be nonnegative.")
            # dealing with lognormal distribution
            # (adding instead of multiplying the mean):
            logmean = self.distribution_args[0]
            logsd = self.distribution_args[1]
            sign = "" if self.distribution_type == "lp" else "-1 *"
            self.distribution_args = [
                f"return {sign}rlnorm(1, {logmean} + log(Q), {logsd});"
            ]
            self.distribution_type = "s"
        elif self.distribution_type == "u":
            # Uniform
            if (
                len(self.distribution_args) != 2
                or self.distribution_args[0] > self.distribution_args[1]
            ):
                raise ValueError(
                    "Uniformly-distributed sel. coefs. (distribution_type='u') "
                    "use a (min, max) parameterisation, with min <= max."
                )
            umin, umax = self.distribution_args
            self.distribution_args = [f"return runif(1, Q * {umin}, Q * {umax});"]
            self.distribution_type = "s"
        else:
            raise ValueError(
                f"{self.distribution_type} is not a supported distribution type."
            )

        # The index(s) of the param in the distribution_args list that should be
        # multiplied by Q when using --slim-scaling-factor Q.
        self.Q_scaled_index = {
            "e": [0],  # mean
            "f": [0],  # fixed value
            "g": [0],  # mean
            "n": [0, 1],  # mean and sd
            "w": [0],  # scale
            "s": [],  # script types should just printout arguments
        }[self.distribution_type]

    @property
    def is_neutral(self):
        """
        Tests whether the mutation type is strictly neutral. This is defined here to
        be of type "f" and with fitness effect 0.0, and so excludes other situations
        that also produce only neutral mutations (e.g., exponential with mean 0).
        """
        return self.distribution_type == "f" and self.distribution_args[0] == 0


@attr.s(kw_only=True)
class DFE:
    """
    Class representing a "Distribution of Fitness Effects", i.e., a DFE.
    The class records the different *mutation types*, and the *proportions*
    with which they occur. The overall rate of mutations will be determined
    by the Contig to which the DFE is applied (see :meth:`.Contig.add_dfe`).

    Instances of this class are constructed by DFE implementors, following the
    :ref:`developer documentation <sec_development_dfe_model>`. To instead
    obtain a pre-specified model as listed in the :ref:`sec_catalog`,
    see :meth:`Species.get_dfe`.

    ``proportions`` and ``mutation_types`` must be lists of the same length,
    and ``proportions`` should be nonnegative numbers summing to 1.

    :ivar ~.mutation_types: A list of :class:`.MutationType` objects associated
        with the DFE. Defaults to an empty list.
    :vartype ~.mutation_types: list
    :ivar ~.proportions: A list of the proportions of new mutations that
        fall in to each of the mutation types (must sum to 1).
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
    :ivar citations: A list of :class:`Citations <.Citation>`, that describe the primary
        reference(s) for the DFE model.
    :vartype citations: list of :class:`Citation`
    """

    id = attr.ib()
    description = attr.ib()
    long_description = attr.ib()
    mutation_types = attr.ib(default=None)
    proportions = attr.ib(default=None)
    citations = attr.ib(default=None)
    qc_dfe = attr.ib(default=None)

    def __attrs_post_init__(self):
        self.citations = [] if self.citations is None else self.citations
        self.mutation_types = [] if self.mutation_types is None else self.mutation_types
        if self.proportions is None and len(self.mutation_types) == 0:
            self.proportions = []
        elif self.proportions is None:
            # will error below if this doesn't make sense
            self.proportions = [1]

        if not (isinstance(self.proportions, (collections.abc.Sequence, np.ndarray))):
            raise ValueError("proportions must be a list or numpy array.")

        if not (isinstance(self.mutation_types, list)):
            raise ValueError("mutation_types must be a list.")

        if not (len(self.proportions) == len(self.mutation_types)):
            raise ValueError(
                "proportions and mutation_types must be lists of the same length."
            )

        for p in self.proportions:
            if not isinstance(p, (float, int)) or p < 0:
                raise ValueError("proportions must be nonnegative numbers.")

        if len(self.proportions) > 0:
            sum_p = sum(self.proportions)
            if not np.isclose(sum_p, 1):
                raise ValueError("proportions must sum to 1.0.")

        for m in self.mutation_types:
            if not isinstance(m, MutationType):
                raise ValueError(
                    "mutation_types must be a list of MutationType objects."
                )

    @property
    def is_neutral(self):
        return all([m.is_neutral for m in self.mutation_types])

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

    def register_qc(self, qc_dfe):
        """
        Register a QC model implementation for this DFE.
        """
        if not isinstance(qc_dfe, self.__class__):
            raise ValueError(f"Cannot register non-DFE '{qc_dfe}' as QC DFE.")
        if self.qc_dfe is not None:
            raise ValueError(f"QC DFE already registered for {self.id}.")
        self.qc_dfe = qc_dfe


def neutral_dfe(convert_to_substitution=True):
    id = "neutral"
    description = "neutral DFE"
    long_description = "strictly neutral mutations"
    neutral = MutationType(convert_to_substitution=convert_to_substitution)
    return DFE(
        id=id,
        description=description,
        long_description=long_description,
        mutation_types=[neutral],
        proportions=[1.0],
    )
