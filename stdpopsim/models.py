"""
Common infrastructure for specifying demographic models.
"""
import copy
import textwrap

import msprime
import warnings
import stdpopsim


class Population:
    # TODO deprecate this - we don't need any internal definition of what
    # a population is.
    """
    Class recording metadata representing a population in a simulation.

    :ivar ~.id: The id of the population.
    :vartype ~.id: str
    :ivar ~.description: a short description of the population
    :vartype ~.description: str
    :ivar sampling_time: an integer value which indicates how many
        generations prior to the present individuals should samples should
        be drawn from this population. If `None`, sampling not allowed from this
        population (default = 0).
    :vartype sampling_time: int
    """

    def __init__(self, id, description, sampling_time=0):
        self.id = id
        self.description = description
        self.sampling_time = sampling_time

    @property
    def allow_samples(self):
        return self.sampling_time is not None

    def asdict(self):
        """
        Returns a dictionary representing the metadata about this population.
        """
        return {
            "id": self.id,
            "description": self.description,
            "sampling_time": self.sampling_time,
        }


class DemographicModel:
    """
    Class representing a demographic model.

    Instances of this class are constructed by model implementors, following the
    :ref:`developer documentation <sec_development_demographic_model>`. To instead
    obtain a pre-specified model as listed in the :ref:`sec_catalog`,
    see :class:`Species.get_demographic_model`.

    :ivar ~.id: The unique identifier for this model. DemographicModel IDs should be
        short and memorable, and conform to the stdpopsim
        :ref:`naming conventions <sec_development_naming_conventions>`
        for demographic models.
    :vartype ~.id: str
    :ivar ~.description: A short description of this model as it would be used in
        written text, e.g., "Three population Out-of-Africa". This should
        describe the model itself and not contain author or year information.
    :vartype ~.description: str
    :ivar long_description: A concise, but detailed, summary of the model.
    :vartype long_description: str
    :ivar generation_time: Mean inter-generation interval, in years.
    :vartype generation_time: int
    :ivar citations: A list of :class:`Citation`, that describe the primary
        reference(s) for the model.
    :vartype citations: list of :class:`Citation`
    :ivar mutation_rate: The mutation rate associated with the demographic model.
        If no mutation rate is given, the species' default mutation rate is used.
    :vartype mutation_rate: float
    """

    def __init__(
        self,
        *,
        id,
        description,
        long_description,
        generation_time=None,
        citations=None,
        qc_model=None,
        population_configurations=None,
        demographic_events=None,
        migration_matrix=None,
        population_id_map=None,
        populations=None,
        model=None,
        mutation_rate=None,
    ):
        self.id = id
        self.description = description
        self.long_description = long_description
        self.generation_time = 1 if generation_time is None else generation_time
        self.mutation_rate = mutation_rate
        self.citations = [] if citations is None else citations
        self.qc_model = qc_model

        if model is None:
            assert population_configurations is not None
            assert populations is not None
            assert len(populations) == len(population_configurations)
            population_configurations = copy.deepcopy(population_configurations)
            # Merge the information from the populations into the msprime
            # Demography.
            for pop, pop_config in zip(populations, population_configurations):
                if pop_config.metadata is None:
                    pop_config.metadata = {}
                if pop.id != "" and not pop.id.startswith("qc_"):
                    pop_config.metadata["name"] = pop.id
                    pop_config.metadata["description"] = pop.description
            # This will become a Demes model in the future - for now it's an
            # msprime model.
            self.model = msprime.Demography.from_old_style(
                population_configurations=population_configurations,
                demographic_events=demographic_events,
                migration_matrix=migration_matrix,
                population_map=population_id_map,
            )
            for msp_pop, local_pop in zip(self.model.populations, populations):
                # We use the "allow_samples" attribute in the CLI and else where
                # so we monkey patch this into the msprime Populations for the
                # moment.
                msp_pop.allow_samples = True
                if local_pop.sampling_time is not None:
                    msp_pop.default_sampling_time = local_pop.sampling_time
                else:
                    msp_pop.allow_samples = False
        else:
            assert population_configurations is None
            assert population_id_map is None
            assert populations is None
            self.model = copy.deepcopy(model)
            # See note above. We allow samples in all populations for now.
            for pop in self.model.populations:
                pop.allow_samples = True

    def __str__(self):
        long_desc_lines = [
            line.strip()
            for line in textwrap.wrap(textwrap.dedent(self.long_description))
        ]
        long_desc = "\n║                     ".join(long_desc_lines)
        s = (
            "Demographic model:\n"
            f"║  id               = {self.id}\n"
            f"║  description      = {self.description}\n"
            f"║  long_description = {long_desc}\n"
            f"║  generation_time  = {self.generation_time}\n"
            f"║  mutation_rate    = {self.mutation_rate}\n"
            f"║  citations        = {[cite.doi for cite in self.citations]}\n"
            f"║{self.model}"
        )
        return s

    @property
    def populations(self):
        return self.model.populations

    @property
    def num_populations(self):
        return self.model.num_populations

    @property
    def num_sampling_populations(self):
        return sum(int(pop.allow_samples) for pop in self.populations)

    def register_qc(self, qc_model):
        """
        Register a QC model implementation for this model.
        """
        if not isinstance(qc_model, self.__class__):
            raise ValueError(
                f"Cannot register non-DemographicModel '{qc_model}' as QC model"
            )
        if self.qc_model is not None:
            raise ValueError(f"QC model already registered for {self.id}.")
        self.qc_model = qc_model

    def get_sample_sets(self, individuals_per_population, ploidy=2):
        """
        Returns a list of `msprime.SampleSet` objects, using the number of
        individuals sampled from each population as specified in
        `individuals_per_population`.  For instance,
        ``model.get_sample_sets({"pop_0":5, "pop_1":5}`` would return sample
        sets encompassing 10 individuals, equally split between populations
        named "pop_0" and "pop_1".  Omitted populations are assumed to have
        zero samples.

        :param dict individuals_per_population: A dict of the form
            `{population_name:number_of_samples}` giving sample counts per
            population.
        :param int ploidy: The ploidy of the samples.

        .. note: This interface is intended for internal use only.
        """

        sample_sets = []
        pop_map = {pop.name: pop for pop in self.populations}
        pop_idx = {pop.name: i for i, pop in enumerate(self.populations)}
        for name, num_samples in individuals_per_population.items():
            if name not in pop_map:
                raise ValueError(
                    f"Population {name} is not one of the populations "
                    f"{list(pop_map.keys())} defined for model {self.id}"
                )
            pop = pop_map[name]
            if pop.allow_samples:
                sample_sets.append(
                    msprime.SampleSet(
                        num_samples=num_samples,
                        population=pop_idx[name],
                        time=pop.default_sampling_time,
                        ploidy=ploidy,
                    )
                )
            elif num_samples > 0:
                raise ValueError(
                    f"Individuals requested from non-sampling population {name}"
                )
        sample_sets = sorted(sample_sets, key=lambda x: x.population)
        return sample_sets

    def get_samples(self, *args):
        """
        Returns a list of `msprime.SampleSet` objects, with the number of
        haploids from each population determined by the positional arguments.
        For instance, ``model.get_samples(2, 5, 7)`` would return a list with
        sample sets encompassing 14 samples, two of which are from the model's
        first population (i.e., with population ID ``model.populations[0].id``),
        five are from the model's second population, and seven are from the
        model's third population.  The number of of arguments must be less than
        or equal to the number of "sampling" populations,
        ``model.num_sampling_populations``; if the number of arguments is less
        than the number of sampling populations, then remaining numbers are
        treated as zero.

        .. note: This interface is deprecated. Instead, a dict containing the
            number of individuals per population should be directly provided to
            :meth:`Engine.simulate`.
        """

        warnings.warn(
            stdpopsim.DeprecatedFeatureWarning(
                "The use of `DemographicModel.get_samples` (Python API) and "
                "positional sample counts (CLI) is deprecated. Instead, supply a "
                "{population_name:num_samples} dict to "
                "`Engine.simulate(samples=...)` (Python API); or use the syntax "
                "`stdpopsim SpeciesName population_name:num_samples` (CLI)."
            )
        )

        samples = []

        for pop_index, n in enumerate(args):
            if self.populations[pop_index].allow_samples:
                samples.append(
                    msprime.SampleSet(
                        num_samples=n,
                        population=pop_index,
                        time=self.populations[pop_index].default_sampling_time,
                        ploidy=1,
                    )
                )
            elif n > 0:
                pop_name = self.populations[pop_index].name
                raise ValueError(
                    f"Samples requested from non-sampling population {pop_name}"
                )
        return samples


class PiecewiseConstantSize(DemographicModel):
    """
    Class representing a generic simulation model that can be run to output a
    tree sequence. This is a piecewise constant size model, which allows for
    instantaneous population size change over multiple epochs in a single population.

    :param float N0: The initial effective population size
    :param args: Each subsequent argument is a tuple (t, N) which gives the
        time at which the size change takes place and the population size.

    The usage is best illustrated by an example:

    .. code-block:: python

        model1 = stdpopsim.PiecewiseConstantSize(N0, (t1, N1))  # One change
        model2 = stdpopsim.PiecewiseConstantSize(N0, (t1, N1), (t2, N2))  # Two changes
    """

    def __init__(self, N0, *args):
        model = msprime.Demography.isolated_model(initial_size=[N0])
        for t, N in args:
            model.add_population_parameters_change(time=t, initial_size=N)

        super().__init__(
            id="PiecewiseConstant",
            description="Piecewise constant size population model over multiple epochs.",
            long_description="",
            model=model,
            generation_time=1,
        )


class IsolationWithMigration(DemographicModel):
    """
    Class representing a generic simulation model that can be run to output a tree
    sequence. A generic isolation with migration model where a single ancestral
    population of size NA splits into two populations of constant size N1
    and N2 time T generations ago, with migration rates M12 and M21 between
    the split populations. Sampling is disallowed in population index 0,
    as this is the ancestral population.

    :param float NA: The initial ancestral effective population size
    :param float N1: The effective population size of population 1
    :param float N2: The effective population size of population 2
    :param float T: Time of split between populations 1 and 2 (in generations)
    :param float M12: Migration rate from population 1 to 2
    :param float M21: Migration rate from population 2 to 1

    Example usage:

    .. code-block:: python

        model1 = stdpopsim.IsolationWithMigration(NA, N1, N2, T, M12, M21)

    """

    def __init__(self, NA, N1, N2, T, M12, M21):
        model = msprime.Demography()
        model.add_population(initial_size=N1, name="pop1")
        model.add_population(initial_size=N2, name="pop2")
        model.add_population(initial_size=NA, name="ancestral")

        # FIXME This is BACKWARDS in time, so the rates are the other
        # way around forwards time. We should explain this in the documentation
        # (and probably swap around). Seems like there's not really much
        # good reason to have this model in here any more though - what
        # does it do that wouldn't be better done in demes/msprime?
        model.set_migration_rate(source="pop1", dest="pop2", rate=M12)
        model.set_migration_rate(source="pop2", dest="pop1", rate=M21)
        model.add_population_split(
            time=T, ancestral="ancestral", derived=["pop1", "pop2"]
        )
        long_description = """
            A generic isolation with migration model where a single ancestral
            population of size NA splits into two populations of constant size N1
            and N2 time T generations ago, with migration rates M12 and M21 between
            the split populations.
            """
        super().__init__(
            id="IsolationWithMigration",
            description="Generic IM model",
            long_description=long_description,
            model=model,
            generation_time=1,
        )
