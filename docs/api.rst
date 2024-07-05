.. _sec_api:

===
API
===

***************
Quick reference
***************

Functions to get things:

.. autosummary::

    stdpopsim.get_species
    stdpopsim.Species.get_contig
    stdpopsim.Species.get_demographic_model
    stdpopsim.Species.get_dfe
    stdpopsim.Species.get_annotations
    stdpopsim.Contig.add_dfe
    stdpopsim.get_engine

Classes of objects:

.. autosummary::

   stdpopsim.Species
   stdpopsim.Genome
   stdpopsim.Chromosome
   stdpopsim.Contig
   stdpopsim.Citation
   stdpopsim.Annotation
   stdpopsim.DFE
   stdpopsim.MutationType
   stdpopsim.DemographicModel
   stdpopsim.Population

.. _sec_api_species_definitions:

*******************
Species definitions
*******************

The :ref:`Catalog <sec_catalog>` contains a large number of species and simulation
model definitions, which are built using a number of classes defined here.
These are usually not intended to be instantiated directly, but should be
accessed through the main entry point, :func:`.get_species`.

.. autofunction:: stdpopsim.get_species

.. autoclass:: stdpopsim.Species()
    :members:

.. autoclass:: stdpopsim.Genome()
    :members:

.. autoclass:: stdpopsim.Chromosome()
    :members:

.. autoclass:: stdpopsim.Contig()
    :members:

.. autoclass:: stdpopsim.Citation()
    :members:

.. autoclass:: stdpopsim.Annotation()
    :members:

.. _sec_api_demographic_models:

******************
Demographic Models
******************

.. autoclass:: stdpopsim.DemographicModel()
    :members:

.. autoclass:: stdpopsim.Population()
    :members:


.. _sec_api_gene_conversion:

*******************************************
Gene conversion and bacterial recombination
*******************************************

Some species have estimates of the process of gene conversion.
However, by default gene conversion does *not* happen, because
(unfortunately) it can make simulations take much longer,
and result in much bigger files.
This is because the rate of gene conversion is usually much higher than
the rate of crossing-over, so enabling gene conversion effectively
increases the rate of recombination, which can have a strong effect on runtimes.
To enable gene conversion for a species that has estimates,
pass the `use_species_gene_conversion=True` argument to
:meth:`Engine.simulate`.

Gene conversion attributes are visible through the `gene_conversion_fraction`
and `gene_conversion_length` properties of the species' `genome`.
The "gene conversion fraction" gives the (average) fraction of
double-stranded breaks that resolve as gene conversion
(a.k.a non-crossover) events,
and the "gene conversion length" is the mean length of the gene conversion tracts.
So, if the rate of crossing over (which is often referred to as
"recombination rate") is `r` and the gene conversion fraction is `f`,
then the *total* rate of double-stranded breaks is `r / (1-f)`,
and the gene conversion rate is `r * f / (1-f)`.

A consequence of this is that
the `recombination_map` attribute of a :class:`Contig` is the rate of
double-stranded break initiation, not the rate of crossovers. The terminology
conflicts somewhat with the `recombination_rate` property of
a :class:`Chromosome`, which specifies the rate of crossovers. So, creating
a contig with `use_species_gene_conversion=True` will result in a Contig
with larger "recombination rates" than otherwise, because these rates
include gene conversion as well as crossing over:

.. code-block:: python

    species = stdpopsim.get_species("DroMel")
    contig = species.get_contig("2L")
    contig.recombination_map.mean_rate
    # 2.40462600791e-08
    contig = species.get_contig("2L", use_species_gene_conversion=True)
    contig.recombination_map.mean_rate
    # 1.414485887005882e-07


Different engines implement gene conversion in different ways: at the time
of writing, msprime only allows a *constant* rate of gene conversion,
while SLiM only allows gene conversion to be a fixed fraction of double-stranded
break events. So, in SLiM we use the recombination map for the local
rate of double-stranded breaks along the genome, while in msprime
we specify a constant rate of gene conversion so that the average total number
is the same as we would get from SLiM.

In principle, **bacterial recombination** is a completely different issue,
since bacterial recombination is by horizontal transfer of genomic segments.
However, bacterial recombination is at present implemented using the engine's
gene conversion mechanisms. So, bacterial species have the `bacterial_recombination`
flag set, and need `gene_conversion_length` to be defined. (It could also be called horizontal
transfer or homologous recombination segment length, but that would add another
option.) For such species, gene conversion does not happen, so `gene_conversion_fraction`
must not be set, and the "recombination rate" is the rate of bacterial recombination.


.. _sec_api_dfes:

******************
Selection and DFEs
******************

To allow for different kinds of mutations,
each :class:`Contig` carries a list of :class:`.DFE` objects
(for "Distribution of Fitness Effects"),
and a list of collections of intervals to which these apply.
The :class:`Contig` contains the overall mutation rate,
and the :class:`.DFE` describes the proportions of different types of mutations.
Each :class:`.Contig` comes, by default, with a single, neutral DFE
that applies to the entire Contig.


.. autoclass:: stdpopsim.DFE()
    :members:

.. autoclass:: stdpopsim.MutationType()
    :members:

.. _sec_api_generic_models:


.. _sec_api_sweeps:

********************
Selection and sweeps
********************

:class:`ExtendedEvent` and subclasses may be used to condition on sequences of
events at particular loci using the SLiM engine, by passing lists of events to
the `extended_events` argument in `Engine.simulate`. A simplified API is provided
to construct the necessary events for selective sweeps.

.. autoclass:: stdpopsim.DrawMutation()
    :members:

.. autoclass:: stdpopsim.ChangeMutationFitness()
    :members:

.. autoclass:: stdpopsim.ConditionOnAlleleFrequency()
    :members:

.. autofunction:: stdpopsim.selective_sweep


**************
Generic models
**************

The :ref:`Catalog <sec_catalog>` contains simulation models from the literature
that are defined for particular species. It is also useful to be able
to simulate more generic models, which are documented here.
Please see the :ref:`sec_tutorial_generic_models` for examples of using
these models.

.. autoclass:: stdpopsim.PiecewiseConstantSize

.. autoclass:: stdpopsim.IsolationWithMigration


.. _sec_api_simulation_engines:

******************
Simulation Engines
******************

Support for additional simulation engines can be implemented by subclassing
the abstract :class:`.Engine` class, and registering an instance of the
subclass with :func:`.register_engine`.
These are usually not intended to be instantiated directly, but should be
accessed through the main entry point, :func:`.get_engine`.

.. autofunction:: stdpopsim.get_engine

.. autofunction:: stdpopsim.get_default_engine

.. autofunction:: stdpopsim.register_engine

.. autoclass:: stdpopsim.Engine()
    :members:

.. autoclass:: stdpopsim.engines._MsprimeEngine()
    :show-inheritance:
    :members: id, description, simulate

.. autoclass:: stdpopsim.slim_engine._SLiMEngine()
    :show-inheritance:
    :members: id, description, simulate, recap_and_rescale
