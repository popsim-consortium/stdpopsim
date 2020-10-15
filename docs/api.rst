.. _sec_api:

===
API
===

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

.. autoclass:: stdpopsim.GeneticMap()
    :members:

.. autoclass:: stdpopsim.Citation()
    :members:

.. autoclass:: stdpopsim.Annotation()
    :members:

******************
Demographic Models
******************

.. autoclass:: stdpopsim.DemographicModel()
    :members:

.. autoclass:: stdpopsim.Population()
    :members:

.. _sec_api_generic_models:

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
