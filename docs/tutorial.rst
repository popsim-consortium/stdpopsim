.. _sec_tutorial:

========
Tutorial
========

.. TODO port these old examples.

.. ********
.. Examples
.. ********

.. This content should go into a tutorial or somewhere else, but for now it's
.. useful to keep a couple of short examples here for reference and
.. to motivate the API.

.. Models and genome data are split by species. Information about the genomes
.. of a particular species is held in the ``genome`` class variable. So,
.. suppose we wish to perform a simple simulation of human chromosome 22, we
.. might have:

.. .. code-block:: python

..     import msprime
..     from stdpopsim import homo_sapiens

..     chrom = homo_sapiens.genome.chromosomes["chr22"]
..     ts = msprime.simulate(
..         sample_size=10,
..         recombination_rate=chrom.default_recombination_rate,
..         mutation_rate=chrom.default_mutation_rate,
..         length=chrom.length)

.. (This should be a very quick simulation, and the result will have very few
.. variants, because although it performs a coalescent simulation of
.. a 51,304,566bp chromosome, it does this with the effective population size of
.. `Ne=1`.)

.. The chromosome definitions also aware of recombination maps,
.. which must be first downloaded. The default for ``homo_sapiens`` is ``HapmapII_GRCh37``,
.. which we can find out, then download it as follows
.. (the maps are stored in your ``~/.cache/stdpopsim/`` directory):

.. .. code-block:: python

..    homo_sapiens.genome.default_genetic_map
..    # 'HapmapII_GRCh37'
..    gmap = homo_sapiens.HapmapII_GRCh37()
..    gmap.download()


.. After this has been done (once only), we can run simulations using this genetic map as follows:

.. .. code-block:: python

..     chrom = homo_sapiens.genome.chromosomes["chr22"]
..     ts = msprime.simulate(
..         sample_size=10,
..         mutation_rate=chrom.default_mutation_rate,
..         recombination_map=chrom.recombination_map())

.. Recombination maps will be downloaded on demand and cached in a
.. platform-appropriate user cache directory (e.g., ``$HOME/.cache/stdpopsim`` on
.. Linux). In this example we didn't specify which recombination map we want, and so the
.. API will use the default. We can also ask for specific maps, if we want:

.. .. code-block:: python

..     chrom = homo_sapiens.genome.chromosomes["chr22"]
..     ts = msprime.simulate(
..         sample_size=10,
..         mutation_rate=chrom.default_mutation_rate,
..         recombination_map=chrom.recombination_map("HapmapII_GRCh37"))


.. Demographic models can also be used. For example

.. .. code-block:: python

..     import stdpopsim
..     from stdpopsim import homo_sapiens

..     chrom = homo_sapiens.genome.chromosomes["chr22"]
..     model = homo_sapiens.GutenkunstThreePopOutOfAfrica()
..     # One sample each from YRI, CEU and CHB.
..     samples = [msprime.Sample(population=j, time=0) for j in range(3)]
..     ts = msprime.simulate(
..         samples=samples,
..         recombination_map=chrom.recombination_map(),
..         mutation_rate=chrom.default_mutation_rate,
..         **model.asdict())

.. (This simulation now has a realistic effective population size,
.. so will produce thousands of variant sties, but still runs very fast.)



.. _sec_tutorial_generic_models:

**************
Generic models
**************

The previous sections showed how to simulate population genetic models that
have been published for particular species. It is also sometimes useful
to simulate more generic models. We can still use the information for
particular species to do this; the only difference is that we instantiate
the models directly rather than retrieving them from the catalog.


.. code-block:: python

    species = stdpopsim.get_species("homsap")
    contig = species.get_contig("chr22", length_multiplier=0.1)
    model = stdpopsim.PiecewiseConstantSize(species.population_size)
    samples = model.get_samples(10)
    ts = model.run(contig, samples)


Here, we simulate 10% of human chromosome 22 under a constant size
population model, using the current best estimate of the human
effective population size from the **TODO add crossref to catalog section**

.. :ref:`sec_catalog_homo_sapiens_genome`




