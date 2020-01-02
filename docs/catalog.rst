.. _sec_catalog:

=======
Catalog
=======

With ``stdpopsim``, you can run simulations from a number of demographic models
that were implemented from published demographic histories. These models have been
rigorously checked and tested by multiple people, so you can rest easy knowing that
your simulations are reproducible and bug-free!

This catalog shows you all of the possible options that you can use to configure
your simulation.
It is organised around a number of choices that you'll need to make about the
:class:`.Species` you wish to simulate:

1. Which **chromosome**? (ie. which :class:`.Genome` object?)
2. Which **genetic map**? (ie. which :class:`.GeneticMap` object?)
3. Which **model of demographic history**? (ie. which :class:`.DemographicModel` object)

For instance, suppose you are interested in simulating modern human samples of

1. chromosome 22, using
2. the HapMapII genetic map, under
3. a 3-population Out-of-Africa model.

The following command simulates 2 samples from each of the three populations,
and saves the output to a file called ``test.trees``:

.. code-block:: console

    $ stdpopsim HomSap -c chr22 -o test.trees -g HapMapII_GRCh37 --model OutOfAfrica_3G09 2 2 2

(To learn more about using ``stdpopsim`` via the command-line, read our :ref:`sec_tutorial`
about it.)

Are there other well-known organisms, genetic maps or models that
you'd like to see in ``stdpopsim``? Head to our :ref:`sec_development`
page to learn about the process for adding new items to the catalog.
Then, if you feel ready, make an issue on our
`GitHub page <https://github.com/popgensims/stdpopsim/issues>`_.


.. speciescatalog:: HomSap

.. speciescatalog:: DroMel

.. speciescatalog:: AraTha

.. speciescatalog:: EscCol
