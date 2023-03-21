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
4. Which **distribution of fitness effects** (ie. which :class:`.DFE` object) within
5. which **annotation track**? (ie. which :class:`.Annotation` object)

For instance, suppose you are interested in simulating modern human samples of

1. chromosome 22, using
2. the HapMapII genetic map, under
3. a 3-population Out-of-Africa model, with
4. background selection acting within
5. exons from the Ensembl Havana 104 annotations.

The following command simulates 2 samples from each of the three populations,
(named YRI, CEU, and CHB) and saves the output to a file called ``test.trees``:

.. code-block:: console

    $ stdpopsim -e slim HomSap -c chr22 -o test.trees -g HapMapII_GRCh38 \
    $    --dfe Gamma_K17 --dfe-annotation ensembl_havana_104_exons \
    $    -d OutOfAfrica_3G09 YRI:2 CEU:2 CHB:2


(To learn more about using ``stdpopsim`` via the command-line, read our
:ref:`tutorial <sec_cli_tute>` about it.)

Are there other well-known organisms, genetic maps or models that
you'd like to see in ``stdpopsim``? Head to our :ref:`sec_development`
page to learn about the process for adding new items to the catalog.
Then, if you feel ready, make an issue on our
`GitHub page <https://github.com/popgensims/stdpopsim/issues>`_.

.. speciescatalog:: AedAeg

.. speciescatalog:: AnaPla

.. speciescatalog:: AnoCar

.. speciescatalog:: AnoGam

.. speciescatalog:: ApiMel

.. speciescatalog:: AraTha

.. speciescatalog:: BosTau

.. speciescatalog:: CaeEle

.. speciescatalog:: CanFam

.. speciescatalog:: ChlRei

.. speciescatalog:: DroMel

.. speciescatalog:: DroSec

.. speciescatalog:: EscCol

.. speciescatalog:: GasAcu

.. speciescatalog:: HelAnn

.. speciescatalog:: HelMel

.. speciescatalog:: HomSap

.. speciescatalog:: MusMus

.. speciescatalog:: PanTro

.. speciescatalog:: PapAnu

.. speciescatalog:: PonAbe

.. speciescatalog:: StrAga


Generic models
==============

In addition to the species-specific models listed in this catalog, ``stdpopsim`` offers
a number of generic demographic models that can be run with any species.
These are described in more detail in the :ref:`API <sec_api_generic_models>`.
Simulations using these generic models must be run via the Python interface; see our
:ref:`Python tutorial <sec_python_tute>` to learn how to do this.

 - :meth:`stdpopsim.PiecewiseConstantSize`
 - :meth:`stdpopsim.IsolationWithMigration`
