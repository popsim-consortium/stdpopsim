.. _sec_introduction:

============
Introduction
============

.. note:: This documentation is incomplete and under development. If
    you would like to help, please open an issue or pull request at
    `GitHub <https://github.com/popgensims/stdpopsim>`_.

This is the documentation for ``stdpopsim``, the standard library for population
genetic simulation models.

We designed ``stdpopsim`` to make it easier for you to run reproducible, bug-free
simulations of genetic datasets from published demographic histories.
Under the hood, ``stdpopsim`` relies on
`msprime <https://msprime.readthedocs.io/en/stable/>`_ and
`SLiM 3 <https://messerlab.org/slim/>`_ to generate sample datasets in the
`tree sequence <https://tskit.readthedocs.io/en/latest/>`_ format.


First steps
-----------

 - Head to the :ref:`Installation <sec_installation>` page to get ``stdpopsim`` installed
   on your computer.

 - Skim the :ref:`Catalog <sec_catalog>` to see what simulations are currently supported
   by ``stdpopsim``.

 - Read the :ref:`Tutorials <sec_tutorial>` to see some examples of ``stdpopsim`` in
   action.

Getting involved
----------------

Are there other features, models or organisms that you'd like to see in ``stdpopsim``?
This software is maintained by the PopSim Consortium,
a global community of scientists and developers who are working together to improve
standards for benchmarking and testing in population genetics.
Read the :ref:`Development <sec_development>` page to find out how you can join us.

Citations
---------

If you use ``stdpopsim`` in your work, please cite us:

.. todo::
	Add link to manuscript when it's ready.

Licence and usage
-----------------

``stdpopsim`` is available under the GPLv3 public license.
The terms of this license can be read
`here <https://www.gnu.org/licenses/gpl-3.0.en.html>`_.
