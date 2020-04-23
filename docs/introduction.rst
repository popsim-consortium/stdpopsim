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

If you use ``stdpopsim`` in your work, please cite our
`manuscript <https://doi.org/10.1101/2019.12.20.885129>`_:

Jeffrey R. Adrion, Christopher B. Cole, Noah Dukler, Jared G. Galloway,
Ariella L. Gladstein, Graham Gower, Christopher C. Kyriazis, Aaron P. Ragsdale,
Georgia Tsambos, Franz Baumdicker, Jedidiah Carlson, Reed A. Cartwright,
Arun Durvasula, Bernard Y. Kim, Patrick McKenzie, Philipp W. Messer,
Ekaterina Noskova, Diego Ortega-Del Vecchyo, Fernando Racimo, Travis J. Struck,
Simon Gravel, Ryan N. Gutenkunst, Kirk E. Lohmeuller, Peter L. Ralph,
Daniel R. Schrider, Adam Siepel, Jerome Kelleher, Andrew D. Kern (2019),
*A community-maintained standard library of population genetic models*,
bioRxiv 2019.12.20.885129; doi: https://doi.org/10.1101/2019.12.20.885129


Bibtex record::

    @article {Adrion2019.12.20.885129,
        author = {Adrion, Jeffrey R. and Cole, Christopher B. and Dukler, Noah
            and Galloway, Jared G. and Gladstein, Ariella L. and Gower, Graham and
            Kyriazis, Christopher C. and Ragsdale, Aaron P. and Tsambos, Georgia and
            Baumdicker, Franz and Carlson, Jedidiah and Cartwright, Reed A. and
            Durvasula, Arun and Kim, Bernard Y. and McKenzie, Patrick and Messer, Philipp W.
            and Noskova, Ekaterina and Vecchyo, Diego Ortega-Del and Racimo, Fernando
            and Struck, Travis J. and Gravel, Simon and Gutenkunst, Ryan N. and
            Lohmeuller, Kirk E. and Ralph, Peter L. and Schrider, Daniel R. and Siepel, Adam
            and Kelleher, Jerome and Kern, Andrew D.},
        title = {A community-maintained standard library of population genetic models},
        elocation-id = {2019.12.20.885129},
        year = {2019},
        doi = {10.1101/2019.12.20.885129},
        publisher = {Cold Spring Harbor Laboratory},
        URL = {https://www.biorxiv.org/content/early/2019/12/21/2019.12.20.885129},
        eprint = {https://www.biorxiv.org/content/early/2019/12/21/2019.12.20.885129.full.pdf},
        journal = {bioRxiv}
    }


Licence and usage
-----------------

``stdpopsim`` is available under the GPLv3 public license.
The terms of this license can be read
`here <https://www.gnu.org/licenses/gpl-3.0.en.html>`_.
