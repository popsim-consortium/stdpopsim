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
`msprime <https://tskit.dev/msprime/>`_ and
`SLiM 3 <https://messerlab.org/slim/>`_ to generate sample datasets in the
`tree sequence <https://tskit.dev/learn.html#what>`_ format.


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

``stdpopsim`` **is a community effort, and we welcome YOU to join us!**

Are there other features, models or organisms that you'd like to see in ``stdpopsim``?
This software is maintained by the PopSim Consortium,
an open, global community of scientists and developers who are working together to improve
standards for benchmarking and testing in population genetics.
You can get into contact with the ``stdpopsim`` community by subscribing to our `email list
serve <https://lists.uoregon.edu/mailman/listinfo/popgen_benchmark>`_
and by creating and commenting on
Github `issues <http://github.com/popgensims/stdpopsim/issues>`_.
There is a lot of chatter through
Github, and we’ve been building code
there cooperatively. If you are interested in joining us please read our
`code of conduct <https://github.com/popsim-consortium/stdpopsim/blob/main/CODE_OF_CONDUCT.md>`_.
A particular goal of this project is to build a diverse and supportive community,
and so we especially encourage individuals from underrepresented groups in popgen to join.
Being a member of the Consortium can provide a valuable avenue for career development, through
building networks, contributing to papers and code, as well as through learning programing standards
and best practices.

For further details read the :ref:`Development <sec_development>` page to find out how you can join us.

Citations
---------

If you use ``stdpopsim`` in your work, please cite our
`manuscript <https://doi.org/10.7554/eLife.54967>`_:

Jeffrey R Adrion, Christopher B Cole, Noah Dukler, Jared G Galloway,
Ariella L Gladstein, Graham Gower, Christopher C Kyriazis, Aaron P Ragsdale,
Georgia Tsambos, Franz Baumdicker, Jedidiah Carlson, Reed A Cartwright,
Arun Durvasula, Ilan Gronau, Bernard Y Kim, Patrick McKenzie,
Philipp W Messer, Ekaterina Noskova, Diego Ortega-Del Vecchyo, Fernando Racimo,
Travis J Struck, Simon Gravel, Ryan N Gutenkunst, Kirk E Lohmeuller,
Peter L Ralph, Daniel R Schrider, Adam Siepel, Jerome Kelleher, Andrew D Kern (2020),
*A community-maintained standard library of population genetic models*,
eLife 2020;9:e54967; doi: https://doi.org/10.7554/eLife.54967


Bibtex record::

    @article {10.7554/eLife.54967,
        article_type = {journal},
        title = {A community-maintained standard library of population genetic models},
        author = {Adrion, Jeffrey R and Cole, Christopher B and Dukler, Noah and Galloway, Jared
            G and Gladstein, Ariella L and Gower, Graham and Kyriazis, Christopher C and Ragsdale,
            Aaron P and Tsambos, Georgia and Baumdicker, Franz and Carlson, Jedidiah and Cartwright,
            Reed A and Durvasula, Arun and Gronau, Ilan and Kim, Bernard Y and McKenzie, Patrick and
            Messer, Philipp W and Noskova, Ekaterina and Ortega-Del Vecchyo, Diego and Racimo,
            Fernando and Struck, Travis J and Gravel, Simon and Gutenkunst, Ryan N and Lohmueller,
            Kirk E and Ralph, Peter L and Schrider, Daniel R and Siepel, Adam and Kelleher, Jerome
            and Kern, Andrew D},
        editor = {Coop, Graham and Wittkopp, Patricia J and Novembre, John and Sethuraman, Arun
            and Mathieson, Sara},
        volume = 9,
        year = 2020,
        month = {jun},
        pub_date = {2020-06-23},
        pages = {e54967},
        citation = {eLife 2020;9:e54967},
        doi = {10.7554/eLife.54967},
        url = {https://doi.org/10.7554/eLife.54967},
        keywords = {simulation, reproducibility, open source},
        journal = {eLife},
        issn = {2050-084X},
        publisher = {eLife Sciences Publications, Ltd},
    }


Licence and usage
-----------------

``stdpopsim`` is available under the GPLv3 public license.
The terms of this license can be read
`here <https://www.gnu.org/licenses/gpl-3.0.en.html>`_.
