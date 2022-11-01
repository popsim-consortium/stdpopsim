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
`msprime <https://tskit.dev/software/msprime.html>`_ and
`SLiM 4 <https://messerlab.org/slim/>`_ to generate sample datasets in the
`tree sequence <https://tskit.dev/learn/>`_ format.


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
A primary goal of PopSim Consortium is to be inclusive to the largest number of contributors,
with the most varied and diverse backgrounds possible. As such, we are committed to providing a
friendly, safe and welcoming environment for all, regardless of gender, sexual orientation, ability,
ethnicity, socioeconomic status, and religion (or lack thereof).
Being a member of the Consortium can provide a valuable avenue for career development, through
building networks, contributing to papers and code, as well as through learning programing standards
and best practices.

For further details read the :ref:`Development <sec_development>` page to find out how you can join us.

Citations
---------

If you use ``stdpopsim`` in your work, please cite our
`original manuscript <https://doi.org/10.7554/eLife.54967>`__ and/or the
`followup manuscript <https://doi.org/10.1101/2022.10.29.514266>`__ describing
major additions to the catalog:

  - Jeffrey R Adrion, Christopher B Cole, Noah Dukler, Jared G Galloway,
    Ariella L Gladstein, Graham Gower, Christopher C Kyriazis, Aaron P Ragsdale,
    Georgia Tsambos, Franz Baumdicker, Jedidiah Carlson, Reed A Cartwright,
    Arun Durvasula, Ilan Gronau, Bernard Y Kim, Patrick McKenzie,
    Philipp W Messer, Ekaterina Noskova, Diego Ortega-Del Vecchyo, Fernando Racimo,
    Travis J Struck, Simon Gravel, Ryan N Gutenkunst, Kirk E Lohmeuller,
    Peter L Ralph, Daniel R Schrider, Adam Siepel, Jerome Kelleher, Andrew D Kern (2020),
    *A community-maintained standard library of population genetic models*,
    eLife 2020;9:e54967; doi: https://doi.org/10.7554/eLife.54967

  - M. Elise Lauterbur, Mariz Izabel A. Cavassim, Ariella L. Gladstein, Graham
    Gower, Nathaniel S. Pope, Georgia Tsambos, Jeffrey Adrion, Saurabh Belsare,
    Arjun Biddanda, Victoria Caudill, Jean Cury, Ignacio Echevarria, Benjamin C.
    Haller, Ahmed R. Hasan, Xin Huang, Leonardo Nicola Martin Iasi, Ekaterina
    Noskova, Jana Obšteter, Vitor Antonio Corrêa Pavinato, Alice Pearson, David
    Peede, Manolo F. Perez, Murillo F. Rodrigues, Chris C. R. Smith, Jeffrey P.
    Spence, Anastasia Teterina, Silas Tittes, Per Unneberg, Juan Manuel Vazquez,
    Ryan K. Waples, Anthony Wilder Wohns, Yan Wong, Franz Baumdicker, Reed A.
    Cartwright, Gregor Gorjanc, Ryan N. Gutenkunst, Jerome Kelleher, Andrew D.
    Kern, Aaron P. Ragsdale, Peter L. Ralph, Daniel R. Schrider, Ilan Gronau (2022).
    *Expanding the stdpopsim species catalog, and lessons learned for realistic genome simulations*,
    bioRxiv 2022.10.29.514266; doi: https://doi.org/10.1101/2022.10.29.514266

Bibtex records::

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

    @article {Lauterbur2022.10.29.514266,
        author = {Lauterbur, M. Elise and Cavassim, Mariz Izabel A. and
            Gladstein, Ariella L. and Gower, Graham and Pope, Nathaniel S. and
            Tsambos, Georgia and Adrion, Jeffrey and Belsare, Saurabh and Biddanda,
            Arjun and Caudill, Victoria and Cury, Jean and Echevarria, Ignacio and
            Haller, Benjamin C. and Hasan, Ahmed R. and Huang, Xin and Iasi,
            Leonardo Nicola Martin and Noskova, Ekaterina and Ob{\v s}teter, Jana
            and Pavinato, Vitor Antonio Corr{\^e}a and Pearson, Alice and Peede,
            David and Perez, Manolo F. and Rodrigues, Murillo F. and Smith, Chris
            C. R. and Spence, Jeffrey P. and Teterina, Anastasia and Tittes, Silas
            and Unneberg, Per and Vazquez, Juan Manuel and Waples, Ryan K. and
            Wohns, Anthony Wilder and Wong, Yan and Baumdicker, Franz and
            Cartwright, Reed A. and Gorjanc, Gregor and Gutenkunst, Ryan N. and
            Kelleher, Jerome and Kern, Andrew D. and Ragsdale, Aaron P. and Ralph,
            Peter L. and Schrider, Daniel R. and Gronau, Ilan},
        title = {Expanding the stdpopsim species catalog, and lessons learned
            for realistic genome simulations},
        elocation-id = {2022.10.29.514266},
        year = {2022},
        doi = {10.1101/2022.10.29.514266},
        publisher = {Cold Spring Harbor Laboratory},
        URL = {https://www.biorxiv.org/content/early/2022/10/31/2022.10.29.514266},
        eprint = {https://www.biorxiv.org/content/early/2022/10/31/2022.10.29.514266.full.pdf},
        journal = {bioRxiv}
    }


Licence and usage
-----------------

``stdpopsim`` is available under the GPLv3 public license.
The terms of this license can be read
`here <https://www.gnu.org/licenses/gpl-3.0.en.html>`_.
