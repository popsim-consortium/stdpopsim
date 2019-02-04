.. _sec_development:

===========
Development
===========

If you would like to help with ``stdpopsim`` development, please read the
following. If you think there is anything missing,
please open an `issue <http://github.com/popgensims/stdpopsim/issues>`_ or
`pull request <http://github.com/popgensims/stdpopsim/pulls>`_ on GitHub!

********
Overview
********

The ``stdpopsim`` library requires Python 3.4 or later. For ``pip`` users,
install the packages required for development using::

    $ python3 -m pip install requirements/development.txt

We do require ``msprime``, so please see the the `installation notes
<https://msprime.readthedocs.io/en/stable/installation.html>` if you
encounter problems with it. Conda users should be able to install the
development requirements using::

    $ conda install --file=requirements/development.txt

***************
GitHub workflow
***************

    1. Make your own `fork <https://help.github.com/articles/fork-a-repo/>`_
       of the ``stdpopsim`` repository on GitHub, and
       `clone <https://help.github.com/articles/cloning-a-repository/>`_
       a local copy.
    2. Make sure that your local repository has been configured with an
       `upstream remote <https://help.github.com/articles/configuring-a-remote-for-a-fork/>`_.
    3. Create a "topic branch" to work on. One reliable way to do it
       is to follow this recipe::

        $ git fetch upstream
        $ git checkout upstream/master
        $ git checkout -b topic_branch_name

    4. As you work on your topic branch you can add commits to it. Once you're
       ready to share this, you can then open a `pull request
       <https://help.github.com/articles/about-pull-requests/>`__. Your PR will
       be reviewed by some of the maintainers, who may ask you to make changes.
    5. If your topic branch has been around for a long time and has gotten
       out of date with the main repository, we might ask you to
       `rebase <https://help.github.com/articles/about-git-rebase/>`_. Usually,
       you can do this by::

        $ git fetch upstream
        $ git rebase -i upstream/master

       Rebasing is very powerful (and complicated!), so please see
       `this guide <https://help.github.com/articles/about-git-rebase/>`_
       for more information. Once you have rebased your topic branch, you
       then need to "force push" ::

        $ git push -f origin topic_branch_name


********************
Model review process
********************

.. todo:: Document the review process for adding models.

****************
Coding standards
****************

To ensure that the code in ``stdpopsim`` is as readable as possible
and follows a reasonably uniform style, we require that all code follows
the `PEP8 <https://www.python.org/dev/peps/pep-0008/>`_ style guide.
Lines of code should be no more than 89 characters.
Conformance to this style is checked as part of the Continuous Integration
testing suite.

**********
Unit tests
**********

All code added to ``stdpopsim`` should have
`unit tests <https://en.wikipedia.org/wiki/Unit_testing>`_. These are typically
simple and fast checks to ensure that the code makes basic sense (the
entire unit test suite should not require more than a few seconds to run).
Test coverage is checked using `CodeCov <https://codecov.io/gh/popgensims/stdpopsim>`_,
which generates reports about each pull request.

It is not practical to test the statistical properties of simulation models
as part of unit tests.

The unit test suite is in the ``tests`` directory. Tests are run using the
`nose <https://nose.readthedocs.io/en/latest/>`_ module. Use::

    $ python3 -m nose tests/

from the project root to run the full test suite.

*************
Documentation
*************

Documentation is written using `reStructuredText <http://docutils.sourceforge.net/rst.html>`_
markup and the `sphinx <http://www.sphinx-doc.org/en/master/>`_ documentation system.
It is defined in the ``docs`` directory.


