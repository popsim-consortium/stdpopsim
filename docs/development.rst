.. _sec_development:

===========
Development
===========

If you would like to help with ``stdpopsim`` development, please read the
following. If you think there is anything missing,
please open an `issue <http://github.com/popgensims/stdpopsim/issues>`_ or
`pull request <http://github.com/popgensims/stdpopsim/pulls>`_ on GitHub!

************
Installation
************

Before installing, be sure to make a fork of the repo and clone it locally
following the instructions in the `GitHub Workflow`_.

The ``stdpopsim`` library requires Python 3.4 or later.

For ``pip`` users, install the packages required for development using::

    $ python3 -m pip install -r requirements/development.txt



For ``conda`` users, you will need to add the conda-forge channel to your conda
environment and then should be able to install the development requirements using::

    $ conda config --add channels conda-forge
    $ conda install --file=requirements/development.txt


We do require ``msprime``, so please see the the `installation notes
<https://msprime.readthedocs.io/en/stable/installation.html>`_ if you
encounter problems with it.

.. Note:: If you have trouble installing any of the requirements, your ``pip`` may be the wrong version.
    Try ``pip3 install -r requirements/development.txt``

.. Warning:: The dependency ``daiquiri`` is not currently a conda package.
    So, ``conda install`` will fail.
    See the `GitHub issue
    <https://github.com/popgensims/stdpopsim/issues/161>`_

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
       `rebase <https://help.github.com/articles/about-git-rebase/>`_. Please
       see the next section on how to rebase.

--------
Rebasing
--------

Rebasing is used for two basic tasks we might ask for during review:

1. Your topic branch has gotten out of date with the tip of ``upstream/master``
   and needs to be updated.
2. Your topic branch has lots of messy commits, which need to be cleaned up
   by "squashing".

`Rebasing <https://help.github.com/articles/about-git-rebase/>`_ in git
basically means changing where your branch forked off the main code
in ``upstream/master``. A good way of visualising what's happening is to
look at the `Network <https://github.com/popgensims/stdpopsim/network>`_ view on
GitHub. This shows you all the forks and branches that GitHub knows about
and how they relate to the main repository. Rebasing lets you change where
your branch splits off.

To see this for your local repo
on your computer, you can look at the Git graph output via the command line::

    $  git log --decorate --oneline --graph

This will show something like:

.. code-block:: none

    |*   923ab2e Merge pull request #9 from mcveanlab/docs-initial
    |\
    | * 0190a92 (origin/docs-initial, docs-initial) First pass at development docs.
    | * 2a5fc09 Initial outline for docs.
    | * 1ccb970 Initial addition of docs infrastructure.
    |/
    *   c49601f Merge pull request #8 from mcveanlab/better-genomes
    |\
    | * fab9310 (origin/better-genomes, better-genomes) Added pongo tests.
    | * 62c9560 Tidied up example.
    | * 51e21e8 Added basic tests for population models.
    | * 6fff557 Split genetic_maps into own module.
    | * 90d6367 Added Genome concept.
    | * e2aaf95 Changed debug to info for logging on download.
    | * 2fbdfdc Added badges for CircleCI and CodeCov.
    |/
    *   c66b575 Merge pull request #5 from mcveanlab/tests-ci
    |\
    | * 3ae454f (origin/tests-ci, tests-ci) Initial circle CI config.
    | * c39415a Added basic tests for genetic map downloads.
    |/
    *   dd47000 Merge pull request #3 from mcveanlab/recomb-map-infrastructure
    |\

This shows a nice, linear git history: we can see four pull requests, each of
which consists of a small number of meaningful commits. This is the ideal that
we're aiming for, and git allows us to achieve it by *rewriting history* as
much as we want within our own forks (we never rewrite history in the
``upstream`` repository, as this would cause problems for other developers).
Having a clean, linear git history is a good idea for lots of reasons, not
least of which is making `git bisect <https://git-scm.com/docs/git-bisect>`_
easier.

One of the most useful things that we can do with rebasing is to "squash" commits
so that we remove some noise from the git history. For example, this PR
(on the branch ``topic_branch_name``) currently looks like:

.. code-block:: none

    $  git log --decorate --oneline --graph

    * 97a9458 (HEAD -> topic_branch_name) DONE!!!
    * c9c4a28 PLEASE work, CI!
    * ad4c807 Please work, CI!
    * 0fe6dc4 Please work, CI!
    * 520e6ac Add documentation for rebasing.
    *   20fb835 (upstream/master) Merge pull request #22 from mcveanlab/port-tennyson
    |\
    | * b3d45ea (origin/port-tennyson, port-tennyson) Quickly port Tennesen et al model.
    |/
    *   79d26b4 Merge pull request #20 from andrewkern/fly_model
    |\

Here, in my initial commit (520e6ac) I've added some updated documentation for rebasing.
Then, there's four more commits where I'm trying
to get CI pass. History doesn't need to know about this, so I can rewrite it
using rebase:

.. code-block:: none

    $ git fetch upstream
    $ git rebase -i upstream/master

We first make sure that we're rebasing against the most recent version of the
upstream repo. Then, we ask git to perform an interactive rebase against
the ``upstream/master`` branch. This starts up your editor, showing something
like this::

    pick 520e6ac Add documentation for rebasing.
    pick 0fe6dc4 Please work, CI!
    pick ad4c807 Please work, CI!
    pick c9c4a28 PLEASE work, CI!
    pick 97a9458 DONE!!!

    # Rebase 20fb835..97a9458 onto 20fb835 (5 commands)
    #
    # Commands:
    # p, pick = use commit
    # r, reword = use commit, but edit the commit message
    # e, edit = use commit, but stop for amending
    # s, squash = use commit, but meld into previous commit
    # f, fixup = like "squash", but discard this commit's log message
    # x, exec = run command (the rest of the line) using shell
    # d, drop = remove commit
    #
    # These lines can be re-ordered; they are executed from top to bottom.
    #
    # If you remove a line here THAT COMMIT WILL BE LOST.
    #
    # However, if you remove everything, the rebase will be aborted.
    #
    # Note that empty commits are commented out

We want git to squash the last five commits, so we edit the rebase instructions
to look like:

.. code-block:: none

    pick 520e6ac Add documentation for rebasing.
    s 0fe6dc4 Please work, CI!
    s ad4c807 Please work, CI!
    s c9c4a28 PLEASE work, CI!
    s 97a9458 DONE!!!

    # Rebase 20fb835..97a9458 onto 20fb835 (5 commands)
    #
    # Commands:
    # p, pick = use commit
    # r, reword = use commit, but edit the commit message
    # e, edit = use commit, but stop for amending
    # s, squash = use commit, but meld into previous commit
    # f, fixup = like "squash", but discard this commit's log message
    # x, exec = run command (the rest of the line) using shell
    # d, drop = remove commit
    #
    # These lines can be re-ordered; they are executed from top to bottom.
    #
    # If you remove a line here THAT COMMIT WILL BE LOST.
    #
    # However, if you remove everything, the rebase will be aborted.
    #
    # Note that empty commits are commented out

After performing these edits, we then save and close. Git will try to do
the rebasing, and if successful will open another editor screen that
lets you edit the text of the commit message:

.. code-block:: none

    # This is a combination of 5 commits.
    # This is the 1st commit message:

    Add documentation for rebasing.

    # This is the commit message #2:

    Please work, CI!

    # This is the commit message #3:

    Please work, CI!

    # This is the commit message #4:

    PLEASE work, CI!

    # This is the commit message #5:

    DONE!!!

    # Please enter the commit message for your changes. Lines starting
    # with '#' will be ignored, and an empty message aborts the commit.
    #
    # Date:      Tue Mar 5 17:00:39 2019 +0000
    #
    # interactive rebase in progress; onto 20fb835
    # Last commands done (5 commands done):
    #    squash c9c4a28 PLEASE work, CI!
    #    squash 97a9458 DONE!!!
    # No commands remaining.
    # You are currently rebasing branch 'topic_branch_name' on '20fb835'.
    #
    # Changes to be committed:
    #       modified:   docs/development.rst
    #
    #

We don't care about the commit messages for the squashed commits, so we
delete them and end up with:

.. code-block:: none

    Add documentation for rebasing.

    # Please enter the commit message for your changes. Lines starting
    # with '#' will be ignored, and an empty message aborts the commit.
    #
    # Date:      Tue Mar 5 17:00:39 2019 +0000
    #
    # interactive rebase in progress; onto 20fb835
    # Last commands done (5 commands done):
    #    squash c9c4a28 PLEASE work, CI!
    #    squash 97a9458 DONE!!!
    # No commands remaining.
    # You are currently rebasing branch 'topic_branch_name' on '20fb835'.
    #
    # Changes to be committed:
    #       modified:   docs/development.rst

After saving and closing this editor session, we then get something like:

.. code-block:: none

    [detached HEAD 6b8a2a5] Add documentation for rebasing.
    Date: Tue Mar 5 17:00:39 2019 +0000
    1 file changed, 2 insertions(+), 2 deletions(-)
    Successfully rebased and updated refs/heads/topic_branch_name.

Finally, after a successful rebase, you **must force-push**! If you try to
push without specifying ``-f``, you will get a very confusing and misleading
message:

.. code-block:: none

    $ git push origin topic_branch_name
    To github.com:jeromekelleher/stdpopsim.git
    ! [rejected]        topic_branch_name -> topic_branch_name (non-fast-forward)
    error: failed to push some refs to 'git@github.com:jeromekelleher/stdpopsim.git'
    hint: Updates were rejected because the tip of your current branch is behind
    hint: its remote counterpart. Integrate the remote changes (e.g.
    hint: 'git pull ...') before pushing again.
    hint: See the 'Note about fast-forwards' in 'git push --help' for details.

**DO NOT LISTEN TO GIT IN THIS CASE!** Git is giving you is **terrible advice**
which will mess up your branch. What we need to do is replace the state of
the branch ``topic_branch_name`` on your fork on GitHub (the ``upstream`` remote)
with the state of your local branch, ``topic_branch_name``. We do this
by "force-pushing":

.. code-block:: none

    $ git push -f origin topic_branch_name
    Counting objects: 4, done.
    Delta compression using up to 4 threads.
    Compressing objects: 100% (4/4), done.
    Writing objects: 100% (4/4), 4.33 KiB | 1.44 MiB/s, done.
    Total 4 (delta 2), reused 0 (delta 0)
    remote: Resolving deltas: 100% (2/2), completed with 2 local objects.
    To github.com:jeromekelleher/stdpopsim.git
     + 6b8a2a5...d033ffa topic_branch_name -> topic_branch_name (forced update)

Success! We can check the history again to see if everything looks OK:

.. code-block:: none

    $  git log --decorate --oneline --graph

    * d033ffa (HEAD -> topic_branch_name, origin/topic_branch_name) Add documentation for rebasing.
    *   20fb835 (upstream/master) Merge pull request #22 from mcveanlab/port-tennyson
    |\
    | * b3d45ea (origin/port-tennyson, port-tennyson) Quickly port Tennesen et al model.
    |/
    *   79d26b4 Merge pull request #20 from andrewkern/fly_model
    |

This looks just right: we have one commit, pointing to the head of ``upstream/master``
and have successfully squashed and rebased.

------------------------
When rebasing goes wrong
------------------------

Sometimes rebasing goes wrong, and you end up in a frustrating loop of making and
undoing the same changes over and over again. In this case, it can be simplest to
make a diff of your current changes, and apply these in a single commit. First
we take the diff between the current state of the files in our branch and
``upstream/master`` and save it as a patch::

    $ git diff upstream/master > changes.patch

After that, we can check out a fresh branch and check if everything works
as its supposed to::

    $ git checkout -b test_branch upstream/master
    $ patch -p1 < changes.patch
    $ git commit -a
    # check things work

After we've verified that everything works, we then checkout the original
topic branch and replace it with the state of the ``test_branch``, and
finally force-push to the remote topic branch on your fork::

    $ git checkout topic_branch_name
    $ git reset --hard test_branch
    $ git push -f origin topic_branch_name

Hard resetting and force pushing are not reversible operations, so please
beware!

********************
Model review process
********************

When Developer A creates a new demographic model on their local fork they must
follow these steps for it to be officially supported by stdpopsim:

    1. Developer A submits PR to add a new model. This must include basic unit
       tests and full documentation (i.e., the documentation that will be
       rendered on rtd). The code is checked for any obvious problems/style
       issues etc by maintainer M and merged when it meets these basic
       standards. The model is considered 'preliminary' and not linked from the
       online documentation and not included in the CLI.

    2. Developer A creates an issue tracking the QC for the model which includes
       information about the primary sources used to create the model and the
       population indices used for their msprime implementation. Developer B is
       then assigned/volunteers to do a blind implementation of the model. 

    3. M creates an issue for the CLI implementation of the model.

    4. Developer B creates a blind implementation of the model in the
       ``qc/species_name_qc.py`` file. Note that if you are adding a new species
       you will have to add a new import to ``qc/__init__.py``.
    
    5. Developer B adds the automatic checking of this model for
       equality with the production model to the suite of unit tests in for the
       demograpic model in ``tests/test_species_name_.py`` following the
       template below::

        def test_qc_model_equal(self):
            model = homo_sapiens.BrowningAmerica()
            self.assertTrue(model.equals(homo_sapiens_qc.BrowningAmerica()))

       Developer B then creates a PR, and all being good, this PR is merged and
       the QC issue is closed. 

    6. Someone then makes a PR updating the CLI, checking that the
       documentation, citations etc all work properly, and adds the model to
       the list of documented models.

------------------------
Arbitration
------------------------

When developers A and B disagree on the model implementation, the process is to:

    1. Try to hash out the details between them on the original issue thread

    2. If this fails, contact the authors of the original publication to resolve
       ambiguities.

    3. If changes have to be made to the production model Developer A submits a
       PR with the hotfix for the production model. Developer B then rebases
       the branch containing their PR against master to check for model
       equality. Repeat steps 1-3 until this is achieved. If changes have to be
       made to the QC model they are committed to the branch where the QC PR
       originates from.

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

It's useful to run the ``flake8`` CI tests *locally* before pushing a commit.
To set this up use either ``pip`` or ``conda`` to install ``flake8``

To run the test simply use::

    $ flake8 --max-line-length 89 stdpopsim tests

If you would like to automatically run this test before a commit is permitted,
add the following line in the file ``stdpopsim/.git/hooks/pre-commit.sample``::

    exec flake8 --max-line-length 89 setup.py stdpopsim tests 
 
before::

    # If there are whitespace errors, print the offending file names and fail.
    exec git diff-index --check --cached $against --

Finally, rename ``pre-commit.sample`` to simply ``pre-commit``

*************
Documentation
*************

Documentation is written using `reStructuredText <http://docutils.sourceforge.net/rst.html>`_
markup and the `sphinx <http://www.sphinx-doc.org/en/master/>`_ documentation system.
It is defined in the ``docs`` directory.

To build the documentation type ``make`` in the ``docs`` directory. This should build
HTML output in the ``_build/html/`` directory.


