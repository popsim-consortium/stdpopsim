.. _sec_development:

===========
Development
===========

``stdpopsim`` **is a community effort, and we welcome YOU to join us!**

We envision at least three main types of ``stdpopsim`` developers:

1. Contributors of new species, demographic models, and other features
   (such as recombination maps, annotations, DFE models)
2. API developers
3. Documentation and tutorial curators

Contributors add simulation code for new or existing species in the catalog.
This is the main way we envision biologists to continually add
to the catalog of available models, and it is a great first step for new
contributors to learn the ins and outs of ``stdpopsim`` development.
See the appropriate sections below:

* `Adding a new species`_
* `Adding a new demographic model`_
* `Adding a genetic map`_
* `Adding a DFE model`_

`API developers` work on infrastructure development for the PopSim Consortium,
which could include improvements and additions to the internal code base of
``stdpopsim``, establishment of benchmarking pipelines,
and new projects that align with consortium goals.

`Documentation and tutorial curators` help maintain the documentation and tutorials.
This can be as easy as pointing out confusing bits of the documentation in a
GitHub `issue <http://github.com/popgensims/stdpopsim/issues>`__., or adding or editing
the documentation. See the section `Documentation`_.

Get into contact with the ``stdpopsim`` community by subscribing to our
`email list serve <https://lists.uoregon.edu/mailman/listinfo/popgen_benchmark>`__
and by creating and commenting on
Github `issue <http://github.com/popgensims/stdpopsim/issues>`__.
There is a lot of chatter through
`Github <http://github.com/popgensims/stdpopsim>`__, and we’ve been building code
there cooperatively.
If you want to help out and don't know where to start, you can look through the
list of
`Good first issues
<https://github.com/popgensims/stdpopsim/issues?q=is%3Aopen+is%3Aissue+label%3A%22
good+first+issue%22>`__
or
`Help wanted issues
<https://github.com/popgensims/stdpopsim/issues?q=is%3Aopen+is%3Aissue+label%3A%22
help+wanted%22>`__


To get started helping with ``stdpopsim`` development, please read the
following sections to learn how to contribute.
And, importantly, please have a look at our
`code of conduct <https://github.com/popsim-consortium/stdpopsim/blob/main/CODE_OF_CONDUCT.md>`__.

.. _sec_development_installation:

************
Installation
************

Before installing, be sure to make a fork of the repo and clone it locally
following the instructions in the `GitHub Workflow`_.

The ``stdpopsim`` library requires Python 3.4 or later.

For ``pip`` users, install the packages required for development using::

    $ python3 -m pip install -r requirements/development.txt

You can then install the development version of ``stdpopsim`` like this::

    $ python3 setup.py install

For ``conda`` users, you will need to add the conda-forge channel to your conda
environment and then should be able to install the development requirements using::

    $ conda config --add channels conda-forge
    $ conda install --file=requirements/development.txt


We do require ``msprime``, so please see the the `installation notes
<https://tskit.dev/msprime/docs/stable/installation.html>`__ if you
encounter problems with it.

.. Note:: If you have trouble installing any of the requirements, your ``pip`` may be the wrong version.
    Try ``pip3 install -r requirements/development.txt``

---------------------------
Using a Virtual Environment
---------------------------

We encourage the use of a virtual environment.

For ``pip``, you can use ``venv``.

First, create the virtual environment (You only need to do this once)::

    $ python3 -m venv stdpopsim_env

Next, activate the virtual environment::

    $ source stdpopsim_env/bin/activate

You will then see the virtual environment in your prompt. Like so::

    (stdpopsim_env) $

Once the virtual environment is activated, install the requirements::

    (stdpopsim_env) $ python3 -m pip install -r requirements/development.txt

You can then run any of the code in the virtual environment with the packages installed,
without conflicting with other packages in your local environment.
To deactivate the virtual environment::

    (stdpopsim_env) $ deactivate


***************
GitHub workflow
***************

    1. Make your own `fork <https://help.github.com/articles/fork-a-repo/>`__
       of the ``stdpopsim`` repository on GitHub, and
       `clone <https://help.github.com/articles/cloning-a-repository/>`__
       a local copy.
    2. Install the pre-commit hooks with::

        $ pre-commit install

    2. Make sure that your local repository has been configured with an
       `upstream remote <https://help.github.com/articles/configuring-a-remote-for-a-fork/>`__.
    3. Create a "topic branch" to work on. One reliable way to do it
       is to follow this recipe::

        $ git fetch upstream
        $ git checkout upstream/main
        $ git checkout -b topic_branch_name

    4. As you work on your topic branch you can add commits to it. Once you're
       ready to share this, you can then open a `pull request
       <https://help.github.com/articles/about-pull-requests/>`__. Your PR will
       be reviewed by some of the maintainers, who may ask you to make changes.
    5. If your topic branch has been around for a long time and has gotten
       out of date with the main repository, we might ask you to
       `rebase <https://help.github.com/articles/about-git-rebase/>`__. Please
       see the next section on how to rebase.

-----------------
Pre-commit checks
-----------------

On each commit a `pre-commit hook <https://pre-commit.com/>`__  will run
that checks for violations of code style and other common problems.
Where possible, these hooks will try to fix any problems that they find (including reformatting
your code to conform to the required style). In this case, the commit
will *not complete* and report that "files were modified by this hook".
To include the changes that the hooks made, ``git add`` any
files that were modified and run ``git commit`` (or, use ``git commit -a``
to commit all changed files.)

If you would like to run the checks without committing, use ``pre-commit run``
(but, note that this will *only* check changes that have been *staged*;
do ``pre-commit run --all`` to check unstaged changes as well).
To bypass the checks (to save or get feedback on work-in-progress) use
``git commit --no-verify``

--------
Rebasing
--------

Rebasing is used for two basic tasks we might ask for during review:

1. Your topic branch has gotten out of date with the tip of ``upstream/main``
   and needs to be updated.
2. Your topic branch has lots of messy commits, which need to be cleaned up
   by "squashing".

`Rebasing <https://help.github.com/articles/about-git-rebase/>`__ in git
basically means changing where your branch forked off the main code
in ``upstream/main``. A good way of visualising what's happening is to
look at the `Network <https://github.com/popgensims/stdpopsim/network>`__ view on
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
least of which is making `git bisect <https://git-scm.com/docs/git-bisect>`__
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
    *   20fb835 (upstream/main) Merge pull request #22 from mcveanlab/port-tennyson
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
    $ git rebase -i upstream/main

We first make sure that we're rebasing against the most recent version of the
upstream repo. Then, we ask git to perform an interactive rebase against
the ``upstream/main`` branch. This starts up your editor, showing something
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

**DO NOT LISTEN TO GIT IN THIS CASE!** Git is giving you **terrible advice**
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
    *   20fb835 (upstream/main) Merge pull request #22 from mcveanlab/port-tennyson
    |\
    | * b3d45ea (origin/port-tennyson, port-tennyson) Quickly port Tennesen et al model.
    |/
    *   79d26b4 Merge pull request #20 from andrewkern/fly_model
    |

This looks just right: we have one commit, pointing to the head of ``upstream/main``
and have successfully squashed and rebased.

------------------------
When rebasing goes wrong
------------------------

Sometimes rebasing goes wrong, and you end up in a frustrating loop of making
and undoing the same changes over and over again. First, here's an explanation
of what's going on. Let's say that the branch we're working on (and trying to
rebase) is called ``topic_branch``, and it branched off from ``upstream/main``
at some point in the past::

         A1---A2---A3  (topic_branch)
        /
    ---M---o---o---o---o---B  (upstream/main)

So, what we'd really like to do is to take the commits ``A1``, ``A2``, and
``A3`` and apply them to the current state of the ``upstream/main`` branch,
i.e., on top of commit ``B``. If we just do ``git rebase upstream/main``
then git will try to first apply ``A1``; then ``A2``; and finally ``A3``.
If there's conflicts, this is painful, so we might want to *first* squash
the three commits together into one commit, and then rebase that single commit.
Then we'll only have to resolve conflicts once. Said another way: we often
use ``git rebase -i upstream/main`` to both squash *and* rebase; but
it may be easier to squash first then rebase after.

We'll be doing irreversible changes, so first we should make a backup copy of
the branch::

    $ git checkout topic_branch  # make sure we're on the right branch
    $ git checkout -b topic_backup # make the backup
    $ git checkout topic_branch  # go back to the topic branch

Next, we take the diff between the current state of the files and the place
where your changes last diverged from ``upstream/main`` (the commit labelled
``M`` in the diagram above), and save it as a patch. To do this, make sure
you are in the root of the git directory, and::

    $ git diff --merge-base upstream/main > changes.patch

After that, we can check out a fresh branch and check if everything works
as it's supposed to::

    $ git checkout -b test_branch upstream/main
    $ patch -p1 < changes.patch
    $ git commit -a
    # check things work

After we've verified that everything works, we then checkout the original
topic branch and replace it with the state of the ``test_branch``, and
finally force-push to the remote topic branch on your fork::

    $ git checkout topic_branch
    $ git reset --hard test_branch
    $ git push -f origin topic_branch

Hard resetting and force pushing are not reversible operations, so please
beware! After you've done this, you can go make sure nothing bad happened
by checking that the only changes listed under "files changed" in the github
pull request are changes that you have made. For more on finding the fork
point, with diagrams, and an alternative workflow, see `the git docs
<https://git-scm.com/docs/git-merge-base>`__.

.. _sec_development_demographic_model:


********************
Adding a new species
********************

---------------------------------------------------
Which information do I need to have for my species?
---------------------------------------------------

In ``stdpopsim``, we aim to be inclusive and facilitate adding a diverse range of species.
That said, there are certain basic requirements we have
for every species added to the :ref:`sec_catalog`.
We specify these requirements below.
If you are unsure whether your species satisfies these baseline requirements,
but you still think it will be useful to add it to ``stdpopsim``,
then we encourage you to `open an issue <http://github.com/popgensims/stdpopsim/issues/new>`__
on the GitHub repository to discuss this.
Others researchers in the community may be able to help you fill in the missing details
or find other solutions.

Every species added to ``stdpopsim`` should have the following information available:

1. A chromosome-level genome assembly
2. Mutation rate (per generation)
3. Recombination rate (per generation)
4. A characteristic population size
5. An average generation time

Of course, many species do not have precise estimates of each of these
(e.g., mutation rates are usually not known).
So, in practice we often have to use approximate estimates.
We provide below a set of guidelines for each of the five components,
with a brief discussion of possible courses of action to take when components have incomplete information.

1. The **genome assembly** should consist of a list of chromosomes or scaffolds and their lengths.
   Having a good quality assembly with complete chromosomes, or at least very long scaffolds,
   is essential for chromosome-level simulations produced by ``stdpopsim``.
   Species with less complete genome builds typically do not have genetic maps
   or good estimates of recombination rates,
   making chromosome-level simulation much less useful.
   Thus, currently, ``stdpopsim`` only supports adding species with near-complete
   chromosome-level genome assemblies (i.e., close to one contig per chromosome).

2. An **average mutation rate**
   should be specified for each chromosome (per generation per bp).
   The mutation rate estimate can be based on sequence data from pedigrees, mutation accumulation studies,
   or comparative genomic analysis calibrated by fossil data (i.e., phylogenetic estimates).
   If there is no information on the variation of mutation rates across chromosomes,
   the average genome-wide mutation rate can be specified for all chromosomes.
   Finally, if your species of interest does not have direct estimates of mutation rates,
   we recommend using estimates for some other species (hopefully closely related).

3. An **average recombination rate**
   should be specified for each chromosome (per generation per bp).
   Ideally, one would want to specify a fine-scale chromosome-level **recombination map**,
   since the recombination rate is known to vary widely across chromosomes.
   If a recombination map exists for your species,
   you may specify it separately (see `Adding a genetic map`_).
   Nonetheless, you should specify a default (average) recombination rate for each chromosome.
   As with mutation rates, if there is no information on the variation of recombination rates
   across chromosomes, the average genome-wide recombination rate can be specified for all chromosomes.
   Furthermore, if your species of interest does not have direct estimates of recombination rates,
   we recommend using estimates for some other species (hopefully closely related).

4. The **effective population size** should represent the historical average effective population size,
   and should ideally produce simulated data that matches the average observed genetic diversity in that species.
   Population size is defined as the number of individuals, regardless of ploidy.
   However, this will often not capture features of genetic variation that are caused by recent changes in population size and the presence of population structure.
   To capture those, one should also provide a demographic model (or multiple models) for the species
   (see `Adding a new demographic model`_).

5. The **average generation time** is an important part of the species' natural history,
   but its value does not directly affect the simulation, since
   the ``SLiM`` and ``msprime`` simulation engines operate in time units of generations.
   Thus, the average generation time is only currently used to convert time units to years,
   which is useful when comparing different demographic models.

All values used in the species model should be based on current knowledge for a typical population
in that species, as represented in the literature.
Before you add your species to ``stdpopsim``, see that you can collect the values
mentioned above from the literature.
You will later need to specify these citations in your code files
(see `Coding the species parameters`_).
If you are unsure whether your species of interest satisfies the base requirements above
(such as a near-complete genome assembly), or have questions about how to set some parameters,
feel free to `open an issue <http://github.com/popgensims/stdpopsim/issues/new>`__
on the GitHub repository to get assistance from other members of the ``stdpopsim`` community.


-----------------------------------
Getting set up to add a new species
-----------------------------------

If this is your first time adding a species to ``stdpopsim``, it's a good
idea to take some time browsing the :ref:`sec_catalog`
to see how existing species are typically specified and documented. If you have
any questions or confusion about the required code, please
don't hesitate to
`open a new issue <https://github.com/popsim-consortium/stdpopsim/issues/new>`__.
We're more than happy to answer any questions and help get you up and running.
Before you add any code, be sure to have forked the ``stdpopsim`` repository
and cloned it locally, following the instructions in the `GitHub Workflow`_ section.


After you collected the relevant parameters from the literature (see list above),
the first step is to create a new subdirectory devoted to the new species,
which you should name using the six-character species identifier
(see `Naming conventions`_ for more details).
All code associated with simulation of this species should go into this directory,
unless explicitly specified otherwise
(code for documentation and testing  is written in other directories).
For example, the simulation code for *D. melanogaster* resides in directory
``stdpopsim/catalog/DroMel/`` in the repository.

Once the species directory is set up, you may use the ``maintenance`` utility
of ``stdpopsim`` to generate template files where you can enter
all relevant information for your species.
The ``maintenance`` utility downloads useful information on a genome build published
in `Ensembl <https://www.ensembl.org/index.html>`__,
and uses it to generate initial versions of the required source files.
A partial list of the
genomes housed on Ensembl can be found `here <https://metazoa.ensembl.org/species.html>`__.
To use this utility, execute the ``maintenance`` command with the Ensembl species ID;
replace spaces in the Ensembl ID with ``_`` characters.
For example, the template files for *A. thaliana* were generated by executing this command:

.. code-block:: shell

    $ python -m maintenance add-species arabidopsis_thaliana

The ``maintenance`` utility generates three new files inside the species directory
(``stdpopsim/catalog/<SPECIES_ID>/``):

* ``__init__.py``: a  script that loads all the relevant libraries for your species.
  It should be edited only when you add components to your species, such as demographic models,
  genetic maps, or DFE models.

* ``genome_data.py``: a file that contains information on the physical map of the genome.
  This file is generated automatically by the ``maintenance`` utility with a data dictionary
  which has slots for the assembly accession number, the assembly name,
  and a dictionary representing the chromosome names and their associated lengths.
  If synonyms are defined (e.g., chr2L for 2L) then those are given in the list that follows.
  You should double-check the downloaded values, but there is probably no reason to edit this file
  after it has been generated by the ``maintenance`` utility.

* ``species.py``: a file containing information about the species' mutation and recombination rates,
  effective population size, and the average generation time,
  along with all accompanying citations
  (see details in `Which information do I need to have for my species?`_).
  The following section provides detailed instructions on how to code information in this file,
  including some specific examples.

.. note::

      The ``maintenance`` utility also generates test code for your species in
      the file ``tests/test_<SPECIES_ID>.py``.
      This is used later for your local tests and in the review process
      (see `Testing your species model and submitting a PR`_
      and `Implementing tests for the review of new species`_).

.. note::

    If your species of interested does not have a published genome in Ensembl,
    you may manually create and edit the three files described above.
    Try to follow an example from the catalog that was downloaded from Ensembl
    to maintain a consistent format.

-----------------------------
Coding the species parameters
-----------------------------

Information about a species' mutation and recombination rates,
effective population size, and the average generation time,
is all summarized in the ``species.py`` file,
along with all accompanying citations
(see details in `Which information do I need to have for my species?`_).
The initial version of the file generated by the ``maintenance`` utility
contains commented instructions to help you figure out where everything goes.
Essentially, the information in this file is recorded in two main objects: ``_genome`` and ``_species``.
The ``_genome`` object contains chromosome-level information, such as
**chromosome ids**, **lengths**, **mutation and recombination rates**, and **ploidy**.
The ``_species`` object contains the remaining information about the species,
including its **full name**, **abbreviated name**, **id**, **effective population size**
and **average generation time**.
Each value specified in these two object should be accompanied by a
``stdpopsim.Citation`` object indicating the publication from which it was derived.
Each ``stdpopsim.Citation`` object is initialized with the following information:

* author (`string`): abbreviated author list in a single string,
  such as `"1000GenomesConsortium"` or `"Huber et al."`.
* year   (`int`): year of publication.
* doi (`string`): a URL for the `doi.org <https://doi.org/>`__ webpage of the publication.
* reasons (list of ``stdpopsim.CiteReason``):
  possible reasons to include a citation in ``species.py`` are:

  * ``stdpopsim.CiteReason.ASSEMBLY``
  * ``stdpopsim.CiteReason.REC_RATE``
  * ``stdpopsim.CiteReason.MUT_RATE``
  * ``stdpopsim.CiteReason.POP_SIZE``
  * ``stdpopsim.CiteReason.GEN_TIME``

To demonstrate how the ``_genome`` and ``_species`` objects are set,
we provide below a detailed example for *A. thaliana*
(see also ``stdpopsim/catalog/AraTha/species.py``).

We start by defining auxiliary objects that specify the recombination rate,
mutation rate, and ploidy for each chromosome.
In the case of *A. Thaliana*, these objects are defined to associate
the mitochondrial and plastid genomes (chromsoomes `Mt` and `Pt`)
with ploidy of 1 and recombination rate of 0.
All other chromosomes are associated with a ploidy of 2 and the
genome-wide average recombination rate.
The genome-wide mutation rate is associated with all chromosomes.

.. code-block:: python

  # genome-wide recombination rate from Huber et al 2014 MBE
  # associated with all recombining chromosomes
  _rho = 200 / 1e6  # 200/Mb
  _Ne = 124000
  _mean_recombination_rate = _rho / (2 * _Ne)
  _recombination_rate = {str(j): _mean_recombination_rate for j in range(1, 6)}
  _recombination_rate["Mt"] = 0
  _recombination_rate["Pt"] = 0

  # genome-wide average mutation rate from Ossowski 2010 Science
  # associated with all chromosomes
  _mean_mutation_rate = 7e-9
  _mutation_rate = {str(j): _mean_mutation_rate for j in range(1, 6)}
  _mutation_rate["Mt"] = _mean_mutation_rate
  _mutation_rate["Pt"] = _mean_mutation_rate

  # species ploidy and chromosome-specific ploidy
  _species_ploidy = 2
  _ploidy = {str(j): _species_ploidy for j in range(1, 6)}
  _ploidy["Mt"] = 1
  _ploidy["Pt"] = 1


The ``_genome`` object is then defined by calling the ``stdpopsim`` function
``stdpopsim.Genome.from_data``.
This functions generates the genome object based on information from the
``data`` object defined in the ``genome_data.py`` file,
the auxiliary objects defined above,
and a list of ``stdpopsim.Citation`` objects.

.. code-block:: python

  _genome = stdpopsim.Genome.from_data(
      genome_data.data,
      recombination_rate=_recombination_rate,
      mutation_rate=_mutation_rate,
      ploidy=_ploidy,
      citations=[
          stdpopsim.Citation(
              author="Ossowski et al.",
              year=2010,
              doi="https://doi.org/10.1126/science.1180677",
              reasons={stdpopsim.CiteReason.MUT_RATE},
          ),
          stdpopsim.Citation(
              author="Huber et al.",
              year=2014,
              doi="https://doi.org/10.1093/molbev/msu247",
              reasons={stdpopsim.CiteReason.REC_RATE},
          ),
          stdpopsim.Citation(
              doi="https://doi.org/10.1093/nar/gkm965",
              year=2007,
              author="Swarbreck et al.",
              reasons={stdpopsim.CiteReason.ASSEMBLY},
          ),
      ],
  )



The ``_species`` object contains a reference to the ``_genome`` object and
the remaining information about the species,
including the **effective population size** and **average generation time**,
accompanied by the appropriate ``stdpopsim.Citation`` objects.

.. code-block:: python

    _species = stdpopsim.Species(
        id="AraTha",
        ensembl_id="arabidopsis_thaliana",
        name="Arabidopsis thaliana",
        common_name="A. thaliana",
        genome=_genome,
        generation_time=1.0,
        population_size=10 ** 4,
        ploidy=_species_ploidy,
        citations=[
            stdpopsim.Citation(
                doi="https://doi.org/10.1890/0012-9658(2002)083[1006:GTINSO]2.0.CO;2",
                year=2002,
                author="Donohue",
                reasons={stdpopsim.CiteReason.GEN_TIME},
            ),
            stdpopsim.Citation(
                doi="https://doi.org/10.1016/j.cell.2016.05.063",
                year=2016,
                author="1001GenomesConsortium",
                reasons={stdpopsim.CiteReason.POP_SIZE},
            ),
        ],
    )


Once these two objects (``_genome`` and ``_species``) are specified in the ``species.py`` file,
you should be able to load and simulate the newly added species using ``stdpopsim``.

----------------------------------------------
Testing your species model and submitting a PR
----------------------------------------------

The ``maintenance`` utility that generated the three species template files
in the species directory (``stdpopsim/catalog/<SPECIES_ID>/``)
also generates test code for the species in a separate file, ``tests/test_<SPECIES_ID>.py``.
The tests in this file are executed as follows
(where ``<SPECIES_ID>`` is the six-character species id):

.. code-block:: shell

   $ python -m pytest tests/test_<SPECIES_ID>.py

The tests already implemented in this file when it is generated
check for basic formatting and missing information.
For example, there is a test checking that the citation year is of type `int`
rather than `string` (e.g. 2004 and not `"2014"`).
Other tests in this file are generated by the ``maintenance`` utility
as blank and disabled.
These tests should **not** be filled out by the person who writes the code in
the ``species.py`` file,
but rather by someone else, as part of the **review process** (see below).
Once your code passes the basic tests implemented in the automatically generated
version of the test file,
you should submit a pull request (PR) with your changes to the catalog.
See the `GitHub workflow`_ for more details about this process.

At this point, most of your work is done.
**You have officially joined the** ``stdpopsim`` **development team. Welcome!!**
Your code still needs to undergo review by another member (or members)
of the development team before it is fully incorporated into ``stdpopsim``.
This will likely require additional feedback from you,
so, stay tuned for discussion during the review process.

----------------------------------------
Overview of the stdpopsim review process
----------------------------------------

We provide here a general outline for the review process we use in ``stdpopsim``,
including guidelines for how to settle discrepancies that are found during review
(see Step 6 below).
The seven steps described below should be followed whenever a **new species** is added,
or when components such as **demographic models** are added to a species
already in the catalog.

1. After the original contributor submitted a PR with their new code,
   the code is checked by one of the core maintainers of
   ``stdpopsim`` for basic problems or style issues.
   Once the code meets the basic standards, the maintainer merges the PR,
   and the newly added code is considered **provisional**.

2. The original contributor then opens a new **QC issue** on GitHub
   to track the progress of the review.
   One simple way to do this is to use one of the `template issues
   <https://github.com/popsim-consortium/stdpopsim/issues/new/choose>`__
   we provide.
   For example, the ``Species QC issue template`` should be used when adding
   a new species and the ``Model QC issue template`` should be used when adding
   a new demographic model.
   Simply press  `Get started` for the appropriate template,
   and fill in the required details.
   If you don't find an appropriate template for your purpose,
   you should simply `open a new blank issue
   <https://github.com/popsim-consortium/stdpopsim/issues/new>`__
   and add the relevant details manually.
   Make sure to include information about the primary sources (citations)
   you used as well as other considerations you made in your code.
   The **QC issue** contains a checklist and all the items on this list
   should be checked off for the review process to complete.

3. A different member of the ``stdpopsim`` community volunteers to review the
   newly added demographic model.
   If you volunteer to review a model, you should state your intention on the
   **QC issue**, so we don't duplicate effort.
   Typically, there will be one reviewer assigned to every **QC issue**.
   However, sometimes multiple reviewers may wish to partition tasks between them.
   For examples, when reviewing a new species, one reviewer may wish to test the
   recombination rates, and another may wish to test the effective population size.
   Some aspects of the review, such as examining citations, involve checking the
   code of the original contributor.
   However, most of the review involves implementing tests
   based on the reviewer's understanding of
   the source publications and additional documentation
   specified by the original contributor in the **QC issue**.
   Ideally, the code for these tests should be written by the reviewer
   **without looking at the original contributor's code**.
   If the reviewer is uncertain about some aspects of the implementation,
   they can discuss this with the original contributor in the **QC issue**.
   Different types of tests are involved when you are reviewing a **new species**
   added to ``stdpopsim`` or when you are reviewing a **demographic model**
   added to an existing species.
   See the appropriate sections below for specific instructions on how to
   implement the different tests.
   The reviewer should write the testing code on their own fork of the repository,
   as outlined in the `GitHub workflow`_.

4. After writing the appropriate code,
   the reviewer should execute it by running the `Unit tests`_.
   The unit tests will produce error messages if
   inconsistencies are found between the original contributor's implementation
   and the tests written by the reviewer.

5. Once the reviewer is confident in their tests,
   they should submit a PR with their test code.
   The reviewer may choose to do so even if some tests fail,
   to facilitate discussion with the original contributor (see Step 6 below).

6. If the tests written by the reviewer produce error messages,
   the differences between the implementation of the original contributor and
   the blind tests of the reviewer need to be resolved through discussion
   between the two of them.
   This discussion can take place either in the **review PR** submitted in Step 5,
   or in the **QC issue** opened in Step 2.
   Differences between the two implementations can indicate an error,
   but very often they are a result of different interpretations of the
   data presented in the source publications.
   For example, there might be different mutation rates estimated for a given species
   from two different groups of samples.
   The original contributor and reviewer should reach an agreement
   as to the best (or at least a reasonable) interpretation of the published data.
   If they cannot reach an agreement,
   then the discussion on GitHub should be opened to others in the community.
   It may also be useful to contact the authors of the original publication
   to resolve some of these ambiguities.
   After each difference is resolved, the final decision should be clearly
   noted in the discussion on GitHub,
   and the code should be modified accordingly.
   This could be either the code written by the original contributor or the
   test code written by the reviewer (or both in some cases).
   Since at this point the **review PR** submitted in Step 5 is still open (not merged),
   then we recommend making the code changes using additional commits in this PR.
   In case the review process found different possible interpretations
   of the published data,
   the rationale behind the final (consensus) interpretation should be clearly
   specified in comments above the relevant block of code.
   This documentation will help future contributors in resolving
   ambiguities in similar cases.

7. Once the **review PR** submitted in Step 5 passes all unit tests,
   it is merged, and the **QC issue** opened in Step 2 is closed.
   **The new code is now officially added to the** ``stdpopsim`` **catalog!**


------------------------------------------------
Implementing tests for the review of new species
------------------------------------------------

The tests associated with the review of a new species
should be written by the reviewer in the ``tests/test_<SPECIES_ID>.py`` file
as part of Step 3 of the review process described above.
Recall that this file was generated by the ``maintenance`` utility, with most
of the tests disabled.
The reviewer should enable all the tests and implement them.
For example, the test for the recombination rates is initialized by the
``maintenance`` utility in the following form:

.. code-block:: python

    @pytest.mark.skip("Recombination rate QC not done yet")
    @pytest.mark.parametrize(["name", "rate"], {}.items())
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).recombination_rate)

When writing the tests for the recombination rates, the reviewer should first
delete the ``@pytest.mark.skip`` line to enable the test.
Then, they should specify inside the ``{ }`` a valid dictionary:
a list of ``key``:``value`` with the name and average
recombination rate for each chromosome.
We provide an example below from *A. aegypti* (see ``tests/test_AedAeg.py``):

.. code-block:: python

    @pytest.mark.parametrize(
        ["name", "rate"],
        {"1": 0.306e-8, "2": 0.249e-8, "3": 0.291e-8, "MT": 0.0}.items(),
    )
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).recombination_rate)


The tests can be executed by running the complete set of `Unit tests`_,
or by invoking only the tests in ``tests/test_<SPECIES_ID>.py``, as follows:

.. code-block:: shell

   $ python -m pytest tests/test_<SPECIES_ID>.py

The tests compare the values specified in the
test file to the values in the ``species.py`` and ``genome_data.py`` files,
and they produce error messages if differences are found.
Differences should be resolved using the general process outlined in
Step 6 of the `Overview of the stdpopsim review process`_.

******************************
Adding a new demographic model
******************************

A demographic model describes ancestral population sizes, split times,
and migration rates.
Misspecification of the model can generate unrealistic patterns of genetic
variation that will affect downstream analyses.
So, having at least one detailed demographic model is recommended for every species.
A given species might have more than one demographic model,
fit from different data or by different methods.

-----------------------------------
What models are appropriate to add?
-----------------------------------
Any model added to ``stdpopsim`` should be based the **published literature**
and a clear citation to the relevant paper(s) should be provided with the model.
The demographic model should include, at a minimum,
a single population with a series of population sizes changes.
Multi-population models typically include other **demographic events**,
such as population splits and changes in the amount of gene flow between populations.
The values of different parameters should be specified in units of "number of individuals"
(for population sizes) and generations (for times).
Sometimes, you will need to convert values published in the literature
to these units by making some assumptions on the mutation rate;
typically the same assumptions made by the study that published the demographic model.


The ``stdpopsim`` :ref:`sec_catalog` also contains a collection of **generic models**,
which are not associated with a certain species and are primarily used for development
and testing of demographic inference methods.
Due to their nature, the rationale for adding such models is different,
and they are also implemented in a slightly different way.
If you wish to contribute a new **generic model**,
then we suggest that you `open a new issue <http://github.com/popgensims/stdpopsim/issues>`__
to discuss your suggestion with others in the community and decide on the best
way to implement your suggestion.

---------------------------------------------
Getting set up to add a new demographic model
---------------------------------------------

If this is your first time implementing a demographic model in ``stdpopsim``, it's a good
idea to take some time browsing the :ref:`sec_catalog`
to see how existing demographic models are coded and documented.
If you have any questions or confusion about formatting or implementing demographic models, please
don't hesitate to `open a new issue <http://github.com/popgensims/stdpopsim/issues>`__.
We're more than happy to answer any questions and help get you up and running.
Before you add any code, be sure to have forked the ``stdpopsim`` repository
and cloned it locally, following the instructions in the `GitHub Workflow`_ section.


All code for a species' demographic models is written in the ``demographic_models.py``
file in that species directory ``stdpopsim/catalog/<SPECIES_ID>/``
(where ``<SPECIES_ID>`` is the six-character identifier of the species;
e.g., CanFam).
If the species does not currently have any demographic model,
then you should add this file to ``stdpopsim/catalog/<SPECIES_ID>/``,
with the following three lines of code:

.. code-block:: python

  import msprime
  import stdpopsim

  _species = stdpopsim.get_species("<SPECIES_ID>")

Furthermore, to ensure that the demographic model(s) are fully incorporated to the
species' code base, you should add the following import to the ``__init__.py`` file
in the species directory:

.. code-block:: python

  from . import demographic_models

----------------------------
Coding the demographic model
----------------------------

The demographic model should be coded in the ``demographic_models.py`` file
by defining a specialized function, which essentially returns
a ``stdpopsim.DemographicModel`` object initialized with the appropriate values.
This function should then be added to the ``_species`` object using the ``add_demographic_model``
function.
We provide below a template block of code for these two operations:

.. code-block:: python

  def _model_func_name():
      return stdpopsim.DemographicModel(
          id=...,
          description=...,
          long_description=...,
          populations=...,
          citations=...,
          generation_time=...,
          mutation_rate=...,
          population_configurations=...,
          migration_matrix=...,
          demographic_events=...,
      )

      _species.add_demographic_model(_model_func_name())

A demographic model is thus defined using ten different attributes.
The first seven attributes are quite straightforward:

* ``id`` (`string`): A unique, short-hand identifier for this demographic model.
  This id contains a short description written in camel case,
  followed by an underscore, and then four characters:
  (1) a digit character specifying the number of sampled populations;
  (2) the first letter of the name of the first author of the publication;
  (3-4) and two digit characters specifying the year the study was published.
  For example, the "Out of Africa" demographic model for humans published by
  Gutenkunst *et al.* (2009) has the ``id`` "OutOfAfrica_3G09".
  See :ref:`sec_development_naming_conventions` for more details.

* ``description`` (`string`): A brief one-line description of the demographic model.

* ``long_description`` (`string`): A more detailed textual description of the model (short paragraph).

* ``populations``: A list of ``stdpopsim.Population`` objects, which have their own
  ``id`` and ``description``. For example, the Thousand Genomes Project Yoruba panel
  could be defined as ``stdpopsim.Population(id="YRI", description="1000 Genomes YRI
  (Yorubans)")``.

* ``citations``: A list of ``stdpopsim.Citation`` objects for the publications
  from which this model was derived.
  The citation object requires author, year, and doi information, and
  a specified reason for citing this model (see `Coding the species parameters`_).
  The reason associated with demographic model citations will typically be
  ``stdpopsim.CiteReason.DEM_MODEL``.

* ``generation_time`` (`double`): The generation time for the species in years.
  The value of this parameter does not directly affect the simulation,
  since the ``SLiM`` and ``msprime`` simulation engines operate in time units of generations.
  The generation time is only currently used to convert time units to years,
  which is useful when comparing among different demographic models.

* ``mutation_rate`` (`double`): The mutation rate assumed during the inference of this demographic
  model (per bp per generation).
  Most demographic inference methods make some assumption about the average genome-wide
  mutation rate.
  These assumptions are sometimes "baked" into the methods,
  and in other cases are just used to convert parameter values from mutation-scale
  to physical scale (number of individuals for population size and generations for times).
  If you are confident that inference did not make any assumption about mutation rate,
  then set the mutation rate of the demographic model to ``None``.
  However, note that this is quite uncommon, so you should make sure this is the case
  before you set the mutation rate to ``None``.

The final three attributes
(``population_configurations``, ``migration_matrix``, and ``demographic_events``)
describe the inferred demographic history that you wish to code.
This history consists of ancestral population size changes,
migration rates, split times, and admixture events.
Note that population size is defined as the number of individuals, regardless of ploidy.
These attributes should be coded using the standard format of ``msprime``.
If this is your first time specifying a demographic model using ``msprime``,
then we highly recommend that you take some time to read through its
`documentation and tutorials <https://tskit.dev/msprime/docs/stable/quickstart.html>`__.


.. note::

   Most published demographic models provide a range of plausible values for each
   parameter of interest.
   In your coded model, you should use some reasonable point estimate,
   such as the value associated with the the maximum likelihood fit,
   or the mean posterior (for Bayesian methods).

------------------------------------
Adding a parameter table to the docs
------------------------------------

The parameters used in the implementation of the demographic model should
also be specified in the docs in a file  ``docs/parameter_tables/<SPECIES_ID>/<MODEL_ID>.csv``,
where ``<SPECIES_ID>`` is the six-character species id,
and ``<MODEL_ID>`` is the ``id`` of the demographic model.
This provides a straightforward documentation and also helps in the review
process (see below).
Each line in this csv file should have the format::

    Parameter Type (units), Value, Description

You may examine csv files currently in  the ``docs/parameter_tables/`` directory
for useful examples.
Once you completed the csv file,
you can check that the documentation was built properly by running
``make`` in the ``docs/`` directory and opening the Catalog page in the
``docs/_build/`` directory.
See `Documentation`_ for more details.



--------------------------------------------------
Testing your demographic model and submitting a PR
--------------------------------------------------

Once you have written the demographic model function in the ``demographic_models.py`` file,
you should test it locally using the development version of ``stdpopsim``.
First, make sure to install the development version of ``stdpopsim`` and its requirements,
by following the :ref:`sec_development_installation` instructions.
Then, check that your new demographic model function has been imported
by executing the following Python code,
where ``<SPECIES_ID>`` is the six-character species id (e.g., HomSap or AraTha):

.. code-block:: python

  import stdpopsim

  species = stdpopsim.get_species("<SPECIES_ID>")
  for x in species.demographic_models:
      print(x.id)


This prints the identifiers (``id``; see above) for all demographic models defined for the species.
You should make sure that the identifier of your newly added model is printed.

.. note::

    If the identifier of your demographic model is not printed,
    make sure that you included the call ``_species.add_demographic_model(_model_func_name())``
    for your newly defined function ``_model_func_name()``
    in the end of the ``demographic_models.py`` file.

    If you are still having trouble, check the
    `GitHub issues <https://github.com/popsim-consortium/stdpopsim/issues?q=is%3Aissue+adding+demographic+model+>`__
    or `open a new issue <https://github.com/popsim-consortium/stdpopsim/issues/new>`__ to get help from others.

After you confirmed that your demographic model was added to the species code,
you should check that you can successfully simulate it with the Python API.
See :ref:`sec_python_tute` for more details.
Finally, once everything looks okay,
you should submit a pull request (PR) with your changes to the code.
See the `GitHub workflow`_ for more details about this process.

At this point, most of your work is done.
**You have officially joined the** ``stdpopsim`` **development team. Welcome!!**
Your model still needs to undergo review by another member (or members)
of the development team before it is fully incorporated into ``stdpopsim``.
This will likely require additional feedback from you,
so, stay tuned for discussion during the review process.

--------------------------------------------------------
Implementing tests for the review of a demographic model
--------------------------------------------------------

After a contributor submits a PR with a new demographic model,
the code undergoes seven steps of review before it
is officially added to ``stdpopsim`` (see `Overview of the stdpopsim review process`_).
In Step 3 of this process, the reviewer writes testing code for the newly
added demographic model.
This is done in file ``stdpopsim/qc/<SPECIES_ID>.py``
(where ``<SPECIES_ID>`` is the six-character identifier of the species).
If this is the first demographic model added for this species,
the reviewer should create this file and add an import
statement for the species to ``stdpopsim/qc/__init__.py``.

The code written by the reviewer in ``stdpopsim/qc/<SPECIES_ID>.py``
should define a function that returns a
``stdpopsim.DemographicModel`` object, parallel to the function defined
by the original contributor of the demographic model (see `Coding the demographic model`_).
After this function is defined, it should be **registered as the QC function** of the
original function by adding this bit of code to ``stdpopsim/qc/<SPECIES_ID>.py``:

.. code-block:: python

  _species.get_demographic_model(_MODEL_ID_).register_qc(_your_review_function())

Where ``_MODEL_ID_`` is the string specified by the original contributor as the
``id`` of the demographic model, and ``_your_review_function()`` is the function
implemented by the reviewer.

The original demographic model and its registered QC model are compared as part of
the ``stdpopsim`` `Unit tests`_.

********************
Adding a genetic map
********************

Some species have sub-chromosomal recombination maps available. They can be added to
`stdpopsim` by creating a new `GeneticMap` object and providing a formatted file
detailing recombination rates to a designated `stdpopsim` maintainer who then uploads
it to AWS. If there is one for your species that you wish to include, create a space
delimited file with four columns: Chromosome, Position(bp), Rate(cM/Mb), and Map(cM).
Each chromosome should be placed in a separate file and with the chromosome id in the
file name in such a way that it can be programatically parsed out. IMPORTANT: chromosome
ids must match those provided in the genome definition exactly! Below is an example start
to a recombination map file (see `here
<https://tskit.dev/msprime/docs/stable/api.html#msprime.RateMap.read_hapmap>`__
for more details)::

    Chromosome Position(bp) Rate(cM/Mb) Map(cM)
    chr1 32807 5.016134 0
    chr1 488426 4.579949 0

Once you have the recombination map files formatted, tar and gzip them into a single
compressed archive. The gzipped tarball must be FLAT (there are no directories in the
tarball). This file will be sent to one of the `stdpopsim` uploaders for placement in the
AWS cloud once the new genetic map(s) are approved. Finally, you must add a `GeneticMap`
object to the file named for your species in the ``stdpopsim/catalog/<SPECIES_ID>/`` directory
(the one that contains all the simulation code for that species,
see `Getting set up to add a new species`_):

.. code-block:: python

    _genetic_map_citation = stdpopsim.Citation(
        doi="FILL_ME", author="FILL_ME", year=9999, reasons={stdpopsim.CiteReason.GEN_MAP}
    )
    """
    The file_pattern argument is a pattern that matches the recombination map filenames,
    where '{id}' is replaced with the 'id' field of a given chromosome.
    """
    _gm = stdpopsim.GeneticMap(
        species=_species,
        id="FILL_ME",  # ID for genetic map, see naming conventions
        description="FILL_ME",
        long_description="FILL_ME",
        url=("https://stdpopsim.s3-us-west-2.amazonaws.com/genetic_maps/dir/filename"),
        sha256="FILL_ME",
        file_pattern="name_{id}_more_name.txt",
        citations=[_genetic_map_citation],
    )

    _species.add_genetic_map(_gm)

The SHA256 checksum of the the genetic map tarball can be obtained using the
``sha256sum`` command from GNU coreutils. If this is not available on your
system, the following can instead be used:

.. code-block:: sh

   python -c 'from stdpopsim.utils import sha256; print(sha256("genetic_map.tgz"))'

Once all this is done, submit a PR containing the code changes and wait for directions
on whom to send the compressed archive of genetic maps to (currently Andrew Kern is the
primary uploader but please wait to send files to him until directed).

**************************
Lifting over a genetic map
**************************
Existing genetic maps will need to be lifted over to a new assembly, if and when the
current assembly is updated in `stdpopsim`. This process can be partially automated by running
the liftOver maintenance code.

First, you must download and install the ``liftOver`` executable from the
`UCSC Genome Browser Store <https://genome-store.ucsc.edu/>`__.
Next, you must download the appropriate chain files, again from UCSC
(see `UCSC Genome Browser downloads
<http://hgdownload.soe.ucsc.edu/downloads.html#liftover>`__ for more details).
To validate the remapping between assemblies it is required to have chain files
corresponding to both directions of the liftOver
(e.g. `hg19ToHg38.over.chain.gz` and `hg38ToHg19.over.chain.gz`) as in the
example below.

An example of the process for
lifting over the `GeneticMap` ``"HapMapII_GRCh37"`` to the ``"Hg19"`` assembly
is shown below:

.. code-block:: sh

    python /maintenance/liftOver_catalog.py \
        --species HomSap \
        --map HapMapII_GRCh37 \
        --chainFile hg19ToHg38.over.chain.gz \
        --validationChain hg38ToHg19.over.chain.gz \
        --winLen 1000 \
        --useAdjacentAvg \
        --retainIntermediates \
        --gapThresh 1000000

Here, the argument ``"--winLen"`` corresponds to the size of the window over which a weighted
average of recombination rates is taken when comparing the original map with the
back-lifted map (for validation purposes only). The argument ``"--gapThresh"`` is used to select a threshold for
which gaps in the new assembly longer than the ``"--gapThresh"`` will be set with a
recombination rate equal to 0.0000, instead of an average rate. The type of average rate used for gaps
shorter than the ``"--gapThresh"`` is determined either by using the mean rate of two most adjacent windows
or by using the mean rate for the entire chromosome, using options ``"--useAdjacentAvg"`` or
``"--useChromosomeAvg"``` respectively.

Validation plots will automatically be generated in the ``liftOver_validation/``
directory. Intermediate files created by the ``liftOver`` executable will be saved
for inspection in the ``"/liftOver_intermediates/"``, only if the
``"--retainInermediates"`` option is used. Once the user has inspected the validation plots
and deemed the liftOver process to be sufficiently accurate, they can proceed to generating
the SHA256 checksum.

The SHA256 checksum of the new genetic map tarball can be obtained using the
``sha256sum`` command from GNU coreutils. If this is not available on your
system, the following can instead be used:

.. code-block:: sh

   python -c 'from stdpopsim.utils import sha256; print(sha256("genetic_map.tgz"))'

The newly lifted over maps will be formatted in a compressed archive and
automatically named using the assembly name from the chain file.
This file will be sent to one of the `stdpopsim` uploaders for placement in the
AWS cloud, once the new map is approved. Finally, you must add a `GeneticMap`
object to the file named for your species in the `stdpopsim/catalog/<SPECIES_ID>/`
directory, as shown in `Adding a genetic map`_.

Again, once all this is done, submit a PR containing the code changes and wait for
directions on whom to send the compressed archive of genetic maps to
(currently Andrew Kern is the primary uploader but please wait to send files
to him until directed).

.. note::

    The ``GeneticMap`` named ``"ComeronCrossoverV2_dm6"`` for ``"DroMel"``
    was generated by similar code (albeit slightly different
    compared to that shown above) using the following command:

.. code-block:: sh

     python /maintenance/liftOver_comeron2012.py \
         --winLen 1000 \
         --gapThresh 1000000 \
         --useAdjacentAvg \
         --retainIntermediates


.. note::

    The ``GeneticMap`` named ``"SalomeAveraged_TAIR10"`` for ``"AraTha"``
    was generated by aligning the TAIR7 and TAIR10 with ``"minimap2"``,
    and lifting the recombination rates on TAIR7 to TAIR10 with
    ``"paftools.js liftover"``.


.. _sec_development_dfe_model:

******************
Adding a DFE model
******************

A distribution of fitness effects (DFE) describes
the probability distribution of selection coefficients
(deleterious, neutral, and beneficial)
for mutations occurring in a certain set of genomic regions.
This is a central component of the way that ``stdpopsim``
incorporates natural selection in its simulations.
See :ref:`sec_simulating_sel`.
There are various computational methods for estimating DFEs from genomic data
and you may use published DFEs to code a DFE model, as described below.

---------------------------------------------
Getting set up to add a new DFE model
---------------------------------------------

If this is your first time implementing a DFE in ``stdpopsim``, it's a good
idea to take some time browsing the :ref:`sec_catalog`
to see how existing DFE models are coded and documented.
If you have any questions or confusion about formatting or implementing demographic models, please
don't hesitate to `open a new issue <http://github.com/popgensims/stdpopsim/issues>`__.
We're more than happy to answer any questions and help get you up and running.
Before you add any code, be sure to have forked the ``stdpopsim`` repository
and cloned it locally, following the instructions in the `GitHub Workflow`_ section.

The code for for a species' DFE models is written in the ``dfes.py``
file in that species directory ``stdpopsim/catalog/<SPECIES_ID>/``
(where ``<SPECIES_ID>`` is the six-character identifier of the species;
e.g., AraTha).
If the species does not currently have any DFE model,
then you should add this file to ``stdpopsim/catalog/<SPECIES_ID>/``,
with the following two lines of code:


.. code-block:: python

  import stdpopsim

  _species = stdpopsim.get_species("<SPECIES_ID>")

Furthermore, to ensure that the DFE model is fully incorporated to the
species' code base, you should add the following import to the ``__init__.py`` file
in the species directory:

.. code-block:: python

  from . import dfes

--------------------
Coding the DFE model
--------------------

The DFE model should be coded in the ``dfes.py`` file
by defining a specialized function, which essentially returns
a ``stdpopsim.DFE`` object initialized with the appropriate values.
This function should then be added to the ``_species`` object using the ``add_dfe``
function.
We provide below a template block of code for these two operations:

.. code-block:: python

  def _dfe_func_name():
      return stdpopsim.DFE(
          id=...,
          description=...,
          long_description=...,
          citations=...,
          mutation_types=...,
          proportions=...,
      )

      _species.add_dfe(_dfe_func_name())

A DFE model is thus defined using six different attributes.

* ``id`` (`string`): A unique, short-hand identifier for this DFE model.
  This id contains a short description of the distribution written in camel case,
  (such as `"LogNormal"` or `"Gamma"`),
  followed by an underscore, and then three characters:
  (1) the first letter of the name of the first author of the publication;
  (2-3) and two digit characters specifying the year the study was published.
  For example, the DFE inferred by Kim *et al.* (2017) has ``id`` set to `"Gamma_K17"`.
  See :ref:`sec_development_naming_conventions` for more details.

* ``description`` (`string`): A brief one-line description of the demographic model.

* ``long_description`` (`string`): A more detailed textual description of the model (short paragraph).

* ``citations``: A list of ``stdpopsim.Citation`` objects for the publications
  from which this model was derived.
  The citation object requires author, year, and doi information, and
  a specified reason for citing this model (see `Coding the species parameters`_).
  The reason associated with demographic model citations will typically be
  ``stdpopsim.CiteReason.DFE``.

* ``mutation_types``: A list of ``stdpopsim.MutationType`` objects corresponding to different
  mutation types (such as negative, neutral, or positive).
  For more details, see the example below and the documentation of :class:`stdpopsim.MutationType`

* ``proportions``: A list of positive numbers that sum to 1 of the same length as ``mutation_types``.
  This list specifies the proportion of each mutation type.

For example, the code block below demonstrates a DFE model
with three mutation types: neutral, negative, and positive.
In this model, negative mutations are assumed to have
dominance coefficience of ``0.5`` and a selection
coefficients distributed according to a Gamma distribution
with mean ``-0.0004`` and shape ``0.27``.
The positive mutations also have a dominance coefficience of ``0.5``,
but they have a fixed selection coefficient of ``0.01``.

.. code-block:: python

    def _dfe_func_name():

        # Default mutation type is neutral
        neutral = stdpopsim.MutationType()
        # Negative mutation type with gamma-distributed selection coefficients
        negative = stdpopsim.MutationType(
            dominance_coeff=0.5,
            distribution_type="g",  # gamma distribution
            distribution_args=[-0.0004, 0.27],  # mean and shape of distributoin
        )
        # Positive mutation type with fixed selection coefficient of 0.01
        positive = stdpopsim.MutationType(
            dominance_coeff=0.5,
            distribution_type="f",  # fixed selection coefficient
            distribution_args=[0.01],  # fixed value
        )

        # The proportions of the three mutation types
        p_neutral = 0.7
        p_negative = 0.299
        p_positive = 1 - p_neutral - p_negative

        return stdpopsim.DFE(
            id=...,
            description=...,
            long_description=...,
            citations=...,
            mutation_types=[neutral, negative, positive],
            proportions=[p_neutral, p_negative, p_positive],
        )


    _species.add_dfe(_dfe_func_name())


------------------------------------------
Testing your DFE model and submitting a PR
------------------------------------------

After you finished your implementation, and specified all the
necessary citations,
we recommend that you run some basic local tests to see that
the model was successfully loaded to ``stdpopsim``.
You may follow the process outlined for `Testing your demographic model and submitting a PR`_.

Once you are convinced that the model was accurately implemented and loaded to ``stdpopsim``,
you should submit a pull request (PR) with your changes to the code.
See the `GitHub workflow`_ for more details about this process.

At this point, most of your work is done.
**You have officially joined the** ``stdpopsim`` **development team. Welcome!!**
Your DFE model still needs to undergo review by another member
of the development team before it is fully incorporated into ``stdpopsim``.
This will likely require additional feedback from you,
so, stay tuned for discussion during the review process.

---------------------
Reviewing a DFE model
---------------------

The review process for DFE models is currently being developed.
For now, we suggest that you
`open a new blank issue <https://github.com/popsim-consortium/stdpopsim/issues/new>`__
and specify the following information:

1. **PR for new model:**
2. **Original paper:**
3. **Parameter values:**
4. **Potential issues:**
5. **QC'er requests:**

A reviewer will be assigned to check your implementation and approve it.
All discussion about the review can be conducted in the **QC issue**
mentioned above.

****************
Coding standards
****************

To ensure that the code in ``stdpopsim`` is as readable as possible
and follows a reasonably uniform style, we require that all code follows
the `PEP8 <https://www.python.org/dev/peps/pep-0008/>`__ style guide.
Lines of code should be no more than 89 characters.
Conformance to this style is checked as part of the Continuous Integration
testing suite.

.. _sec_development_naming_conventions:

******************
Naming conventions
******************

To ensure uniformity in naming schemes across objects in ``stdpopsim``
we have strict conventions for species, genetic maps, and demographic
models.

Species names follow a ``${first_3_letters_genus}${first_3_letters_species}``
convention with capitilization such that Homo sapiens becomes "HomSap". This
is similar to the UCSC Genome Browser naming convention and should be familiar.

Genetic maps are named using a descriptive name and the assembly version according
to ``${CamelCaseDescriptiveName}_${Assembly}``. e.g., the HapMap phase 2 map on
the GRCh37 assembly becomes HapMapII_GRCh37.

Demographic models are named using a combination of a descriptive name,
information about the simulation, and information about the publication it was
presented in. Specifically we use
``${SomethingDescriptive}_${number_of_populations}${first_author_initial}${two_digit_date}``
where the descriptive text is meant to capture something about the model
(i.e. an admixture model, a population crash, etc.) and the number of populations
is the number of populations implemented in the model (not necessarily the number
from which samples are drawn). For author initial we will use a single letter, the 1st,
until an ID collision, in which case we will include the 2nd letter, and so forth.

DFE (Distribution of Fitness Effects) models are similarly named using a string describing
the distribution, and information about the publication:
``${SomethingDescriptive}_${First_authors_last_name_first_letter}{two_digit_date}``.
For instance, if the distribution in question is a lognormal distribution,
then ``LogNormal`` might be the descriptive string.


**********
Unit tests
**********

All code added to ``stdpopsim`` should have
`unit tests <https://en.wikipedia.org/wiki/Unit_testing>`__. These are typically
simple and fast checks to ensure that the code makes basic sense (the
entire unit test suite should not require more than a few seconds to run).
Test coverage is checked using `CodeCov <https://codecov.io/gh/popgensims/stdpopsim>`__,
which generates reports about each pull request.

It is not practical to test the statistical properties of simulation models
as part of unit tests.

The unit test suite is in the ``tests/`` directory. Tests are run using the
`pytest <https://docs.pytest.org/en/stable/>`__ module. Use::

    $ python3 -m pytest

from the project root to run the full test suite. Pytest is very powerful and
has lots of options; please see the `tskit documentation
<https://tskit.dev/tskit/docs/stable/development.html#tests>`__ for help on
how to run pytest and some common options.

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
Code Coverage
*************

As part of the continuous testing suite we have automated checking of how
well the test units cover the source code. As a result it's very helpful
to check locally how well your tests are covering your code by asking
`pytest` for coverage reports. This can be done with::

    $ pytest --cov-report html --cov=stdpopsim tests/

this will output a directory of html files for you to browse test coverage
for every file in `stdpopsim` in a reasonably straightfoward graphical
way. To see them, direct your web browser to the `htmlcov/index.html` file.
You'll be looking for lines of code that are highlighted yellow or red
indicated that tests do not currently cover that bit of code.


*************
Documentation
*************

Documentation is written using `reStructuredText <http://docutils.sourceforge.net/rst.html>`__
markup and the `sphinx <http://www.sphinx-doc.org/en/master/>`__ documentation system.
It is defined in the ``docs/`` directory.

To build the documentation type ``make`` in the ``docs/`` directory. This should build
HTML output in the ``docs/_build/html/`` directory.

.. note::

    You will need ``stdpopsim`` to be installed for the build to work.


********************
Making a new release
********************

Here is a list of things to do when making a new release:

1. Update the changelog and commit
2. Create a release using the GitHub UI
3. `git fetch upstream` on your local branch.
    Then check out `upstream/main` and create a release tarball
    (with `python setup.py sdist`).
    Setuptools_scm will detect the version appopriately.
4. Upload to PyPI: `twine upload dist/{version just tagged}.tar.gz`
5. After the release, if everything looks OK,
   update the symlink for ``stable`` in the
   `stdpopsim-docs <https://github.com/popsim-consortium/stdpopsim-docs>`__
   repository
6. Check on the conda feedstock PR.
