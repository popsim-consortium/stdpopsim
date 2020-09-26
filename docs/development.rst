.. _sec_development:

===========
Development
===========

``stdpopsim`` **is a community effort, and we welcome YOU to join us!**

We envision at least three main types of ``stdpopsim`` developers:

1. Demographic model contributors
2. API developers
3. Documentation and tutorial curators

`Demographic model contributors` add simulation code of published models.
This could be your own published model or any other published model you think
would be useful. This is the main way we envision biologists to continually add
to the catalog of available models, and it is a great first step for new
contributors to learn the ins and outs of ``stdpopsim`` development. See the
sections `Adding a new demographic model`_ and
`Demographic model review process`_ to get started.

`API developers` work on infrastructure development for the PopSim Consortium,
which could include improvements and additions to the internal code base of
``stdpopsim``, establishment of benchmarking pipelines,
and new projects that align with consortium goals.

`Documentation and tutorial curators` help maintain the documentation and tutorials.
This can be as easy as pointing out confusing bits of the documentation in a
GitHub `issue <http://github.com/popgensims/stdpopsim/issues>`_., or adding or editing
the documentation. See the section `Documentation`_.

Get into contact with the ``stdpopsim`` community by subscribing to our email list
`serve <https://lists.uoregon.edu/mailman/listinfo/popgen_benchmark>`_
and by creating and commenting on
Github `issue <http://github.com/popgensims/stdpopsim/issues>`_.
There is a lot of chatter through
`Github <http://github.com/popgensims/stdpopsim>`_, and we’ve been building code
there cooperatively.
If you want to help out and don't know where to start, you can look through the
list of
`Good first issues
<https://github.com/popgensims/stdpopsim/issues?q=is%3Aopen+is%3Aissue+label%3A%22
good+first+issue%22>`_
or
`Help wanted issues
<https://github.com/popgensims/stdpopsim/issues?q=is%3Aopen+is%3Aissue+label%3A%22
help+wanted%22>`_


To get started helping with ``stdpopsim`` development, please read the
following sections to learn how to contribute.
And, importantly, please have a look at our
`code of conduct <https://github.com/popsim-consortium/stdpopsim/blob/master/CODE_OF_CONDUCT.md>`_.

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
<https://msprime.readthedocs.io/en/stable/installation.html>`_ if you
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
Adding a new species
********************
To add a new species to `stdpopsim` several things are required:
1. The genome definition
2. Default species parameters
3. A genetic map with local recombination rates (optional)

Once you have these things the first step is to create a new file in the `catalog`
directory named for the species (see `Naming conventions`_ for more details). All
code described below should go in this file unless explicitly specified otherwise.

--------------------------
Default species parameters
--------------------------

Four default parameters are required to create a new species:
1. Generation time estimate
2. Mutation rate
3. Recombination rate
4. Characteristic population size

These parameters should be based on what values might be drawn from a typical population
as represented in the literature for that species. Consequently one or more citations for
each value are expected and will be required for constructing the species object detailed
below.

-----------------------------------
Adding/Updating a genome definition
-----------------------------------

A genome definition is created with a call to `stdpopsim.Genome()`  which requires a list
of chromosomes and a citation for the assembly. `stdpopsim` has an automated procedure
for obtaining this list from ensembl and saving it for automated parsing. First however
the initial species directory must be created in the `stdpopsim/catalog` directory (e.g.
`stdpopsim/catalog/AraTha`). Once that is done, run the `update_ensembl_data.py` script
present in the top level directory providing the ensembl species id(s) as "_" delimited
name(s) for positional arguments as shown below. If no positional arguments are specified
then all specified registered in `stdpopsim` will be updated.

.. code-block:: shell

    python update_ensembl_data.py arabidopsis_thaliana

This will write/overwrite the `ensembl_info.py` file in the appropriate catalog
subdirectory. Then add the following to the head of `catalog/{species_id}/__init__.py`.

.. code-block:: python

    from . import genome_data

To create the chromosome object that make up a genome add the following code to
`catalog/{species_id}/__init__.py` and supply default mutation and recombination rates
along with citations for the assembly (and additional ones for the mutatation, and
recombination rates if necessary). This is then used to create a `genome` object.

.. code-block:: python

    # A citation for the chromosome parameters. Additional citations may be needed if
    # the mutation or recombination rates come from other sources. In that case create
    # additional citations with the appropriate reasons specified (see API documentation
    # for stdpopsim.citations)

    _assembly_citation = stdpopsim.Citation(
        doi="FILL ME",
        year="FILL ME",
        author="Author et al.",
        reasons={stdpopsim.CiteReason.ASSEMBLY})

    # Parse list of chromosomes into a list of Chromosome objects which contain the
    # chromosome name, length, mutation rate, and recombination rate
    _chromosomes = []

    for name, data in genome_data.data["chromosomes"].items():
        _chromosomes.append(stdpopsim.Chromosome(
            id=name,
            length=data["length"],
            synonyms=data["synonyms"],
            mutation_rate=FILL_ME,
            recombination_rate=FILL_ME
        ))

    # Create a genome object

    _genome = stdpopsim.Genome(
        chromosomes=_chromosomes,
        assembly_citations=[_assembly_citation])

Once you have a genome object you can create a new `Species` object which contains
species identifiers, the genome, and default generation time and population size settings
along with the relevant citation(s). Below is an example species definition for
Arabidopsis thaliana and a final line of code that registers the species in the catalog.

.. code-block:: python

    _gen_time_citation = stdpopsim.Citation(
        doi="https://doi.org/10.1890/0012-9658(2002)083[1006:GTINSO]2.0.CO;2",
        year="2002",
        author="Donohue",
        reasons={stdpopsim.CiteReason.GEN_TIME})

    _pop_size_citation = stdpopsim.Citation(
            doi="https://doi.org/10.1016/j.cell.2016.05.063",
            year="2016",
            author="1001GenomesConsortium",
            reasons={stdpopsim.CiteReason.POP_SIZE})

    _species = stdpopsim.Species(
        id="AraTha",
        name="Arabidopsis thaliana",
        common_name="A. thaliana",
        genome=_genome,
        generation_time=1.0,
        generation_time_citations=[_gen_time_citation],
        population_size=10**4,
        population_size_citations=[_pop_size_citation]
        )

    stdpopsim.register_species(_species)

Once all of this is done, go to the `catalog/__init__.py` file and add a line like the
one below using the six-letter species identifier. Make sure to keep the comment to
prevent linting issues.

.. code-block:: python

    from .catalog import PonAbe  # NOQA

--------------------
Adding a genetic map
--------------------
Some species have sub-chromosomal recombination maps available. They can be added to
`stdpopsim` by creating a new `GeneticMap` object and providing a formatted file
detailing recombination rates to a desginated `stdpopsim` maintainer who then uploads
it to AWS. If there is one for your species that you wish to include, create a space
delimited file with four columns: Chromosome, Position(bp), Rate(cM/Mb), and Map(cM).
Each chromosome should be placed in a seperate file and with the chromosome id in the
file name in such a way that it can be programatically parsed out. IMPORTANT: chromosome
ids must match those provided in the genome definition exactly! Below is an example start
to a recombination map file (see `here
<https://msprime.readthedocs.io/en/stable/api.html#msprime.RecombinationMap.read_hapmap>`_
for more details)::

    Chromosome Position(bp) Rate(cM/Mb) Map(cM)
    chr1 32807 5.016134 0
    chr1 488426 4.579949 0

Once you have the recombination map files formatted, tar and gzip them into a single
compressed archive. The gzipped tarball must be FLAT (there are no directories in the
tarball). This file will be sent to one of the `stdpopsim` uploaders for placement in the
AWS cloud once the new genetic map(s) are approved. Finally, you must add a `GeneticMap`
object to the file named for your species in the `catalog` directory (the same one in
which the genome is defined) as shown below:

.. code-block:: python

    _genetic_map_citation = stdpopsim.Citation(
            doi="FILL_ME",
            author="FILL_ME",
            year=9999,
            reasons={stdpopsim.CiteReason.GEN_MAP})
    """
    The file_pattern argument is a pattern that matches the recombination map filenames,
    where '{id}' is replaced with the 'id' field of a given chromosome.
    """
    _gm = stdpopsim.GeneticMap(
        species=_species,
        id="FILL_ME", # ID for genetic map, see naming conventions
        description="FILL_ME",
        long_description="FILL_ME",
        url=("https://stdpopsim.s3-us-west-2.amazonaws.com/genetic_maps/dir/filename"),
        file_pattern="name_{id}_more_name.txt",
        citations=[_genetic_map_citation])

    _species.add_genetic_map(_gm)

Once all this is done, submit a PR containting the code changes and wait for directions
on whom to send the compressed archive of genetic maps to (currently Andrew Kern is the
primary uploader but please wait to send files to him until directed).

******************************
Adding a new demographic model
******************************

Steps for adding a new demographic model:

1. `Fork the repository and create a branch`_
2. `Write the model function in the catalog source code`_
3. `Write parameter table`_
4. `Test the model locally`_
5. `Submit a Pull Request on GitHub`_

If this is your first time implementing a demographic model in `stdpopsim`, it's a good
idea to take some time browsing the
`Catalog <https://stdpopsim.readthedocs.io/en/latest/catalog.html>`_
and species' demographic models in the
source code to see how existing models are typically written and documented. If you have
any questions or confusion about formatting or implementing demographic models, please
don't hesitate to open an `issue <http://github.com/popgensims/stdpopsim/issues>`_ --
we're more than happy to answer any questions and help get you up and running.

-----------------------------------
What models are appropriate to add?
-----------------------------------

`Stdpopsim` supports any demographic model from the published literature that gives
enough information to be able to define `msprime` demography objects. At a minimum, that
includes population sizes and the timing of demographic events. These values need to
either be given in "physical" units (that is, raw population sizes and time units in
generations), or be able to be converted to physical units using, e.g., mutation rates
used in the published study.

Note that it is not necessary that the demographic model is attached to a particular
species. `Stdpopsim` contains a collection of generic models that are widely used in
developing and testing inference methods. If there is a generic model that does not
currently exist in our catalog but would be useful to include, we also welcome those
contributions. Again, you should provide a citation for a generic models, or it
should be commonly used.

---------------------------------------
Fork the repository and create a branch
---------------------------------------

Before implementing any model, be sure to have forked the `stdpopsim` repository
and cloned it locally, following the instructions in the `GitHub Workflow`_ section.
Models are first implemented and tested locally, and then submitted as a pull request
to the `stdpopsim` repository, at which point it is verified by another developer
before being fully supported within `stdpopsim`.

---------------------------------------------------
Write the model function in the catalog source code
---------------------------------------------------

In the ``stdpopsim`` catalog source code (found in ``stdpopsim/catalog/``),
each species has a module that defines all of the necessary functions to run
simulations for that species, including the demographic model. In each species module,
you will see that each type of function divided by comments, such as::

    ###########################################################
    #
    # Demographic models
    #
    ###########################################################

Go to the ``Demographic models`` section of the source code.
The demographic model function should follow this format:

.. code-block:: python

    def _model_func_name():
        id = "FILL ME"
        description = "FILL ME"
        long_description = """
        FILL ME
        """
        populations = [
            stdpopsim.Population(id="FILL ME", description="FILL ME"),
        ]
        citations = [
            stdpopsim.Citation(
                author="FILL ME",
                year="FILL ME",
                doi="FILL ME",
                reasons={stdpopsim.CiteReason.DEM_MODEL})
        ]

        generation_time = "FILL ME"

        # parameter value definitions based on published values

        return stdpopsim.DemographicModel(
            id=id,
            description=description,
            long_description=long_description,
            populations=populations,
            citations=citations,
            generation_time=generation_time,
            population_configurations=[
            "FILL ME"
            ],
            migration_matrix=[
            "FILL ME"
            ],
            demographic_events=[
            "FILL ME"
            ],
            )


    _species.add_demographic_model(_model_func_name())


The demographic model should include the following:

* ``id``: A unique, short-hand identifier for this demographic model. This ``id``
  contains a short description written in camel case, followed by an underscore, and then
  four characters (the number of sampled populations, the first letter of the name of the
  first author, and the year the study was published). For example, the Gutenkunst et al.
  (2009) Out of Africa demographic model has the ``id`` "OutOfAfrica_3G09". See
  `Naming conventions`_ for more details.
* ``description``: A brief one-line description of the demographic model.
* ``long_description``: A longer description (say, a concise paragraph) that describes
  the model in more detail.
* ``populations``: A list of ``stdpopsim.Population`` objects, which have their own
  ``id`` and ``description``. For example, the Thousand Genomes Project Yoruba panel
  could be defined as ``stdpopsim.Population(id="YRI", description="1000 Genomes YRI
  (Yorubans)")``.
* ``citations``: A list of ``stdpopsim.Citation`` objects for the appropriate citation
  for this model. The citation object requires author, year, and doi information, and
  a specified reason for citing this model.
* ``generation_time``: The generation time for the species in years. If you are
  implementing a generic model, the generation time should default to 1.


Every demographic model has a few necessary features or attributes. First of all,
demographic models are defined by the population sizes, migration rates, split and
admixture times, and generation lengths given in the source publication. We often take
the point estimates for each of the values from the best fit model (for example, the
parameters that give the maximum likelihood fit), which are translated into
`msprime`-formatted demographic inputs.

`Msprime`-defined demographic models are specified through the
``population_configurations``, ``migration_matrix``, and ``demographic_events``. If this
is your first time specifying a model using `msprime`, it's worth taking some time to
read through the `msprime`
`documentation and tutorials <https://msprime.readthedocs.io/en/stable/tutorial.html>`_.


---------------------
Write parameter table
---------------------

The parameters used in the implementation must
also be listed in a csv file in the ``docs/parameter_tables`` directory. This ensures
that the documentation for this model displays the parameters.

Take a look at the csv files currently in ``docs/parameter_tables`` for inspiration.
The csv file should have the format::

    Parameter Type (units), Value, Description


We can check that the documentation builds properly after implementation by running
``make`` in the docs directory and opening the Catalog page from the ``docs/_build/``
directory. See `Documentation`_ for more details.


----------------------
Test the model locally
----------------------

Once you have written the demographic model function, you should test the model locally
with ``stdpopsim``. Follow the development :ref:`sec_development_installation`
instructions to install the development ``stdpopsim`` version along with the
requirements.

Now check that your new demographic model function has been imported:

.. code-block:: python

    import stdpopsim
    species = stdpopsim.get_species("HomSap")
    for x in species.demographic_models:
        print(x.id)

    # OutOfAfrica_3G09
    # OutOfAfrica_2T12
    # Africa_1T12
    # AmericanAdmixture_4B11
    # OutOfAfricaArchaicAdmixture_5R19
    # Zigzag_1S14
    # AncientEurasia_9K19
    # PapuansOutOfAfrica_10J19


The example above lists the imported demographic models for humans.
You should substitute ``"HomSap"`` for which ever species you added your model to.
Your new model should be printed along with currently available demographic models.

.. note::

    If your demographic model does not print, after defining your model function,
    did you include the call ``_species.add_demographic_model(_model_func_name())``,
    where ``_model_func_name()`` is your model function name?

    If you are still having trouble, check the
    `GitHub issues <https://github.com/popsim-consortium/stdpopsim/issues?q=is%3Aissue+adding+demographic+model+>`_,
    or `open an issue <https://github.com/popsim-consortium/stdpopsim/issues/new>`_.

Next, check that you can successfully run a simulation with your new model with the
Python API. See :ref:`sec_python_tute` for more details.

-------------------------------
Submit a Pull Request on GitHub
-------------------------------

Once you have implemented the demographic model locally, including
documentation, the next step is to open a pull request with this addition.
See the `GitHub workflow`_ for more details.

---------------------------------------
So the model is implemented. What next?
---------------------------------------

Now at this point, most of your work is done!  The model is reviewed and
verified following the `Demographic model review process`_ by an independent member
of the development team, and there may be some discussion about formatting and
to clear up any confusing bits of the demographic parameters before the model is
fully incorporated into `stdpopsim`.

Thank you for your contribution, and welcome to the `stdpopsim` development team!

********************************
Demographic model review process
********************************

When Developer A creates a new demographic model on their local fork they must
follow these steps for it to be officially supported by stdpopsim:

    1. Developer A submits a PR to add a new model to the catalog. This includes
       full documentation (i.e., the documentation that will be
       rendered on rtd). The code is checked for any obvious problems/style
       issues etc by a maintainer and merged when it meets these basic
       standards. The new catalog model is considered 'preliminary'.

    2. Developer A creates an `issue
       <https://github.com/popsim-consortium/stdpopsim/issues/new/choose>`__
       tracking the QC for the model which includes information about the
       primary sources used to create the model and the population indices
       used for their msprime implementation. To create a new Model QC issue,
       click "New issue" from the "Issues" tab on GitHub, and click "Get
       started" to use the Model QC issue template. Follow the template to
       include the necessary information in the issue. Developer B is then
       assigned/volunteers to do a blind implementation of the model.

    3. Developer B creates a blind implementation of the model in the
       ``stdpopsim/qc/species_name_qc.py`` file, remembering to register the
       QC model implementation (see other QC models for examples).  Note that
       if you are adding a new species you will have to add a new import to
       ``stdpopsim/qc/__init__.py``.

    4. Developer B runs the units tests to verify the equivalence of the
       catalog and QC model implementations.

    5. Developer B then creates a PR, and all being good, this PR is merged and
       the QC issue is closed.

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

Finally demographic models are named using a combination of a descriptive name,
information about the simulation, and information about the publication it was
presented in. Specifically we use
``${SomethingDescriptive}_${number_of_populations}${first_author_initial}${two_digit_date}``
where the descriptive text is meant to capture something about the model
(i.e. an admixture model, a population crash, etc.) and the number of populations
is the number of populations implemented in the model (not necessarily the number
from which samples are drawn). For author initial we will use a single letter, the 1st,
until an ID collision, in which case we will include the 2nd letter, and so forth.


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

.. note::

    You will need ``stdpopsim`` to be installed for the build to work.

