.. _sec_tutorial:

========
Tutorial
========
There are two main ways of accessing the resources of the ``stdpopsim`` package
that will be detailed in this tutorial. The first is via the command line
interface (CLI). This is useful if you want to do a straightforward run of the
models in the ``stdpopsim`` :ref:`catalog <sec_catalog>`. The other way to
access the ``stdpopsim`` resources is via the python API. This is a bit more
complicated however it allows for more advanced tasks. This tutorial will walk
through both ways of using the ``stdpopsim`` package as well as some more
advanced tasks you may wish to do.

.. sec_cli_tute:

****************************************************
Using the ``stdpopsim`` Command-Line Interface (CLI)
****************************************************
In order to use the ``stdpopsim`` CLI the ``stdpopsim`` package must be
installed (see :ref:`Installation <sec_installation>`). The CLI provides access
to the :ref:`catalog <sec_catalog>` of models that have already been implemented
by ``stdpopsim``. The first step for using the CLI is to select the species that
you are interested in simulating data for. In order to see which species are
available run

.. code-block:: console

    $ stdpopsim --help

This shows the species currently supported by ``stdpopsim``. This means that
``stdpopsim`` knows various traits of these species including chromosome size
and recombination rates. Once we've selected a species, in this case humans, we
can look at the help again as follows.

.. code-block:: console

    $ stdpopsim homsap --help

For conciseness we do not show the output here but this time you should see a
different output which shows options for performing the simulation itself and
the species default parameters. This includes selecting the demographic model,
chromosome, recombination map, and number of samples. 

The most basic simulation we can run is to draw two samples using the species'
defaults as seen in the species help (``stdpopsim homsap --help``). These
defaults include constant size population, uniform recombination map, and a
default chromosome. To save time we will also specify that the simulation use
chromosome 22 with the ``-c`` option. We also specify that the resulting
tree-sequence formated output should be written to the file ``foo.ts`` with the
``-o`` option. For more information on how to use tree-sequence files see
`tskit <https://tskit.readthedocs.io/en/latest/>`_ .

.. code-block:: console

    $ stdpopsim homsap -c chr22 -o foo.ts 2

.. warning:: It's important to remember to either redirect the output of ``stdpopsim``
                to file or to use the ``-o/--output`` option. If you do not, the 
                binary output may mess up your terminal session.

Say we want to use a specific demographic model. We look up the available models
using the ``--help-models`` flag:

.. code-block:: console

    $ stdpopsim homsap --help-models

This gives all of the possible demographic models we could simulate. We choose
the the two population out-of-Africa model from `Tennesen et al. (2012)
<https://doi.org/10.1126/science.1219240>`_ . By looking at the model help we
find that the name for this model is ``ooa_2`` and that we can specify it using
the ``--model`` option. We choose to draw two samples from the African
population and three samples from the European population. To increase
simulation speed we can also chose to simulate a sequence a fraction of the
length of the specified chromosome using the ``-l`` option (e.g. 5%). This is
just specifying a sequence length, not actually selecting a subset of the
chromosome to sequence and as such cannot be used with anything other than a
uniform recombination map. The command now looks like this:

.. code-block:: console

    $ stdpopsim homsap -c chr22 -l 0.05 -o foo.ts --model ooa_2 2 3

Note that there are now two numbers after the model option. This is because the
model simulates two populations so we have to specify a number of samples to
draw from each population (the order of those numbers is the same as the order
specified in the model documentation). In this case, we are simulating two
African American samples and three European American samples.

Now we want to add an empirical recombination map to make the simulation more
realistic. We can run ``stdpopsim homsap --help-genetic-maps`` to view the
available recombination maps. In this case we choose the HapmapII map. Empirical
recombination maps cannot be used with length multipliers so we have to remove
the ``-l`` option. (NOTE: this may a minute or so to run).

.. code-block:: console

    $ stdpopsim homsap -g HapmapII_GRCh37 -c chr22 -o foo.ts --model ooa_2 2 3

For reproducibility we can also choose set seed for the simulator using the
``-s`` flag.

.. code-block:: console

    $ stdpopsim homsap -s 1046 -g HapmapII_GRCh37 -c chr22 -o foo.ts --model ooa_2 2 3

Lastly, the CLI also outputs the relevant citations for both the simulator used
and the resources used for simulation scenario.

.. sec_python_tute:

*****************************
The Python interface
*****************************

--------------------------
Running a prexisting model
--------------------------

.. todo::
    Add up to date code

.. TODO port these old examples.

.. ********
.. Examples
.. ********

.. This content should go into a tutorial or somewhere else, but for now it's
.. useful to keep a couple of short examples here for reference and
.. to motivate the API.

.. Models and genome data are split by species. Information about the genomes
.. of a particular species is held in the ``genome`` class variable. So,
.. suppose we wish to perform a simple simulation of human chromosome 22, we
.. might have:

.. .. code-block:: python

..     import msprime
..     from stdpopsim import homo_sapiens

..     chrom = homo_sapiens.genome.chromosomes["chr22"]
..     ts = msprime.simulate(
..         sample_size=10,
..         recombination_rate=chrom.default_recombination_rate,
..         mutation_rate=chrom.default_mutation_rate,
..         length=chrom.length)

.. (This should be a very quick simulation, and the result will have very few
.. variants, because although it performs a coalescent simulation of
.. a 51,304,566bp chromosome, it does this with the effective population size of
.. `Ne=1`.)

.. The chromosome definitions also aware of recombination maps,
.. which must be first downloaded. The default for ``homo_sapiens`` is ``HapmapII_GRCh37``,
.. which we can find out, then download it as follows
.. (the maps are stored in your ``~/.cache/stdpopsim/`` directory):

.. .. code-block:: python

..    homo_sapiens.genome.default_genetic_map
..    # 'HapmapII_GRCh37'
..    gmap = homo_sapiens.HapmapII_GRCh37()
..    gmap.download()


.. After this has been done (once only), we can run simulations using this genetic map as follows:

.. .. code-block:: python

..     chrom = homo_sapiens.genome.chromosomes["chr22"]
..     ts = msprime.simulate(
..         sample_size=10,
..         mutation_rate=chrom.default_mutation_rate,
..         recombination_map=chrom.recombination_map())

.. Recombination maps will be downloaded on demand and cached in a
.. platform-appropriate user cache directory (e.g., ``$HOME/.cache/stdpopsim`` on
.. Linux). In this example we didn't specify which recombination map we want, and so the
.. API will use the default. We can also ask for specific maps, if we want:

.. .. code-block:: python

..     chrom = homo_sapiens.genome.chromosomes["chr22"]
..     ts = msprime.simulate(
..         sample_size=10,
..         mutation_rate=chrom.default_mutation_rate,
..         recombination_map=chrom.recombination_map("HapmapII_GRCh37"))


.. Demographic models can also be used. For example

.. .. code-block:: python

..     import stdpopsim
..     from stdpopsim import homo_sapiens

..     chrom = homo_sapiens.genome.chromosomes["chr22"]
..     model = homo_sapiens.GutenkunstThreePopOutOfAfrica()
..     # One sample each from YRI, CEU and CHB.
..     samples = [msprime.Sample(population=j, time=0) for j in range(3)]
..     ts = msprime.simulate(
..         samples=samples,
..         recombination_map=chrom.recombination_map(),
..         mutation_rate=chrom.default_mutation_rate,
..         **model.asdict())

.. (This simulation now has a realistic effective population size,
.. so will produce thousands of variant sties, but still runs very fast.)



.. _sec_tutorial_generic_models:

**************
Generic models
**************

The previous sections showed how to simulate population genetic models that have
been published for particular species. It is also sometimes useful to simulate
more generic models. We can still use the information for particular species to
do this; the only difference is that we instantiate the models directly rather
than retrieving them from the :ref:`catalog <sec_catalog>`.


.. code-block:: python

    species = stdpopsim.get_species("homsap")
    contig = species.get_contig("chr22", length_multiplier=0.1)
    model = stdpopsim.PiecewiseConstantSize(species.population_size)
    samples = model.get_samples(10)
    engine = stdpopsim.get_default_engine()
    ts = engine.simulate(model, contig, samples)


Here, we simulate 10% of human chromosome 22 under a constant size
population model, using the current best estimate of the human
effective population size from the **TODO add crossref to catalog section**

.. :ref:`sec_catalog_homo_sapiens_genome`




