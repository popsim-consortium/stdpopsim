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


.. _sec_tutorial_generic_models:

-------------------------------------------------
Running a generic model and outputting a vcf file
-------------------------------------------------

In this example, we will use the ``stdpopsim`` API to simulate a generic
model for a particular species. We will use the information for a particular
species and instantiate the model directly.

Here, we will simulate 10% of human chromosome 22 under a constant size
population model, using the current best estimate of the human
effective population size from the :ref:`sec_catalog`.


1. Import the necessary packages

.. code-block:: python

    >>> import stdpopsim

2. Get the particular species information. In this case, we are using
`Homo sapiens`, which has the id "homsap".
But, you could use any species from the :ref:`sec_catalog`.

.. code-block:: python

    >>> species = stdpopsim.get_species("homsap")

3. Set the contig length. We are simulating 0.1 x chromosome 22,
which is about 5Mb. Again, you could use a fraction of any of the
chromosomes listed in the :ref:`sec_catalog`, keeping in mind that
larger contigs will take longer to simulate.

.. code-block:: python

    >>> contig = species.get_contig("chr22", length_multiplier=0.1)

4. Set the model as the generic piecewise constant size model, using the
predefined human effective population size (see :ref:`sec_catalog`).
Since we are providing one effective population size, the model is constant
population size for one population over time.

.. code-block:: python

    >>> model = stdpopsim.PiecewiseConstantSize(species.population_size)

5. Set the number of samples and set the simulation engine.
In this case we will simulate 10 samples and use the default simulator,
which is `msprime`. But, you can go crazy with the sample size!
`msprime` is great at simulating large samples!

.. code-block:: python

    >>> samples = model.get_samples(10)
    >>> engine = stdpopsim.get_default_engine()

6. Simulate the model with the contig length and number of samples we defined above.
We capture the simulation results in a tree sequence object
(:class:`tskit.TreeSequence`).

.. code-block:: python

    >> ts = engine.simulate(model, contig, samples)

7. We can now do some simple checks that our simulation worked with
`tskit
<https://tskit.readthedocs.io>`__.

.. code-block:: python

    >>> ts.num_samples
    10
    >>> ts.num_populations
    1
    >>> ts.num_mutations
    6197
    >>> ts.num_trees
    6863

As expected, there are 10 samples in one population. We can also see that 6197 mutations
and 6863 trees were simulated (since we are not using a seed here, the number of mutations
and trees will be slightly different for you). Try running the simulation again, and notice
that the number of samples and populations stays the same, while the number of mutations
and trees changes.

8. In addition to working directly with the simulated tree squence, we can also output
other common formats used for population genetics analyses.
We can use ``tskit`` to convert the tree sequence to a vcf file called "foo.vcf".
See the tskit documentation (:meth:`tskit.TreeSequence.write_vcf`) for more information.

.. code-block:: python

    >>> with open("foo.vcf", "w") as vcf_file:
    >>>    ts.write_vcf(vcf_file)

Taking a look at the vcf file, we see something like this:

.. code-block::

    ##fileformat=VCFv4.2
    ##source=tskit 0.2.2
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##contig=<ID=1,length=5130457>
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	tsk_0	tsk_1	tsk_2	tsk_3	tsk_4	tsk_5	tsk_6	tsk_7	tsk_8	tsk_9
    1	96	.	0	1	.	PASS	.	GT	0	0	1	0	1	0	0	0	1	0
    1	129	.	0	1	.	PASS	.	GT	0	0	0	0	0	0	0	0	1	0
    1	436	.	0	1	.	PASS	.	GT	0	0	0	0	0	1	0	0	0	0
    1	466	.	0	1	.	PASS	.	GT	0	0	1	0	1	0	0	0	0	0
    1	558	.	0	1	.	PASS	.	GT	0	0	0	0	0	0	0	0	1	0
    1	992	.	0	1	.	PASS	.	GT	1	1	0	1	0	1	1	1	0	1


