.. _sec_tutorial:

=========
Tutorials
=========
There are two main ways of accessing the resources of the ``stdpopsim`` package
that will be detailed in this tutorial. The first is via the command line
interface (CLI). This is useful if you want to do a straightforward run of the
models in the ``stdpopsim`` :ref:`catalog <sec_catalog>`. The other way to
access the ``stdpopsim`` resources is via the python API. This is a bit more
complicated however it allows for more advanced tasks. This tutorial will walk
through both ways of using the ``stdpopsim`` package as well as some more
advanced tasks you may wish to do.

*********************
Running ``stdpopsim``
*********************

.. sec_cli_tute:

The command-line interface (CLI)
********************************

In order to use the ``stdpopsim`` CLI the ``stdpopsim`` package must be
installed (see :ref:`Installation <sec_installation>`). The CLI provides access
to the :ref:`catalog <sec_catalog>` of models that have already been implemented
by ``stdpopsim``. The first step for using the CLI is to select the species that
you are interested in simulating data for. In order to see which species are
available run

.. command-output:: stdpopsim --help

This shows the species currently supported by ``stdpopsim``. This means that
``stdpopsim`` knows various traits of these species including chromosome size
and recombination rates. Once we've selected a species, in this case humans, we
can look at the help again as follows.

.. command-output:: stdpopsim homsap --help
    :ellipsis: 20

For conciseness we do not show all the output here but this time you should see a
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

.. command-output:: stdpopsim homsap --help-models
    :ellipsis: 20

This gives all of the possible demographic models we could simulate. We choose
the two population out-of-Africa :ref:`model <sec_catalog_homsap_models_outofafrica_2t12>`
from `Tennesen et al. (2012) <https://doi.org/10.1126/science.1219240>`_ .
By looking at the model help we
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
available recombination maps. In this case we choose the
:ref:`sec_catalog_homsap_genetic_maps_hapmapii_grch37` map. Empirical
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

.. code-block:: none

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


************************************
Example analyses with ``stdpopsim``
************************************

.. sec_tute_divergence:

Calculating genetic divergence
******************************

In this tutorial, we will simulate some samples of human chromosomes
from different populations,
and then estimate the genetic divergence between each population pair.

-------------------------
1. Simulating the dataset
-------------------------

First, let's use the ``--help-models`` option to see the selection of demographic
models available to us:

.. command-output:: stdpopsim homsap --help-models
    :ellipsis: 20

This prints detailed information about all of the available models to
the terminal.
In this tutorial, we will use the model of African-American admixture from
`2011 Browning et al <http://dx.doi.org/10.1371/journal.pgen.1007385>`_.
From the help output (or the :ref:`Catalog <sec_catalog>`),
we can see that this model has id ``america``,
and allows samples to be drawn from 4 contemporary populations representing African,
European, Asian and African-American groups.

Using the ``--help-genetic-maps`` option, we can also see what recombination maps
are available:

.. command-output:: stdpopsim homsap --help-genetic-maps

Let's go with ``HapmapII_GRCh37``.
The next command simulates 4 samples of chromosome 1 from each of the four
populations, and saves the output to a file called ``afr-america-chr1.trees``.
For the purposes of this tutorial, we'll also specify a random seed using the
``-s`` option.
(Note: This took around 8 minutes to run on a laptop.)

.. code-block:: console

    $ stdpopsim homsap -c chr1 -o afr-america-chr1.trees -s 13 -g HapmapII_GRCh37\
    --model america 4 4 4 4

--------------------------
2. Calculating divergences
--------------------------

We should now have a file called ``afr-america-chr1.trees``.
Our work with ``stdpopsim`` is done; we'll now switch to a Python console and import
the ``tskit`` package to load and analyse this simulated tree sequence file.

.. code-block:: python

    >>> import tskit
    >>> ts = tskit.load("afr-america-chr1.trees")

Recall that `genetic divergence` is the probability that two randomly sampled
chromosomes differ at a nucleotide base.
For a given pair of populations, a pair-specific divergence value is obtained
by randomly sampling one chromosome from each population.
These quantities can be estimated directly from our sample using tskit's
inbuilt :meth:`tskit.TreeSequence.diversity` method.

By looking at
`the documentation <https://tskit.readthedocs.io/en/latest/python-api.html#tskit.TreeSequence.divergence>`_
for this method, we can see that we'll need two inputs: ``sample_sets`` and
``indexes``.
Let's think about what these inputs are, and how we can obtain them with
Python commands.
In our case, the sample sets correspond to the lists
of sample chromosomes (nodes) from each separate population.
We can obtain the necessary list of lists like this:

.. code-block:: python

    >>> sample_list = []
    >>> for pop in range(0, ts.num_populations):
    ...     sample_list.append(ts.samples(pop).tolist())
    >>> print(sample_list)
    [[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11], [12, 13, 14, 15]]

Note that the samples with node IDs 0 - 3 are from population 0,
samples with node IDs 4 - 7 are from population 1 and so on.

The indexes are the pairs of integer indexes corresponding to the populations
that we wish to compare.
We can do this quickly with the ``itertools`` module:

.. code-block:: python

    >>> import itertools
    >>> inds = itertools.combinations_with_replacement(range(0, ts.num_populations), 2)
    >>> inds = list(inds)
    >>> print(inds)
    [(0, 0), (0, 1), (0, 2), (0, 3), (1, 1), (1, 2), (1, 3), (2, 2), (2, 3),
     (3, 3)]

We are now ready to calculate the genetic divergences.

.. code-block:: python

    >>> divs = ts.divergence(sample_sets=sample_list, indexes=inds)
    >>> print(divs)
    array([0.00035424, 0.0003687 , 0.00036707, 0.0003705 , 0.00026696,
       0.00029148, 0.00029008, 0.00025767, 0.0002701 , 0.00028184])

---------------------------
3. Plotting the divergences
---------------------------

The output lists the divergences of all population pairs that are specified in
``indexes``, in the same order.
However, instead of simply printing these values to the console, it might be nicer
to create a heatmap of the values.
Here is some (more advanced) code that does this.
It relies on the ``numpy``, ``seaborn`` and ``matplotlib`` packages.

.. code-block:: python

    >>> import numpy as np
    >>> import seaborn
    >>> import matplotlib.pyplot as plt
    >>> import matplotlib.ticker as ticker
    >>> div_matrix = np.zeros((ts.num_populations, ts.num_populations))
    >>> for pair in range(0, len(inds)):
    ...     pop0, pop1 = inds[pair]
    ...     div_matrix[pop0, pop1] = divs[pair]
    ...     div_matrix[pop1, pop0] = divs[pair]
    >>> seaborn.heatmap(div_matrix, vmin=0, vmax=0.0005, square=True)
    >>> ax = plt.subplot()
    >>> plt.title("Genetic divergence")
    >>> plt.xlabel("Populations", fontweight="bold")
    >>> plt.ylabel("Populations", fontweight="bold")
    >>> ax.set_xticks([0,1,2,3], minor=True)
    >>> ax.set_xticklabels(['AFR', 'EUR', 'ASI', 'ADM'], minor=False)
    >>> ax.tick_params(which='minor', length=0)
    >>> ax.set_yticks([0,1,2,3], minor=True)
    >>> ax.set_yticklabels(['AFR', 'EUR', 'ASI', 'ADM'], minor=False)
    >>> ax.tick_params(which='minor', length=0)

.. image:: _static/tute-divergence.png
    :width: 400px
    :align: center
    :height: 265px
    :alt: Heatmap of divergence values.

These values make sense given the model of demography we have specified:
the highest divergence estimates were obtained when African samples were
compared with samples from other populations, and the lowest divergence
estimates were obtained when Asian samples were compared with themselves.
However, the overwhelming sameness of the sample chromosomes is also evident:
on average, any two sample chromosomes differ at less than 0.04% of positions,
regardless of the populations they come from.
