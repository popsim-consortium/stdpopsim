.. _sec_installation:

============
Installation
============

You can get started right away, without installing ``stdpopsim`` locally, by using the
`Jupyter Binder <https://mybinder.org/v2/gh/popsim-consortium/stdpopsim/main?filepath=stdpopsim_example.ipynb>`_. |binder|

.. |binder| image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/popsim-consortium/stdpopsim/main?filepath=stdpopsim_example.ipynb


There are two basic options for installing ``stdpopsim`` and its dependencies: either
using :ref:`sec_installation_conda` or :ref:`sec_installation_pip`.
We recommend conda for most users (particularly those using OSX or Windows),
although pip can be more convenient in certain cases.

.. _sec_installation_requirements:

************
Requirements
************

Library requirements for ``stdpopsim`` should be installed automatically when using
either conda or pip.

``stdpopsim`` requires Python 3.5 or later.


.. _sec_installation_conda:

*****
Conda
*****

Pre-built binary packages for ``stdpopsim`` are available through
`conda <https://conda.io/docs/>`_, and built using `conda-forge <https://conda-forge.org/>`_.
Packages for recent version of Python are available for Linux, OSX and Windows. Install
using::

    $ conda install -c conda-forge stdpopsim


+++++++++++
Quick Start
+++++++++++

1. Install ``conda`` using `miniconda <https://conda.io/miniconda.html>`_.
   Make sure you follow the instructions to fully activate your ``conda``
   installation!
2. Set up the `conda-forge channel <https://conda-forge.org/>`_ using
   ``conda config --add channels conda-forge``.
3. Install stdpopsim: ``conda install stdpopsim``.
4. Try it out: ``stdpopsim --version``.


There are several different ways to obtain ``conda``. Please see the
`anaconda installation documentation <https://docs.anaconda.com/anaconda/install/>`_
for full details.


.. _sec_installation_pip:

***
Pip
***

Installing using pip is usually as simple as::

    $ python -m pip install stdpopsim --user

This will install ``stdpopsim`` into your local user Python packages
(on some systems you will need to use ``python3`` rather than
``python``). Please see the Python `package installation
<https://packaging.python.org/tutorials/installing-packages/>`_
tutorial for more details on the various installation options. In particular,
we recommend using a `virtual environment
<https://packaging.python.org/tutorials/installing-packages/#creating-virtual-environments>`_
to improve reproducibility.

It may also be necessary to update your PATH to run the command
line interface. See `here
<https://packaging.python.org/tutorials/installing-packages/#installing-to-the-user-site>`_
for details on what this means and how to do it.

We use `msprime <https://tskit.dev/msprime>`_ as the
default simulation engine, which has some system level dependencies
and requires a functioning compiler. Please see the msprime
`installation documentation
<https://tskit.dev/msprime/docs/stable/installation.html>`_ for
instructions if you encounter errors during installation.

.. _sec_installation_running_cli:

***************
Running the CLI
***************

After installation is complete it should be possible to run the
:ref:`command line interface <sec_cli_args>`. This can be done in one
of two ways:

1. The most reliable way is to use ::

       $ python -m stdpopsim

   Once the ``python`` executable is the same one as was used when installing
   ``conda`` or ``pip``, this is guaranteed to work.

2. It is also possible to run ``stdpopsim`` like a regular Unix program
   using::

    $ stdpopsim

   However, this requires that your `PATH <https://en.wikipedia.org/wiki/PATH_(variable)>`_
   environment variable contains the directory where conda or pip installed the
   executable. Please see the specific documentation on these methods above for
   details on how to do this.
