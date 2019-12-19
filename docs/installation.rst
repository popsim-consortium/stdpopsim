.. _sec_installation:

============
Installation
============

There are two basic options for installing ``stdpopsim`` and its dependencies: either
using :ref:`sec_installation_conda` or :ref:`sec_installation_pip`.
We recommend conda for most users (particularly those using OSX or Windows),
although pip can be more convenient in certain cases.

.. _sec_installation_requirements:

************
Requirements
************

Library requirements for stdpopsim should be installed automatically when using
either conda or pip.

Stdpopsim requires Python 3.5 or later.


.. _sec_installation_conda:

*****
Conda
*****

.. todo:: The conda package for stdpopsim is currently under development and will
    be available shortly. Please use pip for now.


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

We use `msprime <https://msprime.readthedocs.io/>`_ as the
default simulation engine, which has some system level dependencies
and requires a functioning compiler. Please see the msprime
`installation documentation
<https://msprime.readthedocs.io/en/stable/installation.html>`_ for
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


