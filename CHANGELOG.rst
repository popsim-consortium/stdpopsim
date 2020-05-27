--------------------
[0.1.2] - 2020-XX-XX
--------------------


**Bug fixes**:


**New features**:

**Breaking changes**:

- Move the --quiet/-q command line option to the top-level. Previously
  we would write ``stdpopsim HomSap -q 10`` whereas we now write
  ``stdpopsim -q HomSap``. (:user:`jeromekelleher`, :issue:`515`, :pr:`547`)

- The long form ``--verbosity`` argument has been changed to ``--verbose``
  (:pr:`547`).

--------------------
[0.1.1] - 2020-01-02
--------------------

Bugfix release. Fixes some distribution issues and temporarily removes the
PonPyg species.

**Bug fixes**:

- Pin the msprime and attrs packages to resolve some distribution problems
  (:issue:`366`; :user:`jgallowa07` and :user:`gtsambos`).

**New features**:

- Provide citations for the genome assembly (:issue:`359`, :pr:`360`;
  :user:`andrewkern` and :user:`grahamgower`).

**Breaking changes**:

- Temporarily remove the PonPyg species from the catalog to provide time
  to fix issues with genomes and multi-species models (:issue:`365`).

--------------------
[0.1.0] - 2019-12-18
--------------------

Initial release.
