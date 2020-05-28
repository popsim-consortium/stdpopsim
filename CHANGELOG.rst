--------------------
[0.1.2] - 2020-05-29
--------------------

Important bugfix and feature release, recommended for all users.

Significant errors in the HomSap/OutOfAfrica_3G09 and HomSap/OutOfAfrica_2T12
models have been fixed. **We recommend repeating any analyses performed using
these models**. See `here
<https://github.com/jeromekelleher/msprime-model-errors>`__ for more details on
the error in the three population Out of Africa model and analysis of the
differences from the correct model.

The recombination rate for AraTha was also off by a factor of 10.
**We recommend repeating any analyses performed using this species**.

**Bug fixes**:

- Fix error in HomSap/OutOfAfrica_3G09 model, in which migration between
  ancestral African and European populations was allowed to continue in the
  most ancient time period (:user:`apragsdale`, :pr:`496`, :issue:`516`).

- Fix similar error in HomSap/OutOfAfrica_2T12 model
  (:user:`ndukler`, :pr:`520`, :issue:`516`).

- Fix recombination rate estimate for AraTha (:user:`grahamgower`,
  :issue:`537`, :pr:`527`), which was off by a factor of 10.

- Require attrs >=19.10 (:user:`grahamgower`, :pr:`399`, :issue:`394`)

**New species**:

- Canis familiaris (:user:`grahamgower`, :pr:`375`).

- Pongo abelii (:user:`apragsdale`, :pr:`363`).

**New models**:

- HomSap/PapuansOutOfAfrica_10J19 model (:user:`grahamgower`, :pr:`372`).
  QC'd by :user:`noscode`, :pr:`387`.

- HomSap/AshkSub_7G19 model (:user:`agladstein`, :pr:`494`).
  QC'd by :user:`ndukler`, :pr:`536`.

**New features**:

- SLiM simulation engine (:user:`grahamgower`, :pr:`409`, plus numerous others.
  See e.g. :issue:`132` and :issue:`133` for background.)

- Support for DTWF, SMC, and SMC' models in msprime engine
  (:user:`grahamgower`, :pr:`398`, :issue:`392`).

- Warnings for users running simulations on non-autosomes
  (:user:`grahamgower`, :pr:`407`).

- Migrate all genetic map data to AWS (:user:`ndukler`, :pr:`514`, :issue:`335`)

- Warnings for users running simulations on non QC'd models
  (:user:`grahamgower`, :pr:`525`).

- Add `generation_time` (default=1) attribute to generic models
  (:user:`grahamgower`, :pr:`477`, :issue:`471`).

- Various documentation and citation improvements.

**Breaking changes**:

- Move the --quiet/-q command line option to the top-level. Previously
  we would write ``stdpopsim HomSap -q 10`` whereas we now write
  ``stdpopsim -q HomSap``. (:user:`jeromekelleher`, :issue:`515`, :pr:`547`)

- The long form ``--verbosity`` argument has been changed to ``--verbose``
  (:pr:`547`).

- Removed DroMel chrM (:user:`grahamgower`, :pr:`528`, :issue:`405`).

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
