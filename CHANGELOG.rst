--------
Upcoming
--------

--------------------
[0.3.0] - 2025-03-18
--------------------

Major feature release adding documented, stable support for simulating selection (via SLiM),
as well as new species and models.

**Bug fixes**:

- Incorrect scaling of the DroMel/LognormalPlusPositive_R16 DFE led to negative selection
  coefficients that were too large (many on the order of -1e3); scaling was fixed for the
  log scale (:user:`clararehmann`, :pr:`1699`)

- Factor of two error fixed in HomSap/Gamma_K17 DFE (:user:`RyanGutenkunst`, :pr:`1478`)

- Discretization for exponential growth using the SLiM engine led to the final size
  in growing/declining populations not matching (by a small amount) what is listed
  in the model; improved discretization scheme (:user:`petrelharp`, :pr:`1622`)

- Updates to SLiM support: updated the `active` flag in the SLiM code to be integer.`

**Breaking changes**:

- The 2018 Browning et al HomSap demographic model previously named
  "AmericanAdmixture_4B11" is now named "AmericanAdmixture_4B18"; the old name
  works but throws a deprecation warning. (:user:`petrelharp`, :pr:`1603`)

- The `time_units` attribute of tree sequence metadata is now set to "generations"
  for SLiM output, avoiding the "time units mismatch warning".
  (:user:`nspope`, :pr:`1567`)

- Previously, Contigs specified with `left` and `right` arguments produced
  simulations having coordinates shifted so they were relative to `left`. This
  is no longer the case: coordinates in resulting simulations retain the
  original chromosome's coordinates. As a result, regions outside `[left,
  right)` will now have missing data. To produce simulation output with the
  previous behavior, use the `.trim()` method of the resulting TreeSequence.
  So, the Contig's `.original_coordinates` attribute is now called
  `.coordinates`, and several methods (e.g., `Contig.basic_contig()`) now take
  `left` and `right` arguments as well as `length`. Finally, arguments such as
  `relative_coordinates=True/False` to `Contig.dfe_breakpoints()` are no longer
  necessary and are deprecated. (:user:`petrelharp`, :pr:`1570`)

- To add the possibility for a relationship between dominance and selection
  coefficient, now each stdpopsim MutationType might have more than one
  underlying SLiM mutation type; so, where this is recorded in top-level
  metadata (under `ts.metadata["stdpopsim"]["DFEs"]`) is now a list
  instead of a single value. This will not affect anyone who is not
  parsing the metadata related to DFEs.

- SLiM extended events and selective sweep infrastructure have been
  moved from the `stdpopsim.ext` namespace into ``stdpopsim`` proper.

- The `length_multiplier` option to `Species.get_contig` is deprecated and
  prints a warning. The options `left` and `right` should be used to truncate a
  contig, instead. (:user:`nspope`, :pr:`1605`)

- Added `assembly_source` and `assembly_build_version` properties to
  `Genome` objects, so that different genomes could be based on
  different sources (previously, all were with reference to the Ensembl
  release 103). (:user:`andrewkern`, :pr:1646)

**New features**:

- Relationship between dominance and selection coefficient:
    Added the `dominance_coeff_list` argument to `MutationType`, allowing
    for DFEs with a discretized relationship between h and s.
    (:user:`petrelharp`, :pr:`1498`)

- Model specific recombination rates:
    Added the `recombination_rate` attribute to demographic model,
    allowing for model specific recombination rates in addition to
    the species default rates.
    (:user:`gregorgorjanc`, :pr:`1591`)

- Added support for the SLiM engine on Windows. (:user:`petrelharp`, :pr:`1571`)


**New species**:

- Mus musculus (:user:`peterdfields`, :pr:`1437`).
  QC'd by :user:`igronau`, :pr:`1454`.

- Orzya sativa (:user:`ornobalam`, :pr:`1453`).
  QC'd by :user:`minesrebollo`, :pr:`1461`

- Phocoena sinus (:user:`igronau`, :pr:`1514`).
  QC'd by :user:`ckyriazis`, :pr:`1538`

- Gorila gorila (:user:`ChristianHuber`, :pr:`1517`)
  QC'd by :user:`andrewkern`, :pr:`1701`

- Rattus norvegicus (:user:`kevinkorfman`, :pr:`1700`)
  QC'd by :user:`andrewkern`, :pr:`1706`

- Sus scrofa (:user:`aprilyuzhang`, :pr:`1672`)
  QC'd by :user:`gregorgorjanc`, :pr:`1689`


**New DFEs**:

- Generic "uniform" DFE (:user:`petrelharp`, :pr:`1492`)

- HomSap/LogNormal_H17 (:user:`RyanGutenkunst`, :pr:`1480`)
  QC'd by :user:`clararehmann`, :pr:`1698`

- HomSap/Mixed_K23 (:user:`chriscrsmith`, :pr:`1505`)
  QC'd by :user:`clararehmann`, :pr:`1696`

- PhoSin/Gamma_R22 (:user:`igronau`, :pr:`1547`)
  QC'd by :user:`ckyriazis`, :pr:`1560`

- AraTha/Gamma_H18 (:user:`chriscrsmith`, :pr:`1324`)
  QC'd by :user:`clararehmann`, :pr:`1647`

- HomSap/GammaPos_H17 and DroMel/GammaPos_Z21 (:user:`petrelharp`, :pr:`1656`)
  QC'd by :user:`clararehmann`, :pr:`1694` & :pr:`1695`


**New demographic models**:

- MusMus/DomesticusEurope_1F22 (:user:`peterdfields`, :pr:`1485`)
  QC'd by :user:`igronau`, :pr:`1531`

- MusMus/MusculusKorea_1F22 (:user:`peterdfields`, :pr:`1485`)
  QC'd by :user:`igronau`, :pr:`1531`

- MusMus/CastaneusIndia_1F22 (:user:`peterdfields`, :pr:`1485`)
  QC'd by :user:`igronau`, :pr:`1531`

- OrySat/BottleneckMigration_3C07 (:user:`ornobalam`, :pr:`1453`)
  QC'd by :user:`minesrebollo`, :pr:`1466` and :user:`petrelharp`, :pr:`1524`

- PhoSin/Vaquita2Epoch_1R22 (:user:`igronau`, :pr:`1526`)
  QC'd by :user:`ckyriazis`, :pr:`1538`

- CanFam/EarlyWolfAdmixture_6F14 (:user:`agladstein`, :pr:`1632`)
  QC'd by :user:`gregorgorjanc`, :pr:`1644`

- BosTau/HolsteinFriesian_1B16, BosTau/Angus_1B16, BosTau/Fleckvieh_1B16, and
  BosTau/Jersey_1B16 (:user:`gregorgorjanc`, :pr:`1593`)

**New annotations**:

- PhoSin exons and CDS (:user:`chriscrsmith`, :pr:`1520`)

**New QCs**:

- DroMel/Gamma_H17 (:user:`clararehmann`, :pr:`1668`)

- DroMel/LognormalPlusPositive_R16 (:user:`clararehmann`, :pr:`1699`)

- HomSap/Gamma_K17 (:user:`clararehmann`, :pr:`1687`).

- MusMus/Gamma_B21 (:user:`clararehmann`, :pr:`1669`)


---------------------
[0.2.1a] - 2024-07-06
---------------------

This was the alpha release for 0.3.0 (before releasing, we decided to change
the version number to 0.3.0).

--------------------
[0.2.0] - 2022-11-01
--------------------

Major feature release adding many new species and models, as well as support
for simulating selection via SLiM.

**Bug fixes**:

- Parameters in the HomSap/Zigzag_1S14 model now match those in Schiffels &
  Durbin (2014) (:user:`grahamgower`, :pr:`750`).

- Recombination rate for DroMel chr4 changed to 0
  (:user:`izabelcavassim`, :pr:`1092`).

- Per-chromosome mean recombination rates for HomSap were incorrectly
  calculated (:user:`nspope`, :pr:`1345`).

**Breaking changes**:

- Removed `GeneticMap` class from public API (:user:`jeromekelleher`, :pr:`713`).

- Samples are now specified via population/individual pairs, using
  species/chromosome ploidy.  The old API for specifying haploid samples via
  population index has been retained, but is deprecated and will be
  removed at some point (:user:`nspope`, :pr:`1361`).

**New species**:

- Aedes aegypti (:user:`manolofperez`, :pr:`871`).
  QC'd by :user:`petrelharp`, :pr:`893`.

- Anas platyrhynchos (:user:`petrelharp`, :pr:`826`).
  QC'd by :user:`igronau`, :pr:`1070`.

- Anolis carolinensis (:user:`vcaudill`, :pr:`874`).
  QC'd by :user:`andrewkern`, :pr:`896`.

- Anopheles gambiae (:user:`andrewkern`, :pr:`856`).
  QC'd by :user:`petrelharp`, :pr:`906`.

- Apis mellifera (:user:`janaobsteter`, :pr:`1025`).
  QC'd by :user:`manolofperez`, :pr:`1268`.

- Bos taurus (:user:`grahamgower`, :pr:`600`).
  QC'd by :user:`gtsambos`, :pr:`1269`.

- Caenorhabditis elegans (:user:`attrna`, :pr:`910`).
  QC'd by :user:`chriscrsmith`, :pr:`1265`.

- Chlamydomonas reinhardtii (:user:`aays`, :pr:`863`).
  QC'd by :user:`izabelcavassim`, :pr:`1067`.

- Drosophila sechellia (:user:`jradrion`, :pr:`872`).
  QC'd by :user:`vitorpavinato`, :pr:`1264`.

- Gasterosteus aculeatus (:user:`vitorpavinato`, :pr:`1105`).
  QC'd by :user:`manolofperez`, :pr:`1253`.

- Helianthus annuus (:user:`chriscrsmith`, :pr:`1218`).
  QC'd by :user:`xin-huang`, :pr:`1250`.

- Heliconius melpomene (:user:`percyfal`, :pr:`870`).
  QC'd by :user:`noscode`, :pr:`1165`.

- Pan troglodytes (:user:`xin-huang`, :pr:`1215`).
  QC'd by :user:`janaobsteter`, :pr:`1291`.

- Papio anubis (:user:`saurabhbelsare`, :pr:`1216`).
  QC'd by :user:`mufernando`, :pr:`1263`.

- Streptococcus agalactiae (:user:`jeanrjc`, :pr:`854`).
  QC'd by :user:`vitorpavinato`, :pr:`1251`.

**New models**:

- AnaPla/MallardBlackDuck_2L19 (:user:`petrelharp`, :pr:`883`).
  QC'd by :user:`igronau`, :pr:`1021`.

- AnoGam/GabonAg1000G_1A17 (:user:`andrewkern`, :pr:`856`).
  QC'd by :user:`petrelharp`, :pr:`1279`.

- BosTau/HolsteinFriesian_1M13 (:user:`grahamgower` and :user:`iechevarriaz`, :pr:`558` and :pr:`600`).
  QC'd by :user:`igronau`, :pr:`1272`.

- HomSap/OutOfAfricaExtendedNeandertalAdmixturePulse_3I21
  (:user:`leonardolasi`, :pr:`1066`).
  QC'd by :user:`awohns`, :pr:`1259`.

- HomSap/OutOfAfrica_4J17 (:user:`rwaples`, :pr:`726`).
  QC'd by :user:`jeffspence`, :pr:`1246`.

- HomSap/Africa_1B08 (:user:`izabelcavassim`, :pr:`993`).
  QC'd by :user:`petrelharp`, :pr:`995`.

- HomSap/AncientEurope_4A21 (:user:`alipearson`, :pr:`941`).
  QC'd by :user:`mufernando`, :pr:`1256`.

- PanTro/BonoboGhost_4K19 (:user:`xin-huang`, :pr:`1215`).
  QC'd by :user:`kuhlwilm`, :pr:`1370`.

- PapAnu/SinglePopSMCpp_1W22 (:user:`saurabhbelsare`, :pr:`1216`).
  QC'd by :user:`attrna`, :pr:`1261`.

**New genetic maps**:

- CaeEle/RockmanRIAIL_ce11 (:user:`attrna`, :pr:`910`).

- DroMel/ComeronCrossoverV2_dm6 liftover (:user:`grahamgower`, :pr:`592`).

- HomSap/HapMapII_GRCh38 liftover (:user:`saurabhbelsare`, :pr:`1301`).

- HomSap/DeCodeSexAveraged_GRCh38 liftover (:user:`saurabhbelsare`, :pr:`1301`).

- HomSap/PyrhoXXX_GRCh38 (:user:`jeffspence`, :pr:`572` and :pr:`575`),
  for XXX in ACB, ASW, BEB, CDX, CEU, CHB, CHS, CLM, ESN, FIN, GBR, GIH, GWD,
  IBS, ITU, JPT, KHV, LWK, MSL, MXL, PEL, PJL, PUR, STU, TSI, and YRI.

- PapAnu/Pyrho_PAnubis1_0 (:user:`saurabhbelsare`, :pr:`1216`)

**New features**:

- Distributions of fitness effects ("DFEs") defined over genomic intervals
  (:user:`mufernando`, :pr:`644`; :user:`izabelcavassim`, :pr:`1002`;
  plus numerous others).

- DFE simulation via SLiM
  (:user:`mufernando`, :pr:`930`; plus numerous others).

- Metadata for tree sequences produced by SLiM
  (:user:`mufernando`, :pr:`1152`).

- Per-generation fitness statistics for SLiM simulations
  (:user:`petrelharp`, :pr:`1200`).

- Selective sweep simulation and allele frequency conditioning via SLiM
  (:user:`grahamgower`, :pr:`462`; :user:`nspope`, :pr:`1341`).

- Gene conversion simulation via msprime and SLiM
  (:user:`fbaumdicker`, :pr:`1106`; :user:`petrelharp`, :pr:`1355`).

- Genome annotation tracks
  (:user:`andrewkern`, :pr:`560` and :pr:`960`).

- Masking intervals in simulated data
  (:user:`apragsdale`, :pr:`664`).

- Method to get generic contig of arbitrary length for a species
  (:user:`apragsdale`, :pr:`664`).

- Method to get contig from a segment of a named chromosome
  (:user:`nspope`, :pr:`1348`).

- Pass keyworded arguments from simulation engine to msprime
  (:user:`awohns`, :pr:`736`).

- Use msprime 1.0 for simulation from msprime engine
  (:user:`jeromekelleher`, :pr:`764`).

- Use SLiM 4.0 for simulation from SLiM engine
  (:user:`petrelharp`, :pr:`1326`).

- Mutation rates can be stored in catalog models
  (:user:`apragsdale`, :pr:`839`).

- Ploidy is a species and chromosome attribute
  (:user:`nspope`, :pr:`1361`).

- Mutations from SLiM simulations converted to nucleotides
  (:user:`nspope`, :pr:`1356`).

- Various improvements and fixes to the documentation and error messaging.

**Additions to CLI**:

- Sample specification has switched from positional and haploid (e.g.
  ``stdpopsim HomSap -d OutOfAfrica_3G09 6 0 10``) to named with species-specific
  ploidy (equivalent to ``stdpopsim HomSap -d OutOfAfrica_3G09 YRI:3 CEU:0
  CHB:5``). Positional sample specification is still supported but will raise a
  deprecation warning.

- Arguments ``--dfe``, ``--dfe-interval``, ``--dfe-bed-file``, ``--help-dfe``
  for specifying DFEs (:user:`izabelcavassim`, :pr:`1052`).

- Arguments ``--help-annotations``, ``--dfe-annotation`` for associating annotation
  tracks with DFEs (:user:`andrewkern`, :pr:`1117`).

- Argument ``--length`` for simulating from a generic contig
  (:user:`apragsdale`, :pr:`664`).

- Arguments ``--inclusion-mask``, ``--exclusion-mask`` for masking simulated sequences
  (:user:`apragsdale`, :pr:`664`).

- Arguments ``--left`` and ``--right`` for simulating an interval on a named chromosome
  (:user:`nspope`, :pr:`1348`)

- Argument ``--keep-mutation-ids-as-alleles`` retains SLiM mutation IDs for
  allele codes instead of converting these to nucleotides (:user:`nspope`, :pr:`1356`).

**Catalog maintenance infrastructure**:

- Quality control infrastructure for DFEs
  (:user:`xin-huang`, :pr:`1292`).

- Pull species information from NCBI
  (:user:`andrewkern`, :pr:`875`).

- Automated species addition to catalog
  (:user:`jeromekelleher`, :pr:`790`).

- Github issue template for requesting addition of species
  (:user:`petrelharp`, :pr:`772`).

- Tools for assembly liftover
  (:user:`jradrion`, :pr:`574`).

- Pull genome data from Ensembl
  (:user:`jeromekelleher`, :pr:`563`).

**New annotations**:

- AraTha/araport_11 (:user:`andrewkern`, :pr:`1327`).

- DroMel/FlyBase_BDGP6.32.51 (:user:`andrewkern`, :pr:`1042`).

- HomSap/ensembl_havana_104 (:user:`andrewkern`, :pr:`960`).

**New DFEs**:

- DroMel/Gamma_H17 (:user:`izabelcavassim`, :pr:`1046`).

- DroMel/LognormalPlusPositive_R16 (:user:`apragsdale`, :pr:`1178`).

- HomSap/Gamma_K17 (:user:`izabelcavassim`, :pr:`1002`).

- HomSap/Gamma_H17 (:user:`chriscrsmith`, :pr:`1099`).

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
