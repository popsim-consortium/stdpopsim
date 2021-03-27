---
name: New species
about: Adding a new species to the stdpopsim catalogue.
title:
labels: enhancement
assignees: ''

---

- **Name:** scientific name of the species
- **Common name:**

**Description:** any other discussion or useful background? (optional)

*Note: Maybe you don't know a lot of the things below - but, that's OK!
The criteria are that (a) the choice is reasonable, and (b) the rationale is clear.
What we want is a "reasonable estimate" (or, "best guess") - if you were to use a number in a paper,
what would you use and how would you justify it to your collaborators?
For instance, for "generation time", best would be to cite a paper that actually measures
or estimates generation time. But most species don't have this;
next best would be to give a number used in the literature, and provide a citation that used it
(and hopefully the publication gives a justification also).
Lacking this, you might provide a number from a related species (and a citation for this).*

*Here is the checklist of things that we need to add a species to the catalogue.
Each thing should be provided, with a justification (maybe short) and a citation.*

**Demographic information:**

- [] generation time in years
- [] "default" population size


**Chromosome structure:**

*Note:* The assembly should be *chromosome-level*, i.e., not composed of thousands of scaffolds.

- [] ensembl assembly ID *or*, if the assembly is not in Ensembl, a list of chromosomes with *name* and *length* (in bp)

*Note: species in Ensembl can be found in one of these lists:
[vertebrates](https://uswest.ensembl.org/info/about/species.html),
["metazoa"](https://metazoa.ensembl.org/species.html),
[plants](https://plants.ensembl.org/species.html),
[fungi](https://fungi.ensembl.org/species.html),
[protists](http://protists.ensembl.org/species.html),
or [bacteria](http://bacteria.ensembl.org/species.html).*


**Recombination rates:**

- [] genetic map of recombination rates (as a hapmap or csv file)
- [] genome-wide mean recombination rate (which will be used as the default)


**Mutation rate:**

- [] genome-wide mean mutation rate


**Demographic model:**

- [] as a list of population sizes, growth rates, migration rates, etcetera. **(optional)**
- [] citation


**Other information:**

These are things we don't currently use, but will want to use in the future:

- [] sex determination system, and which chromosomes are the sex chromosome(s)
- [] ploidy
