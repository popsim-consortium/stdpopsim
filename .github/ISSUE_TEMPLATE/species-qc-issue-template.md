---
name: Species QC issue template
about: Quality control process for species addition
title: QC for {species}
labels: Species QC
assignees: ''

---

**PR for new species:** {link to Pull Request}

If you volunteer to QC this species, please use the checklist below.
While this list is intended to be comprehensive, it may not be exhaustive.
Where relevant, the QC reviewer should identify that parameter values match
those given in the linked citation(s).

- [ ] Recombination rate.
  - This might be genome-wide, or per-chromosome. Both are fine.
  - Check there's a comment describing where it came from, and/or how calculated.
  - From a publication? Check the value(s) match the publication.
  - Calculated somehow? Average over a recombination map? Redo the calculation.
- [ ] Mutation rate.
  - This might be genome-wide, or per-chromosome. Both are fine.
  - Check there's a comment describing where it came from, and/or how calculated.
  - From a publication? Check the value(s) match the publication.
  - Calculated somehow? Redo the calculation.
- [ ] Recombination map (if present).
  - Does it match the assembly? Liftover is fine, if clearly stated.
  - Is the description/long_description a good summary of how the map was created?
- [ ] Population size.
- [ ] Generation time.

For each citation, check:
- Doi link.
- Is publication a preprint? Is there a peer-reviewed publication instead?
- Is the year correct.
- Is the author correct (spelling, accents/ligatures/etc.)

Citations are required for:
- [ ] Genome reference assembly.
- [ ] Mutation rate.
- [ ] Recombination rate.
- [ ] Recombination map(s) (if relevant).
- [ ] Population size.
- [ ] Generation time.
