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

The QC reviewer should start a pull request that fills out the test stubs
with independently obtained values. The reviewer may look at the python code
for rationale provided in comments, but should ignore the actual code
as much as possible - comments in the code should give enough information
that it's obvious how to get the correct value from the provided references.
(In particular, we shouldn't copy-paste the value from the code into the test!)

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

The final PR should:
- [ ] fill out the test stubs
- [ ] delete the `pytest.mark.skip` lines that make the tests not run
- [ ] make sure they pass, talking to the original author to figure out discrepancies
