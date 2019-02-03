"""
Example of using the stdpopsim library with msprime.
"""
import msprime

# import daiquiri
# # Print out log messages for recombination map retrieval
# daiquiri.setup(level="DEBUG")

import stdpopsim
from stdpopsim import homo_sapiens

chrom = homo_sapiens.genome.chromosomes["chr22"]
recomb_map = chrom.recombination_map()

model = homo_sapiens.GutenkunstThreePopOutOfAfrica()
model.debug()

# One sample each from YRI, CEU and CHB. There's no point in pushing
# the sampling strategy into the model generation
samples = [
    msprime.Sample(population=0, time=0),
    msprime.Sample(population=1, time=0),
    msprime.Sample(population=2, time=0)]

ts = msprime.simulate(
    samples=samples,
    recombination_map=chrom.recombination_map(),
    mutation_rate=chrom.mean_mutation_rate,
    **model.asdict())
print("simulated:", ts.num_trees, ts.num_sites)
