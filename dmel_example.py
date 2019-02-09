"""
Example of using the stdpopsim library with msprime.
"""
import msprime


import stdpopsim
from stdpopsim import drosophila_melanogaster

chrom = drosophila_melanogaster.genome.chromosomes["chrX"]
recomb_map = chrom.recombination_map()

model = drosophila_melanogaster.SheehanSongThreeEpoch()
model.debug()

samples = [msprime.Sample(population=0, time=0),msprime.Sample(population=0,
                                                               time=0)]

ts = msprime.simulate(
    samples=samples,
    recombination_map=chrom.recombination_map(),
    mutation_rate=chrom.mean_mutation_rate,
    **model.asdict())
print("simulated:", ts.num_trees, ts.num_sites)
