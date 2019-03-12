"""
Example of using the stdpopsim library with msprime.
"""
import msprime


import stdpopsim
from stdpopsim import e_coli

chrom = e_coli.genome.chromosomes["chr"]
#recomb_map = chrom.recombination_map()

model = e_coli.LapierreConstant()
model.debug()

samples=[msprime.Sample(population=0, time=0)]*100

ts = msprime.simulate(
    samples=samples,
    mutation_rate=chrom.mean_mutation_rate,
    **model.asdict())
print("simulated:", ts.num_trees, ts.num_sites)
print("mean mut rate", chrom.mean_mutation_rate)
