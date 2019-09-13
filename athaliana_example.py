"""
A. thaliana msprime example
"""
import msprime

from stdpopsim import arabidopsis_thaliana

chrom = arabidopsis_thaliana.genome.chromosomes["chr1"]
recomb_map = chrom.recombination_map()

model = arabidopsis_thaliana.Durvasula2017MSMC()
model.debug()

samples = [
    msprime.Sample(population=0, time=0), msprime.Sample(population=0, time=0)]

ts = msprime.simulate(
    samples=samples,
    recombination_map=recomb_map,
    mutation_rate=chrom.default_mutation_rate,
    **model.asdict())
print("simulated:", ts.num_trees, ts.num_sites)
