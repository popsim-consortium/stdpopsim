"""
Example of using the stdpopsim library with msprime.
"""
import msprime


import stdpopsim
from stdpopsim import drosophila_melanogaster

chrom = drosophila_melanogaster.genome.chromosomes["chr2R"]
recomb_map = chrom.recombination_map()

print("Testing Li and Stephan (2006) model")
model = drosophila_melanogaster.LiStephanTwoPopulation()
model.debug()

samples = [msprime.Sample(population=0, time=0),
           msprime.Sample(population=0, time=0),
           msprime.Sample(population=1, time=0),
           msprime.Sample(population=1, time=0)]


ts = msprime.simulate(
    samples=samples,
    # recombination rate placeholder for quick runtime
    recombination_rate=1e-09,
    mutation_rate=chrom.mean_mutation_rate,
    **model.asdict())
print("simulated:", ts.num_trees, ts.num_sites)

print("====================\n")
print("Testing Sheehan and Song (2016) model")
model = drosophila_melanogaster.SheehanSongThreeEpoch()

model.debug()

samples = [msprime.Sample(population=0, time=0),
           msprime.Sample(population=0, time=0)
          ]

ts = msprime.simulate(
    samples=samples,
    # recombination rate placeholder for quick runtime
    recombination_rate=1e-09,
    mutation_rate=chrom.mean_mutation_rate,
    **model.asdict())
print("simulated:", ts.num_trees, ts.num_sites)

print("====================\n")
print("Testing Sheehan and Song (2016) model on chrom2R map.")
print("This will take a while- perhaps 2 hours\n")
model = drosophila_melanogaster.SheehanSongThreeEpoch()

model.debug()

samples = [msprime.Sample(population=0, time=0),
           msprime.Sample(population=0, time=0)
          ]

ts = msprime.simulate(
    samples=samples,
    # actual recombination map; slow
    recombination_map=chrom.recombination_map(),
    mutation_rate=chrom.mean_mutation_rate,
    **model.asdict())
print("simulated:", ts.num_trees, ts.num_sites)

