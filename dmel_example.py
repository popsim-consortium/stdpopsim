"""
Example of using the stdpopsim library with msprime.
"""
import msprime


import stdpopsim
from stdpopsim import drosophila_melanogaster

chr_str = "chr2L"
chrom = drosophila_melanogaster.genome.chromosomes[chr_str]
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
    mutation_rate=chrom.default_mutation_rate,
    **model.asdict())
print("simulated:", ts.num_trees, ts.num_sites)

print("====================\n")
print(f"Testing Li and Stephan (2006) model- with on {chr_str} map")
model = drosophila_melanogaster.LiStephanTwoPopulation()
model.debug()

samples = [msprime.Sample(population=0, time=0)] * 10 + [msprime.Sample(population=1, time=0)] * 10


ts = msprime.simulate(
    samples=samples,
    recombination_map=chrom.recombination_map(),
    mutation_rate=chrom.default_mutation_rate,
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
    mutation_rate=chrom.default_mutation_rate,
    **model.asdict())
print("simulated:", ts.num_trees, ts.num_sites)

print("====================\n")
print(f"Testing Sheehan and Song (2016) model on {chr_str} map.")
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
    mutation_rate=chrom.default_mutation_rate,
    **model.asdict())
print("simulated:", ts.num_trees, ts.num_sites)

