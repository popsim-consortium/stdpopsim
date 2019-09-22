"""
Example of using the stdpopsim library with msprime.
"""
import msprime
import stdpopsim
from stdpopsim import homo_sapiens

chrom = homo_sapiens.genome.chromosomes["chr22"]
recomb_map = chrom.recombination_map()


model = homo_sapiens.GenericConstantSize()
print("Human constant size debug")
model.debug()

# three samples
samples = [
    msprime.Sample(population=0, time=0),
    msprime.Sample(population=0, time=0),
    msprime.Sample(population=0, time=0)]

ts = msprime.simulate(
    samples=samples,
    recombination_map=chrom.recombination_map(),
    mutation_rate=chrom.default_mutation_rate,
    **model.asdict())
print("human generic simulated:", ts.num_trees, ts.num_sites)
# GenericTwoEpoch
model = stdpopsim.homo_sapiens.GenericTwoEpoch()
print("Human two epoch debug")
model.debug()

# two samples
samples = [
    msprime.Sample(population=0, time=0),
    msprime.Sample(population=0, time=0)]

ts = msprime.simulate(
    samples=samples,
    recombination_map=chrom.recombination_map(),
    mutation_rate=chrom.default_mutation_rate,
    **model.asdict())
print("Human GenericTwoEpoch simulated:", ts.num_trees, ts.num_sites)

# do generic thing with arabidopsis
chrom = stdpopsim.arabidopsis_thaliana.genome.chromosomes["chr5"]
recomb_map = chrom.recombination_map()

model = stdpopsim.arabidopsis_thaliana.GenericConstantSize()
print("Arabidopsis constant size debug")
model.debug()

# two samples
samples = [
    msprime.Sample(population=0, time=0),
    msprime.Sample(population=0, time=0)]

ts = msprime.simulate(
    samples=samples,
    recombination_map=chrom.recombination_map(),
    mutation_rate=chrom.default_mutation_rate,
    **model.asdict())
print("Arabidopsis generic simulated:", ts.num_trees, ts.num_sites)

