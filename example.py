"""
Example of using the stdpopsim library with msprime.
"""
import msprime

# import daiquiri
# # Print out log messages for recombination map retrieval
# daiquiri.setup(level="DEBUG")

import stdpopsim
from stdpopsim import homo_sapiens

# print(homo_sapiens.genome)
# print(homo_sapiens.genome.chromosomes["chr22"])
chrom = homo_sapiens.genome.chromosomes["chr22"]
recomb_map = chrom.recombination_map()
# print(recomb_map)

# ts = msprime.simulate(10, recombination_map=recomb_map)
# print("ts = ", ts)
#     # recombination_map=h_sap.chr22.recombination_map(),
#     # length=h_sap.chr22.length,
#     # recombination_rate=h_sap.chr22.mean_recombination_rate,
#     # mutation_rate=h_sap.chr22.mean_mutation_rate,
#     # **model.asdict())





# import stdpopsim.pongo as pongo
# model = pongo.LockeEtAlPongoIM()
# model.debug()

# genetic_map = stdpopsim.get_genetic_map("h_sapiens", "HapmapII_GRCh37")
# recomb_map = genetic_map.get_chromosome_map("chr20")
# print(recomb_map)


model = homo_sapiens.GutenkunstThreePopOutOfAfrica()
model.debug()

# One sample each from YRI, CEU and CHB. There's no point in pushing
# the sampling strategy into the model generation
samples = [
    msprime.Sample(population=0, time=0),
    msprime.Sample(population=1, time=0),
    msprime.Sample(population=2, time=0)]

# recombination_map=h_sap.chr22.recombination_map()
# print(recombination_map)

ts = msprime.simulate(
    samples=samples,
    # recombination_map=chrom.recombination_map(),
    length=chrom.length,
    recombination_rate=chrom.mean_recombination_rate,
    mutation_rate=chrom.mean_mutation_rate,
    **model.asdict())

# # print(ts.tables)
# print("simulated:", ts.num_trees, ts.num_sites)

