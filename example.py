"""
Example of using the stdpopsim library with msprime.
"""
import msprime
import stdpopsim.h_sapiens as h_sap

# import stdpopsim.pongo as pongo
# model = pongo.LockeEtAlPongoIM()
# model.debug()

model = h_sap.GutenkunstThreePopOutOfAfrica()
model.debug()

# One sample each from YRI, CEU and CHB. There's no point in pushing
# the sampling strategy into the model generation
samples = [
    msprime.Sample(population=0, time=0),
    msprime.Sample(population=1, time=0),
    msprime.Sample(population=2, time=0)]

ts = msprime.simulate(
    samples=samples,
    recombination_map=h_sap.chr22.recombination_map(),
    # length=h_sap.chr22.length,
    # recombination_rate=h_sap.chr22.mean_recombination_rate,
    mutation_rate=h_sap.chr22.mean_mutation_rate,
    **model.asdict())

# print(ts.tables)
print("simulated:", ts.num_trees, ts.num_sites)

