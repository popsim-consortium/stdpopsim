"""
Example of using the stdpopsim library with msprime.
"""
# import daiquiri
# # Print out log messages for recombination map retrieval
# daiquiri.setup(level="DEBUG")

import stdpopsim

species = stdpopsim.get_species("homo_sapiens")
contig = species.get_contig("chr22", length_multiplier=0.1)

model = species.get_model("constant")
samples = model.get_samples(2)
ts = model.run(contig, samples)

model = species.get_model("2_epoch")
samples = model.get_samples(2)
ts = model.run(contig, samples)

model = species.get_model("ooa", 3)
# 2 sampes from pop 0, 1 from pop 2
samples = model.get_samples(2, 0, 1)
ts = model.run(contig, samples)

