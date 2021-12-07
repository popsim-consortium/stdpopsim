import stdpopsim

from . import genome_data

# ref for mutation, assembly, pop size
_DaCunha_et_al = stdpopsim.Citation(
    author="Da Cunha et al",
    year=2014,
    doi="https://doi.org/10.1038/ncomms5544",
)

# Per site per generation
# Estimated in Da Cunha et al. as 0.56 SNP/Mb/Year
# See mutation rate for CC17 at Table 2 (SNP/Mb/Year).
# Given that the generation time is estimated at 1 generation per day we have :
# mu =  0.56/(1e6*365)
_mutation_rate = {"1": 1.53e-9}

# no cross-over recombination in bacteria
_recombination_rate = {"1": 0}

# gene conversion rate (HGT rate):
# hard to estimate precisely. What is often estimated is the r/m or rho/theta ratio.
# in Oliveira et. al, PNAS, 2016, they estimate those rates on coregenome
# (conservative estimates) it was found for the species (not CC17 in particular) that
# r/m = 7.42% and rho/theta = 0.964%
# In Almeida et al, 2017, mSystems, they identified for a population of CC17,
# 502 homoplasic SNPs over 16572 SNPs, which gives a r/m ratio of 3%
# In Lefebure et al, Genome Biology 2007, they found that between 3% and 18% of the genes
#  were recombinant (depending on the methods used). This illustrates also that this
#  species do not recombine much
# In Da Cunha et al, 2016, Nat Comms, CC17 is shown to have a cumulated size of
# recombination fragment that is about 10% that of other recombining region of other
# clonal complexes of the same species, so CC17 does not recombine much in a species
# that is already not recombining a lot.
# Overall, as most of the estimation might be conservatives one can use a relative
# rate rho/theta = 10%.
# _gene_conversion_rate = {"1": _mutation_rate/10}

# Mean gene conversion  tract length
# In Brochet et al., 2008, PNAS, they estimate the size of transferred DNA experimentally
#  ([33, 334, 252, 28, 21, 105, 99, 14, 286, 49] kb), which leads to an estimation of a
# mean tract length of 122kb, which is plausible given the size that was identified in
# Almeida et al. 2017, mSystems. (25, 634, 778 and 880kb)
# This estimate (as well as others) was used in Cury et al, 2021, PCI Evol Biol, and
# discussed in more detail in response to a reviewer, that can be seen there :
# https://evolbiol.peercommunityin.org/articles/rec?id=356
#      > Revision round #3 > Author's reply (pdf).
# _mean_gene_conversion_tract_length = 120000


_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
    citations=[
        _DaCunha_et_al.because(stdpopsim.CiteReason.ASSEMBLY),
        _DaCunha_et_al.because(stdpopsim.CiteReason.MUT_RATE),
    ],
)
stdpopsim.utils.append_common_synonyms(_genome)

# Generation time :
# No estimation was done for S agalactiae, but we can use estimates from E coli (which
# grows a bit faster in the lab than S agalactiae), which ranges from around
# 0.5 to 20 generations per day.
# ref for lower range estimates : Michael A. Savageau, « Escherichia coli habitats,
# cell types, and molecular mechanisms of gene control », The american naturalist 122,
#  nᵒ 6 (1983): 732–744.
# ref for higher range estimates : Mohamed Ghalayini et al., « Evolution of a Dominant
#  Natural Isolate of Escherichia coli in the Human Gut over the Course of a Year
# Suggests a Neutral Evolution with Reduced Effective Population Size »,
# Applied and Environmental Microbiology 84, nᵒ 6 (5 janvier 2018): e02377-17,
#  https://doi.org/10.1128/AEM.02377-17.
# This is likely over estimated as it is only in the gut,
# where there is plenty of nutrient for bacteria to grow,
# unlike outside of it where growth might be slower.
# So we use the value of 1 generation per day.

# Population size
# We estimate it from the Watterson estimator :
# theta = 2.Ne.mu = S / sum_{i=1}^{k=n-1}(1/k)
# With n the number of samples and S the number of segregating sites.
# From Da Cunha et al, we have 3,922 polymorphic SNPs in 1,860Kbp sites.
# which means that they found that S = 3922 / 1.86Mb = 2.1×10−3 SNP/bp.
# Given k = 79 (number of isolates for CC17 in Table ), theta =  4.2e-4$.
# For mu=1.53×10−9 SNP/bp/generation and theta, Ne ~ 140000.


_species = stdpopsim.Species(
    id="StrAga",
    ensembl_id="NA",
    name="Streptococcus agalactiae",
    common_name="Group B Streptococcus",
    genome=_genome,
    generation_time=1 / 365,  # year / generations
    population_size=140000,
    citations=[
        _DaCunha_et_al.because(stdpopsim.CiteReason.POP_SIZE),
        stdpopsim.Citation(
            author="Savageau M.A.",
            year=1983,
            doi="https://doi.org/10.1086/284168",
            reasons={stdpopsim.CiteReason.GEN_TIME},
        ),
    ],
)

stdpopsim.register_species(_species)
