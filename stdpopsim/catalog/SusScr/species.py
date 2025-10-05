import stdpopsim

from . import genome_data

# Revisiting the Evolutionary History of Pigs via De Novo Mutation Rate
# Estimation in a Three-Generation Pedigree
_ZhangEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1016/j.gpb.2022.02.001",
    year=2022,
    author="Zhang et al.",
    reasons={
        # stdpopsim.CiteReason.MUT_RATE,  # on page 1042
        # (based on a pedigree with 9 individuals)
        # replaced by the more recent citation _RochusEtAl
        stdpopsim.CiteReason.POP_SIZE,  # on page 1048
    },
)

# Influence of harvesting pressure on demographic tactics: implications
# for wildlife management
_ServantyEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1111/j.1365-2664.2011.02017.x",
    year=2011,
    author="Servanty et al.",
    reasons={
        stdpopsim.CiteReason.GEN_TIME,  # on page 835
    },
)

# Genetic variation in recombination rate in the pig
_JohnssonEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1186/s12711-021-00643-0",
    year=2021,
    author="Johnsson et al.",
    reasons={stdpopsim.CiteReason.REC_RATE},
    # in Additional file 4: Table S3 (sex-averaged)
    # Additional file 2: Table S1 holds   male-specific values
    # Additional file 3: Table S2 holds female-specific values
    # (based on multigenerational pedigrees from 9 populations,
    # total 145,763 informative individuals, and
    # up to 80K SNP array genotypes)
)

# An improved pig reference genome sequence to enable pig genetics and
# genomics research
_WarrEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1093/gigascience/giaa051",
    year=2020,
    author="Warr et al.",
    reasons={stdpopsim.CiteReason.ASSEMBLY},
)

# Estimating mutation rate and characterising single nucleotide de novo
# mutations in pigs
_RochusEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1186/s12711-025-00967-1",
    year=2025,
    author="Rochus et al.",
    reasons={stdpopsim.CiteReason.MUT_RATE},  # Table 1, page 4, also page 5
    # (based on 46 trios from two populations)
)

# Ploidy - pig is diploid with XY sex chromosome system and has MT
_ploidy = 2
_ploidy_contig = {str(i): _ploidy for i in range(1, 19)}
_ploidy_contig["X"] = _ploidy
_ploidy_contig["Y"] = 1
_ploidy_contig["MT"] = 1

# De-novo mutation rate per site per generation
# from Rochus et al. (2025) Table 1, page 4, also page 5
_mutation_rate = 6.3e-9
_mutation_rate_contig = {str(i): _mutation_rate for i in range(1, 19)}
_mutation_rate_contig["X"] = _mutation_rate
_mutation_rate_contig["Y"] = _mutation_rate
_mutation_rate_contig["MT"] = _mutation_rate

# Recombination rate
# This mean rate for each chromosome was calculated from
# Johnsson et al. (2021) in Additional file 4: Table S3
# (summed rates along chromosome and divided by the reported physical length)
_tmp = [
    5.335247608983185e-09,
    8.6057260999254e-09,
    9.751141303348547e-09,
    9.79460049114614e-09,
    1.1992632296111e-08,
    8.786173521918499e-09,
    1.0887772285689678e-08,
    8.923031414594611e-09,
    9.337506174170301e-09,
    1.654233205589542e-08,
    1.2145120574910172e-08,
    1.6905062453467094e-08,
    6.162796630328349e-09,
    8.721569375548237e-09,
    8.114468762954893e-09,
    1.1407464686635653e-08,
    1.3779949516098884e-08,
    1.3679923888648167e-08,
]
_recombination_rate_contig = {str(i): _tmp[i - 1] for i in range(1, 19)}
# Setting X-chromosome recombination rate as an average across autosomes.
_total_length = sum(
    genome_data.data["chromosomes"][str(i)]["length"] for i in range(1, 19)
)
_weighted_recombination = sum(
    _tmp[i - 1] * genome_data.data["chromosomes"][str(i)]["length"]
    for i in range(1, 19)
)
_recombination_rate_autosome_mean = _weighted_recombination / _total_length
_recombination_rate_contig["X"] = _recombination_rate_autosome_mean
_recombination_rate_contig["Y"] = 0.0
_recombination_rate_contig["MT"] = 0.0

_genome = stdpopsim.Genome.from_data(
    genome_data=genome_data.data,
    ploidy=_ploidy_contig,
    mutation_rate=_mutation_rate_contig,
    recombination_rate=_recombination_rate_contig,
    citations=[
        _WarrEtAl,  # ASSEMBLY
        _JohnssonEtAl,  # REC_RATE
        _RochusEtAl,  # MUT_RATE
    ],
)

_species = stdpopsim.Species(
    id="SusScr",
    ensembl_id="sus_scrofa",
    name="Sus scrofa",
    common_name="Pig",
    separate_sexes=True,
    genome=_genome,
    ploidy=_ploidy,
    # Servanty et al. (2011) on page 837 write:
    # "We Calculated generation time as the inverse of the relative elasticity of
    # the population growth rate to a change in all recruitment parameters."
    # For comparison, Groenen et al. (2012, 10.1038/nature11622) used a best guess
    # of 5 years (supplement page 56) and this value is cited very often
    # (e.g., Wang et al., 2025, 10.1016/j.xgen.2025.100954).
    # To improve upon this best guess, Zhang et al. (2022) assumed a generation time
    # of 3 years, as the age of parents at the first litter, but pigs have
    # multiple parities in their lifetime, so 3 years must be an underestimate.
    # Servanty et al. (2011) seems to be the most reliable source with some
    # concrete data, though from a lightly hunted population, so this estimate
    # is also likely a slight underestimate, but for the deep coalescent simulations
    # we also expect that there has been some predation in nature in the past.
    generation_time=3.6,
    # Zhang et al. (2022) on page 1045 write:
    # "... Ne of pigs (the maximum Ne was ~4x10^4; see Figure S6D) ..."
    # and on page 1048 write:
    # "The new estimated mutation rate also revealed a maximum Ne of
    # 2.7x10^5 in pigs, ~6 times larger than that estimated previously
    # [2,3] (Figure 3A, Figure S9)"
    population_size=270_000,
    citations=[
        _ZhangEtAl,  # POP_SIZE
        _ServantyEtAl,  # GEN_TIME
    ],
)

stdpopsim.register_species(_species)
