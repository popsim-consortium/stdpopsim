import math
import stdpopsim

_species = stdpopsim.get_species("HomSap")

###########################################################
#
# DFEs
#
###########################################################


def _KimDFE():
    id = "Gamma_K17"
    description = "Deleterious Gamma DFE"
    long_description = """
    Deleterious, gamma-distributed DFE estimated using 'Fit∂a∂i' from the nonsynonymous
    SFS in exons from a dataset of European humans by Kim et al. (2017).
    """
    # The DFE was inferred assuming synonymous variants are neutral and a relative
    # mutation rate ratio of 2.31 nonsynonymous to 1 synonymous mutation
    # (Huber et al. 2016). Gamma parameters are given as the shape and the mean
    # selection coefficient (E[s]) escaled, such that E[s] = -shape*scale*2/(2Na),
    # where Na is the ancestral population size (12378) as in Kim et al. (2017).
    citations = [
        stdpopsim.Citation(
            author="Kim et al.",
            year=2017,
            doi="https://doi.org/10.1534/genetics.116.197145",
            reasons={stdpopsim.CiteReason.DFE},  # include the dfe_model reason
        )
    ]
    neutral = stdpopsim.MutationType()
    Na = 12378
    gamma_scale = 875
    gamma_shape = 0.186  # shape
    # Extra factor of 2 in mean is to account for difference between tools
    # in how fitness is defined
    # (1+s for homozygote in SLiM versus 1+2s in dadi)
    gamma_mean = (-gamma_shape * gamma_scale * 2) / (2 * Na)  # expected value
    h = 0.5  # dominance coefficient
    negative = stdpopsim.MutationType(
        dominance_coeff=h,
        distribution_type="g",  # gamma distribution
        distribution_args=[gamma_mean, gamma_shape],
    )

    return stdpopsim.DFE(
        id=id,
        description=description,
        long_description=long_description,
        mutation_types=[neutral, negative],
        proportions=[0.3, 0.7],
        citations=citations,
    )


_species.add_dfe(_KimDFE())


def _HuberDFE():
    id = "Gamma_H17"
    description = "Deleterious Gamma DFE"
    long_description = """
    Deleterious, gamma-distributed DFE estimated from the SFS of YRI samples in
    the 1000 genomes project by Huber et al. (2017). DFE parameters are from
    the "full" model in Table S2, in which all data (none of the listed
    filters) were used.
    """
    # [this was a different filtering scheme than the Huber et al DroMel DFE ]
    citations = [
        stdpopsim.Citation(
            author="Huber et al.",
            year=2017,
            doi="https://doi.org/10.1073/pnas.1619508114",
            reasons=stdpopsim.CiteReason.DFE,
        )
    ]
    neutral = stdpopsim.MutationType()
    gamma_shape = 0.19  # shape
    # Extra factor of 2 in mean is to account for difference between tools
    # in how fitness is defined
    # (1+s for homozygote in SLiM versus 1+2s in dadi)
    gamma_mean = -0.014 * 2  # expected value
    h = 0.5  # dominance coefficient
    negative = stdpopsim.MutationType(
        dominance_coeff=h,
        distribution_type="g",  # gamma distribution
        distribution_args=[gamma_mean, gamma_shape],
    )

    return stdpopsim.DFE(
        id=id,
        description=description,
        long_description=long_description,
        mutation_types=[neutral, negative],
        proportions=[
            0.3,
            0.7,
        ],  # [0.3 and 0.7 were used in Xin's analysis,
        #  but I couldn't find these values in the Huber paper]
        citations=citations,
    )


_species.add_dfe(_HuberDFE())


def _HuberLogNormalDFE():
    id = "LogNormal_H17"
    description = "Deleterious Log-normal DFE"
    long_description = """
    Deleterious, log-normal-distributed DFE estimated from the SFS of YRI
    samples in the 1000 genomes project by Huber et al. (2017).
    Parameters as shown for the "full" model in Table S3.
    """
    citations = [
        stdpopsim.Citation(
            author="Huber et al.",
            year=2017,
            doi="https://doi.org/10.1073/pnas.1619508114",
            reasons={stdpopsim.CiteReason.DFE},
        )
    ]
    neutral = stdpopsim.MutationType()
    # The log(2) term doubles all selection coefficients, to convert
    # from dadi convention to SLiM convention.
    mulog = -7.37 + math.log(2)
    sigmalog = 4.58
    h = 0.5  # dominance coefficient
    negative = stdpopsim.MutationType(
        dominance_coeff=h,
        distribution_type="ln",  # negative logNormal distribution
        distribution_args=[mulog, sigmalog],
    )

    return stdpopsim.DFE(
        id=id,
        description=description,
        long_description=long_description,
        mutation_types=[neutral, negative],
        proportions=[0.3, 0.7],
        citations=citations,
    )


_species.add_dfe(_HuberLogNormalDFE())


def _KyriazisDFE():
    id = "Mixed_K23"
    description = "Deleterious Gamma DFE with additional lethals"
    long_description = """
    The DFE estimated from human data recommended in Kyriazis et al.  (2023),
    for general use.  This model is similar to the Kim et al. (2017) DFE based
    on human genetic data, modified to include the dominance distribution from
    Henn et al. (2016).  The model is also augmented with an additional
    proportion of 0.3% of recessive lethals, based on the analysis of
    Wade et al. (2023).
    """
    citations = [
        stdpopsim.Citation(
            author="Kyriazis et al.",
            year=2023,
            doi="https://doi.org/10.1086/726736",
            reasons={stdpopsim.CiteReason.DFE},
        )
    ]
    neutral = stdpopsim.MutationType()
    gamma_mean = -0.0131
    gamma_shape = 0.186
    coefs = [0.45, 0.2, 0.05, 0]
    breaks = [0.001, 0.01, 0.1]
    gamma = stdpopsim.MutationType(
        dominance_coeff_list=coefs,
        dominance_coeff_breaks=breaks,
        distribution_type="g",  # gamma distribution
        distribution_args=[gamma_mean, gamma_shape],
    )
    lethal = stdpopsim.MutationType(
        distribution_type="f",  # fixed value
        distribution_args=[-1],  # fitness in SLiM for homozygotes is multiiplied by 1+s
        dominance_coeff=0,
    )
    proportion_deleterious = 2.31 / (1 + 2.31)
    lethal_prop = proportion_deleterious * 0.003  # 0.3% lethals
    gamma_prop = proportion_deleterious - lethal_prop
    neutral_prop = 1 - proportion_deleterious
    return stdpopsim.DFE(
        id=id,
        description=description,
        long_description=long_description,
        mutation_types=[neutral, gamma, lethal],
        proportions=[neutral_prop, gamma_prop, lethal_prop],
        citations=citations,
    )


_species.add_dfe(_KyriazisDFE())


def _RodriguesDFE():
    id = "PosNeg_R24"
    description = "Deleterious Gamma and Beneficial Exponential DFE"
    long_description = """
    The best-fitting DFE simulated by Rodrigues et al (2024),
    from among 57 simulated scenarios with varying amounts of positive
    and negative selection. Fit was based on similarity of diversity
    and divergence across the great ape clade in simulations; the shape
    of the DFE was not inferred, only the proportion of positive and
    negative mutations; the shape of the deleterious portion of the DFE
    was obtained from Castellano et al (2019),
    https://doi.org/10.1534/genetics.119.302494.
    """
    citations = [
        stdpopsim.Citation(
            author="Rodrigues et al.",
            year=2024,
            doi="https://doi.org/10.1093/genetics/iyae006",
            reasons={stdpopsim.CiteReason.DFE},
        ),
        stdpopsim.Citation(
            author="Castellano et al.",
            year=2019,
            doi="https://doi.org/10.1534/genetics.119.302494",
            reasons={stdpopsim.CiteReason.DFE},
        ),
    ]
    neutral = stdpopsim.MutationType()
    # parameters from row 32 in Table S1, based on
    # identification of rates in Figure 9
    neg_mean = -3e-2
    neg_shape = 0.16
    negative = stdpopsim.MutationType(
        dominance_coeff=0.5,
        distribution_type="g",  # gamma distribution
        distribution_args=[neg_mean, neg_shape],
    )
    pos_mean = 1e-2
    positive = stdpopsim.MutationType(
        dominance_coeff=0.5,
        distribution_type="e",  # exponential distribution
        distribution_args=[pos_mean],
    )
    total_rate = 2e-8
    pos_rate = 1e-12
    neg_rate = 1.2e-8
    prop_pos = pos_rate / total_rate
    prop_neg = neg_rate / total_rate
    neutral_prop = 1 - prop_pos - prop_neg
    return stdpopsim.DFE(
        id=id,
        description=description,
        long_description=long_description,
        mutation_types=[neutral, negative, positive],
        proportions=[neutral_prop, prop_neg, prop_pos],
        citations=citations,
    )


_species.add_dfe(_RodriguesDFE())


def _ZhenDFE():
    id = "GammaPos_Z21"
    description = "Deleterious Gamma DFE with fixed-s beneficials"
    long_description = """
    The DFE estimated from human-chimp data estimated in Zhen et al
    (2021, https://dx.doi.org/10.1101/gr.256636.119).
    This uses the demographic model and deleterious-only DFE from
    Huber et al (2017), and then the number of nonsynonymous differences
    to chimpanzee to estimate a proportion of nonsynonmyous differences
    and (single) selection coefficient. So, this is the "Gamma_H17" DFE,
    with some proportion of positive selection with fixed s.
    """
    citations = [
        stdpopsim.Citation(
            author="Zhen et al.",
            year=2021,
            doi="https://dx.doi.org/10.1101/gr.256636.119",
            reasons={stdpopsim.CiteReason.DFE},
        )
    ]
    # modified from _HuberDFE() above
    neutral = stdpopsim.MutationType()
    gamma_shape = 0.19  # shape
    # Extra factor of 2 in mean is to account for difference between tools
    # in how fitness is defined
    # (1+s for homozygote in SLiM versus 1+2s in dadi)
    gamma_mean = -0.014 * 2  # expected value
    negative = stdpopsim.MutationType(
        dominance_coeff=0.5,
        distribution_type="g",  # gamma distribution
        distribution_args=[gamma_mean, gamma_shape],
    )
    # p. 2 in supplement says that the total sequence length of synonymous sites LS
    # related to the total sequence length of nonsynonymous sites LNS
    # by LNS = 2.31 * LS
    # so, this is 1 / (1 + 2.31) = 0.3021148036253776
    prop_synonymous = 0.3

    sel_coeff = 10 ** (-3.949)
    prop_beneficial = 1.55e-2
    positive = stdpopsim.MutationType(
        dominance_coeff=0.5,
        distribution_type="f",
        distribution_args=[sel_coeff],
    )

    return stdpopsim.DFE(
        id=id,
        description=description,
        long_description=long_description,
        mutation_types=[neutral, negative, positive],
        proportions=[
            prop_synonymous,
            (1 - prop_synonymous) * (1 - prop_beneficial),
            (1 - prop_synonymous) * prop_beneficial,
        ],
        citations=citations,
    )


_species.add_dfe(_ZhenDFE())
