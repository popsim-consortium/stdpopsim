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
    Return neutral and negative MutationType()s representing a human DFE.
    Kim et al. (2017), https://doi.org/10.1534/genetics.116.197145
    The DFE was inferred assuming synonymous variants are neutral and a relative
    mutation rate ratio of 2.31 nonsynonymous to 1 synonymous mutation
    (Huber et al. 2016). Gamma parameters are given as the shape and the mean
    selection coefficient (E[s]) escaled, such that E[s] = -shape*scale*2/(2Na),
    where Na is the ancestral population size (12378) as in Kim et al. (2017).
    """
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
    Return neutral and negative MutationType()s representing a Homo sapiens DFE.
    Huber et al. (2017), https://doi.org/10.1073/pnas.1619508114.
    DFE parameters are based on the Full Model described in Table S2, in which
    All Data (none of the listed filters) were used.
    """  # [this was a different filtering scheme than the Huber et al DroMel DFE ]
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
    Return neutral and negative MutationType()s representing a human DFE.
    Huber et al. (2017), https://doi.org/10.1073/pnas.1619508114.
    DFE parameters are based off the Full model in Table S3, Using recent
    mutation rate estimates.
    Log-normal distribution parameters were given as the
    mean and standard deviation of the log.
    """
    citations = [
        stdpopsim.Citation(
            author="Huber et al.",
            year=2017,
            doi="https://doi.org/10.1073/pnas.1619508114",
            reasons={stdpopsim.CiteReason.DFE},  # include the dfe_model reason
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
    The DFE estimated from human data recommended in Kyriazis et al.
    (2023), https://doi.org/10.1086/726736, for general use.
    This model is similar to the Kim et al. (2017) DFE based on human
    genetic data, modified to include the dominance distribution from
    Henn et al. (2016).
    The model is also augmented with an additional proportion of 0.3% of
    recessive lethals, based on the analysis of Wade et al. (2023).
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
        proportions=[neutral_prop, prop_pos, prop_neg],
        citations=citations,
    )


_species.add_dfe(_RodriguesDFE())
