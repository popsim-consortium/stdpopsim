import stdpopsim
import math

_species = stdpopsim.get_species("DroMel")

###########################################################
#
# DFEs
#
###########################################################


def _HuberDFE():
    id = "Gamma_H17"
    description = "Deleterious Gamma DFE"
    long_description = """
    Deleterious, gamma-distributed DFE estimated from a D. melanogaster SFS in
    Huber et al. (2017).  DFE parameters are based on the "full" model
    described in Table S2, for which singletons were excluded and a recent
    mutation rate estimate is used (3e-9, Keightley 2014).
    """
    # https://doi.org/10.1073/pnas.1619508114
    citations = [
        stdpopsim.Citation(
            author="Huber et al.",
            year=2017,
            doi="https://doi.org/10.1073/pnas.1619508114",
            reasons={stdpopsim.CiteReason.DFE},  # include the dfe_model reason
        )
    ]
    neutral = stdpopsim.MutationType()
    gamma_shape = 0.33  # shape
    gamma_scale = 1.2e-3  # scale
    gamma_mean = gamma_shape * gamma_scale  # expected value
    h = 0.5  # dominance coefficient
    negative = stdpopsim.MutationType(
        dominance_coeff=h,
        distribution_type="g",  # gamma distribution
        # (1+s for homozygote in SLiM versus 1+2s in dadi)
        distribution_args=[-2 * gamma_mean, gamma_shape],
    )
    # p. 2 in supplement says that the total sequence length of synonymous sites LS
    # related to the total sequence length of nonsynonymous sites LNS by LNS = 2.85 * LS;
    # so this is 1 / (1 + 2.85) = 0.2597402597402597
    prop_synonymous = 0.26
    return stdpopsim.DFE(
        id=id,
        description=description,
        long_description=long_description,
        mutation_types=[neutral, negative],
        proportions=[prop_synonymous, 1 - prop_synonymous],
        citations=citations,
    )


_species.add_dfe(_HuberDFE())


def _ZhenDFE():
    id = "GammaPos_H17"
    description = "Deleterious Gamma DFE with fixed-s beneficials"
    long_description = """
    The DFE estimated from D.melanogaster-simulans data estimated in Zhen et al
    (2021, https://dx.doi.org/10.1101/gr.256636.119).
    This uses the demographic model and deleterious-only DFE from
    Huber et al (2017), and then the number of nonsynonymous differences
    to simulans to estimate a proportion of nonsynonmyous differences
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
    neutral = stdpopsim.MutationType()
    gamma_shape = 0.33  # shape
    gamma_scale = 6.01e-4
    gamma_mean = gamma_shape * gamma_scale
    h = 0.5  # dominance coefficient
    negative = stdpopsim.MutationType(
        dominance_coeff=h,
        distribution_type="g",  # gamma distribution
        distribution_args=[round(-2 * gamma_mean, 7), gamma_shape],
    )
    # p. 2 in supplement says that the total sequence length of synonymous sites LS
    # related to the total sequence length of nonsynonymous sites LNS
    # by LNS = 2.85 * LS
    # so this is 1 / (1 + 2.85) = 0.2597402597402597
    prop_synonymous = 0.26

    sel_coeff = 10 ** (-4.801)
    prop_beneficials = 6.75e-4
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
            (1 - prop_synonymous) * (1 - prop_beneficials),
            (1 - prop_synonymous) * prop_beneficials,
        ],
        citations=citations,
    )


_species.add_dfe(_ZhenDFE())


def _RagsdaleDFE():
    id = "LognormalPlusPositive_R16"
    description = "Deleterious log-normal and beneficial mixed DFE"
    long_description = """
    Estimated DFE containing deleterious and beneficial mutations, estimated
    by Ragsdale et al. (2016), with a log-normal distribution for deleterious
    mutations and a single selection coefficient for positive mutations.  The
    DFE was inferred in scaled units from the triallelic frequency spectrum in
    D. melanogaster, here scaled to real units using an effective population
    size of Ne=2.8e6 (Huber et al 2017).
    """
    # DFE parameters are given in Table S1, with deleterious mutations drawn from a
    # log-normal distribution and a point mass of positive selection.
    # Uses only nonsynonymous changes, and assumes a ratio of 2.5 nonsynonymous
    # to 1 synonymous mutation. Results are given as scaled parameters, so that
    # S = 2*N*s, so an estimate of Ne is required to convert the scaled
    # selection coefficients to unscaled coefficients. Because the original
    # study did not report an estimated effective population size, we use the
    # estimate from Huber et al. (2017) of Ne=2.8e6, which was estimated from
    # observed genome-wide synonymous mutation diversity and assuming a
    # mutation rate of 3e-9.
    citations = [
        stdpopsim.Citation(
            author="Ragsdale et al.",
            year=2016,
            doi="https://doi.org/10.1534/genetics.115.184812",
            reasons={stdpopsim.CiteReason.DFE},
        )
    ]
    # the effective population size is used to rescale genetic to physical units
    # and Ne here comes from Huber et al (2017), which used the same dataset to
    # estimate Ne and a different parameterization of the DFE
    Ne = 2.8e6

    # mutation proportions
    p_syn = 1 / (1 + 2.5)  # proportion of neutral
    p_non = 1 - p_syn
    p_positive = 0.0079 * p_non
    p_negative = (1 - 0.0079) * p_non

    # selection strengths
    S_positive = 39.9
    s_positive = S_positive / 2 / Ne
    # double for DaDi
    s_positive = 2 * s_positive
    meanlog = 5.42 - math.log(2 * Ne)
    sdlog = 3.36
    h = 0.5
    # double for DaDi, on the log scale
    meanlog = meanlog + math.log(2)
    synonymous = stdpopsim.MutationType()
    negative = stdpopsim.MutationType(
        dominance_coeff=h,
        distribution_type="ln",  # negative log-normal distribution
        distribution_args=[meanlog, sdlog],
    )
    positive = stdpopsim.MutationType(
        dominance_coeff=h,
        distribution_type="f",
        distribution_args=[s_positive],
    )

    return stdpopsim.DFE(
        id=id,
        description=description,
        long_description=long_description,
        mutation_types=[synonymous, negative, positive],
        proportions=[p_syn, p_negative, p_positive],
        citations=citations,
    )


_species.add_dfe(_RagsdaleDFE())
