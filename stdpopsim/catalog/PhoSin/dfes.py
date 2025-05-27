import stdpopsim

_species = stdpopsim.get_species("PhoSin")

###########################################################
#
# DFEs
#
###########################################################


def _RobinsonDFE():
    id = "Gamma_R22"
    description = "Deleterious Gamma DFE"
    long_description = """
    Gamma-distributed deleterious DFE inferred by Robinson et al. (2022),
    and used in their simulations, from the synonymous and nonsynonymous SFS from
    Vaquita data, following the methods of Huber et al (2017).
    """
    # The DFE was inferred assuming synonymous variants are neutral and a relative
    # mutation rate ratio of 2.31 nonsynonymous to 1 synonymous mutation
    # (Huber et al. 2017). The inference procedure is described in pages 9-10 of
    # the supplementary material of Robinson et al. (2022), and values used in
    # SLiM simulations based on this DFE are specified in pages 11-12.
    citations = [
        stdpopsim.Citation(
            author="Robinson et al.",
            year=2022,
            doi="https://doi.org/10.1126/science.abm1742",
            reasons={stdpopsim.CiteReason.DFE},  # include the dfe_model reason
        )
    ]
    # simple neutral mutation type
    neutral = stdpopsim.MutationType()

    ################################################
    # properties for negatively selected mutations #
    ################################################
    gamma_shape = 0.131  # shape inferred by fitDadi (see top of page 11 in supp)
    # gamma_scale = 0.098  # scale inferred by fitDadi (see top of page 11 in supp)
    # gamma_scale is not used here, since we directly used the gamma_mean reported
    # in the bottom of page 12 in the supp (see below)
    # The inferred shape and scale were used to compute the mean of the gamma
    # distribution using this formula:
    # gamma_mean = (-gamma_shape * gamma_scale * 2) / (2 * Na)
    # [ the extra factor of 2 is to account for difference between how fitness
    #   is defined in fitDadi and in SLiM  ]
    # Since the inferred value of Na is not reported in the supplement, we use
    # here the value of the mean as reported in the supp info
    gamma_mean = -0.0257  # mean of gamma          (see bottom of page 12 in supp)

    # dominance coefficient assumed in fitDadi inference was fully additive (h=0.5).
    # However, in the simulations described by Robinson et al. (2022), they used
    # an ad-hoc inverse relationship between h and s,
    # as described in the top of page 12 in the supp:
    # " For dominance coefficients (h), we assumed an inverse relationship between
    #   h and s (hs relationship) with:
    #   h = 0.0  for very strongly deleterious mutations (s < -0.1),
    #   h = 0.01 for strongly deleterious mutations (-0.1 ≤ s < -0.01),
    #   h = 0.1  for moderately deleterious mutations (-0.01 ≤ s < -0.001), and
    #   h = 0.4  for weakly deleterious mutations (s > -0.001).                  "
    negative = stdpopsim.MutationType(
        dominance_coeff_list=[0.0, 0.01, 0.1, 0.4],
        dominance_coeff_breaks=[-0.1, -0.01, -0.001],
        distribution_type="g",  # gamma distribution
        distribution_args=[gamma_mean, gamma_shape],
    )

    return stdpopsim.DFE(
        id=id,
        description=description,
        long_description=long_description,
        mutation_types=[neutral, negative],
        proportions=[0.3, 0.7],  # due to assumption of ratio of 2.31 nonsynonymous
        # to 1 synonymous mutation (Huber et al. 2017).
        citations=citations,
    )


_species.add_dfe(_RobinsonDFE())
