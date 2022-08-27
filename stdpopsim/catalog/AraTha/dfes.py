import stdpopsim

_species = stdpopsim.get_species("AraTha")

###########################################################
#
# DFEs
#
###########################################################


def _HuberDFE():
    id = "Gamma_H18"
    description = "Deleterious Gamma DFE"
    long_description = """
    Return negative (no neutral, in this case) MutationType()s
    representing an Arabidopsis DFE. From Huber et al. (2018).
    Gamma parameters are based on Supplementary Table 4,
    the genome-wide, additive-only model for A. LYRATA- due to
    challenges with simulating with selfing for A. thaliana.
    The Supplementary Table 4 DFEs are not noted to contain
    any neutral mutations (in contrast, Supplementary Table 3
    notes neutral proportions for two other DFEs), and in the
    main text including neutral mutations seems like a
    supplementary analysis rather than the main strategy.
    """
    citations = [
        stdpopsim.Citation(
            author="Huber et al.",
            year=2018,
            doi="https://doi.org/10.1038/s41467-018-05281-7",
            reasons=stdpopsim.CiteReason.DFE,
        )
    ]
    neutral = stdpopsim.MutationType()
    gamma_shape = 0.27  # shape
    gamma_mean = -0.0004  # expected value
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
            0,
            1.0,
        ],
        citations=citations,
    )


_species.add_dfe(_HuberDFE())
