import stdpopsim

_species = stdpopsim.get_species("AraTha")

###########################################################
#
# DFEs
#
###########################################################


def _HuberDFE():
    id = "GammaAdditive_H18"
    description = "Deleterious additive Gamma DFE"
    long_description = """
    Additive, deleterious DFE from Huber et al. (2018) for Arabidopsis,
    estimated from the SFS of A. lyrata as a Gamma distribution
    of deleterious effects. Parameters are from Supplementary Table 4,
    the "genome-wide, additive-only model for A. lyrata". The DFE for
    A. lyrata (rather than A. thaliana) is provided due to challenges
    with simulating selfing.
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
    gamma_scale = -0.0004  # scale
    gamma_mean = gamma_shape * gamma_scale
    h = 0.5  # dominance coefficient
    negative = stdpopsim.MutationType(
        dominance_coeff=h,
        distribution_type="g",  # gamma distribution
        # extra factor of 2 is to convert dadi to SLiM
        # (1+s for homozygote in SLiM versus 1+2s in dadi)
        distribution_args=[-2 * gamma_mean, gamma_shape],
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
