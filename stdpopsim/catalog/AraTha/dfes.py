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
    Return neutral and negative MutationType()s representing an Arabidopsis DFE.
    From Huber et al. (2018), https://doi.org/10.1038/s41467-018-05281-7.
    Gamma parameters are based on the additive-only model for A. thaliana described in
    Supplementary Table 4.
    """
    citations = [
        stdpopsim.Citation(
            author="Huber et al.",
            year=2018,
            doi="https://doi.org/10.1038/s41467-018-05281-7",
            reasons="to be defined",  # include the dfe_model reason
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
