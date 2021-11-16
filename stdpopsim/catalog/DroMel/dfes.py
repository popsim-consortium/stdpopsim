import stdpopsim

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
    Return neutral and negative MutationType()s representing a drosophila DFE.
    Huber et al. (2017), https://doi.org/10.1073/pnas.1619508114.
    DFE parameters are based on the Full model described in Table S2, in which
    singletons are excluded and a recent mutation rate estimate is used
    (mu=3x10e-9, Keightley 2014).
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
    gamma_shape = 0.33  # shape
    gamma_mean = -3.96e-04  # expected value
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
        proportions=[0.26, 0.74],  # LNS = 2.85 x LS
        citations=citations,
    )


_species.add_dfe(_HuberDFE())
