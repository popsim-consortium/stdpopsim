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
    """
    citations = [
        stdpopsim.Citation(
            author="Kim et al.",
            year=2017,
            doi="https://doi.org/10.1534/genetics.116.197145",
            reasons={stdpopsim.CiteReason.DFE},  # include the dfe_model reason
        )
    ]
    neutral = stdpopsim.ext.MutationType()
    gamma_shape = 0.186  # shape
    gamma_mean = -0.01314833  # expected value
    h = 0.5  # dominance coefficient
    negative = stdpopsim.ext.MutationType(
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
