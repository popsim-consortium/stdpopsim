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
