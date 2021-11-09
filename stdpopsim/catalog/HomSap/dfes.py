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
            reasons="to be defined",  # include the dfe_model reason
        )
    ]
    neutral = stdpopsim.MutationType()
    gamma_shape = 0.186  # shape
    gamma_mean = -0.01314833  # expected value
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
    """ # [this was a different filtering scheme than the Huber et al DroMel DFE ]
    citations = [
        stdpopsim.Citation(
            author="Huber et al.",
            year=2017,
            doi="https://doi.org/10.1073/pnas.1619508114",
            reasons="to be defined",  # include the dfe_model reason
        )
    ]
    neutral = stdpopsim.MutationType()
    gamma_shape = 0.19  # shape
    gamma_mean = -0.014  # expected value      
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
        proportions=[0.3, 0.7], # [0.3 and 0.7 were used in Xin's analysis, but I couldn't find these values in the Huber paper]
        citations=citations,
    )


_species.add_dfe(_HuberDFE()) 







    
