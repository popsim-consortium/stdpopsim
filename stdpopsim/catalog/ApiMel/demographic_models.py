# import msprime
# import numpy as np
import stdpopsim

_species = stdpopsim.get_species("ApiMel")


# Population definitions that are reused - the four major lineager
_apis_cerana = stdpopsim.Population(id="Acer", description="Apis cerana population")

_C_lineage = stdpopsim.Population(id="C", description="Apis mellifera C. lineage")
_M_lineage = stdpopsim.Population(id="M", description="Apis mellifera M. lineage")
_A_lineage = stdpopsim.Population(id="A", description="Apis mellifera A. lineage")
_O_lineage = stdpopsim.Population(id="O", description="Apis mellifera O. lineage")

_WallbergEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1038/ng.3077", year=2014, author="Wallberg et al."
)


# def model1():
#     """
#     Demographic model modelled by Wallberg et al., 2014
#     """
#     times = [3e04, 3e05, 6.6e06]
#     sizes = [500,000]
