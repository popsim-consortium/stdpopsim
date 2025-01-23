import msprime
import numpy as np
import stdpopsim

_species = stdpopsim.get_species("MusMus")

# gen time and mut rate for all three models
spec_generation_time = 1
spec_mutation_rate = 5.7e-9

# Parameters are from Fujiwara et al. 2022 (Figure 3)
# The values themselves were provided by the authors
# directly to Peter Fields who implemented the model.
# The values were visually compared to the curves,
# and after a few tweaks, they were confirmed.
# Here, we use the same values in a separate
# implementation of the three demographic models


def QC_DomesticusEurope():
    # Domesticus model (blue curve in Fig 3 of
    # Fujiwara et al. (2022)
    id = "QC-DomesticusEurope_1F22"
    pop_id = "M_musculus_domesticus"
    time_n_size = np.array(
        [
            (0, 2040),
            (83, 3844),
            (180, 90428),
            (291, 145603),
            (420, 111242),
            (570, 115399),
            (743, 147212),
            (943, 159142),
            (1175, 136620),
            (1443, 97250),
            (1754, 58488),
            (2114, 33028),
            (2530, 18939),
            (3012, 11758),
            (3570, 8463),
            (4216, 7480),
            (4964, 8332),
            (5829, 11240),
            (6831, 16490),
            (7991, 23419),
            (9334, 29931),
            (10889, 34163),
            (12688, 36886),
            (14772, 41195),
            (17183, 50557),
            (19975, 67337),
            (23207, 90926),
            (26948, 115426),
            (31279, 131016),
            (36292, 132063),
            (42096, 121751),
            (48815, 107067),
            (56592, 93046),
            (65596, 81892),
            (76019, 74185),
            (88084, 69939),
            (102052, 69317),
            (118221, 73097),
            (136938, 82953),
            (158606, 101471),
            (183689, 131392),
            (212726, 173264),
            (246340, 222951),
            (285254, 271935),
            (330300, 309961),
            (382446, 327217),
            (442812, 316861),
            (512693, 279833),
            (593589, 227037),
            (687237, 173594),
            (795646, 131050),
            (921140, 98811),
            (1066418, 98811),
            (1234595, 133912),
            (1429281, 133912),
            (1654653, 133912),
        ]
    )

    model = msprime.Demography()
    model.add_population(
        name=pop_id, description=pop_id, initial_size=time_n_size[0][1]
    )
    for j in range(1, time_n_size.shape[0]):
        time, size = time_n_size[j, :]
        model.add_population_parameters_change(
            time,
            initial_size=size,
            population=pop_id,
        )

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=spec_generation_time,
        mutation_rate=spec_mutation_rate,
        model=model,
    )


def QC_MusculusKorea():
    # Musculus Korean population (red curve in Fig 3 of Fujiwara et al. (2022)
    id = "QC-MusculusKorea_1F22"
    pop_id = "M_musculus_musculus"
    time_n_size = np.array(
        [
            (0, 179912),
            (35, 8931),
            (76, 8035),
            (123, 9029),
            (177, 9960),
            (240, 12104),
            (313, 16254),
            (398, 25527),
            (495, 42715),
            (609, 61935),
            (740, 68111),
            (891, 55959),
            (1067, 36220),
            (1270, 20382),
            (1505, 11222),
            (1778, 6695),
            (2093, 4605),
            (2458, 3751),
            (2881, 3643),
            (3370, 4177),
            (3936, 5506),
            (4591, 7990),
            (5350, 12072),
            (6229, 17741),
            (7246, 23546),
            (8423, 26648),
            (9785, 25399),
            (11363, 21219),
            (13189, 16747),
            (15303, 13588),
            (17750, 12259),
            (20583, 13023),
            (23863, 16339),
            (27659, 22556),
            (32054, 30806),
            (37142, 38441),
            (43031, 42857),
            (49849, 43874),
            (57741, 43467),
            (66878, 43933),
            (77455, 47001),
            (89698, 54304),
            (103872, 67725),
            (120280, 88494),
            (139274, 116547),
            (161262, 151909),
            (186716, 194969),
            (216182, 245823),
            (250293, 302950),
            (289781, 359368),
            (335491, 400867),
            (388409, 407105),
            (449667, 407105),
            (520579, 152757),
            (602670, 152757),
            (697702, 152757),
        ]
    )

    model = msprime.Demography()
    model.add_population(
        name=pop_id, description=pop_id, initial_size=time_n_size[0][1]
    )
    for j in range(1, time_n_size.shape[0]):
        time, size = time_n_size[j, :]
        model.add_population_parameters_change(
            time,
            initial_size=size,
            population=pop_id,
        )

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=spec_generation_time,
        mutation_rate=spec_mutation_rate,
        model=model,
    )


def QC_CastaneusIndia():
    # Castaneus Indian model (green curve in Fig 3 of Fujiwara et al. (2022)
    id = "QC-CastaneusIndia_1F22"
    pop_id = "M_musculus_castaneus"
    time_n_size = np.array(
        [
            (0, 64853),
            (887, 5064712),
            (1886, 938111),
            (3011, 291323),
            (4279, 141377),
            (5709, 86633),
            (7319, 60911),
            (9133, 47736),
            (11178, 41845),
            (13481, 41297),
            (16077, 45372),
            (19002, 53747),
            (22297, 66260),
            (26011, 83216),
            (30195, 105533),
            (34909, 134861),
            (40221, 173419),
            (46207, 222464),
            (52951, 279816),
            (60551, 338426),
            (69114, 388298),
            (78762, 420389),
            (89634, 429975),
            (101884, 418882),
            (115687, 394606),
            (131239, 365787),
            (148764, 338818),
            (168510, 317235),
            (190760, 302369),
            (215830, 294117),
            (244079, 291560),
            (275907, 293233),
            (311772, 297138),
            (352184, 300840),
            (397719, 301783),
            (449026, 297772),
            (506839, 287572),
            (571979, 271326),
            (645379, 250544),
            (728084, 227605),
            (821274, 205098),
            (926277, 185260),
            (1044593, 169657),
            (1177907, 159139),
            (1328123, 154098),
            (1497384, 155027),
            (1688102, 163295),
            (1903000, 182054),
            (2145140, 217017),
            (2417965, 275733),
            (2725404, 360433),
            (3071789, 463464),
            (3462105, 463464),
            (3901912, 344802),
            (4397456, 344802),
            (4955842, 344802),
        ]
    )

    model = msprime.Demography()
    model.add_population(
        name=pop_id, description=pop_id, initial_size=time_n_size[0][1]
    )
    for j in range(1, time_n_size.shape[0]):
        time, size = time_n_size[j, :]
        model.add_population_parameters_change(
            time,
            initial_size=size,
            population=pop_id,
        )

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=spec_generation_time,
        mutation_rate=spec_mutation_rate,
        model=model,
    )


_species.get_demographic_model("DomesticusEurope_1F22").register_qc(
    QC_DomesticusEurope()
)

_species.get_demographic_model("MusculusKorea_1F22").register_qc(QC_MusculusKorea())

_species.get_demographic_model("CastaneusIndia_1F22").register_qc(QC_CastaneusIndia())


def QC_GammaB21():
    id = "QC-GammaB21"
    # M. m. castaneus from https://www.biorxiv.org/content/10.1101/2021.06.10.447924v2
    # 0-fold site model with parameters from Tables S1, S2
    neutral = stdpopsim.MutationType()
    gamma_shape = 0.18617976
    # Using polyDFE, the DFE is estimated in terms of
    # the scaled selection coefficient for deleterious mutations,
    # 2Nesd, where sd is the reduction in fitness experienced by
    # an individual homozygous for the mutation
    gamma_mean = -50044.583
    Ne = 420000  # effective population size from methods
    gamma_mean /= Ne
    h = 0.5
    negative = stdpopsim.MutationType(
        dominance_coeff=h,
        distribution_type="g",
        # PolyDFE gives two times the selection coefficient on a homozygote,
        # so divide mean by two
        distribution_args=[gamma_mean/2, gamma_shape],
    )
    tstv = 3.3833
    prop_nonsynonymous = 1 / (1 + tstv)
    return stdpopsim.DFE(
        id=id,
        description=id,
        long_description=id,
        mutation_types=[neutral, negative],
        proportions=[1 - prop_nonsynonymous, prop_nonsynonymous],
    )


_species.get_dfe("Gamma_B21").register_qc(QC_GammaB21())
