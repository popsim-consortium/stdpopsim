import msprime
import stdpopsim
import numpy as np

_species = stdpopsim.get_species("SusScr")

_ZhangEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1016/j.gpb.2022.02.001",
    year=2022,
    author="Zhang et al.",
    reasons={stdpopsim.CiteReason.DEM_MODEL},
    # Figure 3A Demographic history of four wild boar populations using PSMC.
)


def _WildBoar_4Z22_SMW():
    id = "WildBoar_4Z22_SMW"
    description = "Piecewise size model for Sumatran wild boar (Zhang et al. 2022)"
    long_description = """
    This demographic model is a piecewise size model for Sumatran wild boar
    from Zhang et al. (2022).
    """
    populations = [
        stdpopsim.Population(id="SMW", description="Sumatran wild boar"),
    ]

    # Recombination rate used in simulations by Zhang et al. (2022).
    recombination_rate = 1.0e-8

    # Mutation rate as reported by Zhang et al. (2022, p. 1042).
    mutation_rate = 3.6e-9

    # Arrays of time intervals and effective population sizes (Ne)
    # extracted from Figure 3A of Zhang et al. (2022).
    times_SMW = np.array(
        [
            2.73927358e02,
            5.78023037e02,
            9.15644600e02,
            1.29042319e03,
            1.70646249e03,
            2.16836360e03,
            2.68112557e03,
            3.25034434e03,
            3.88228734e03,
            4.58381894e03,
            5.36262421e03,
            6.22720899e03,
            7.18699932e03,
            8.25254043e03,
            9.43542214e03,
            1.07485773e04,
            1.22063812e04,
            1.38247512e04,
            1.56213455e04,
            1.76158620e04,
            1.98300377e04,
            2.22880719e04,
            2.50168500e04,
            2.80461921e04,
            3.14091764e04,
            3.51425621e04,
            3.92871621e04,
            4.38882663e04,
            4.89961137e04,
            5.46665891e04,
            6.09615960e04,
            6.79499271e04,
            7.57079856e04,
            8.43205313e04,
            9.38817004e04,
            1.04495925e05,
            1.16279229e05,
            1.29360367e05,
            1.43882272e05,
            1.60003671e05,
            1.77900673e05,
            1.97768912e05,
            2.19825486e05,
            2.44311369e05,
            2.71494220e05,
            3.01670997e05,
            3.35171512e05,
            3.72361864e05,
            4.13648392e05,
            4.59482304e05,
            5.10364447e05,
            5.66850835e05,
            6.29558737e05,
            6.99173418e05,
            7.76455553e05,
            8.62249706e05,
            9.57493383e05,
        ]
    )

    sizes_SMW = np.array(
        [
            1779489.40969792,
            1779489.40969792,
            13416.26015625,
            13416.26015625,
            4897.87538542,
            4407.55934167,
            4602.38501458,
            4882.16945417,
            5175.03595208,
            5521.00167917,
            5931.53208958,
            6379.94077917,
            6835.21381875,
            7291.75527083,
            7769.6234625,
            8299.03645625,
            8915.17404583,
            9670.28985208,
            10648.26076042,
            11973.82643542,
            13810.97271667,
            16361.98652708,
            19850.04628958,
            24476.4067875,
            30345.97319583,
            37379.65626458,
            45251.64807708,
            53429.50139792,
            61330.81592292,
            68519.9417,
            74824.37462917,
            80318.72721042,
            85281.27919792,
            90121.39207292,
            95327.01915,
            101382.3458125,
            108697.15634583,
            117580.74443333,
            128223.79475417,
            140717.28481667,
            155055.94200417,
            171127.08900625,
            188647.94444792,
            207099.76492292,
            225653.40658958,
            243141.8429,
            258150.88096458,
            269176.84263542,
            274910.67647083,
            274517.54320833,
            267869.51971667,
            255602.00604583,
            238998.64808125,
            219719.2630625,
            177796.99883125,
            177796.99883125,
            177796.99883125,
            177796.99883125,
        ]
    )

    demographic_events = []
    population_configurations = [
        msprime.PopulationConfiguration(
            # SMW
            initial_size=sizes_SMW[0],
            metadata=populations[0].asdict(),
        )
    ]

    for i, t in enumerate(times_SMW):
        curr_time = t
        demographic_events.append(
            msprime.PopulationParametersChange(
                time=curr_time, initial_size=sizes_SMW[i + 1], population_id=3
            )
        )
    citations = [
        _ZhangEtAl,
    ]
    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=3,
        mutation_rate=mutation_rate,
        recombination_rate=recombination_rate,
        population_configurations=population_configurations,
        demographic_events=demographic_events,
    )


_species.add_demographic_model(_WildBoar_4Z22_SMW())


def _WildBoar_4Z22_SCW():
    id = "WildBoar_4Z22_SCW"
    description = "Piecewise size model for South Chinese wild boar (Zhang et al. 2022)"
    long_description = """
    This demographic model is a piecewise size model for South Chinese wild boar
    from Zhang et al. (2022).
    """
    populations = [
        stdpopsim.Population(id="SCW", description="South Chinese wild boar"),
    ]

    # Recombination rate used in simulations by Zhang et al. (2022).
    recombination_rate = 1.0e-8

    # Mutation rate as reported by Zhang et al. (2022, p. 1042).
    mutation_rate = 3.6e-9

    # Arrays of time intervals and effective population sizes (Ne)
    # extracted from Figure 3A of Zhang et al. (2022).
    times_SCW = np.array(
        [
            1.99403694e03,
            4.16635167e03,
            6.53298917e03,
            9.11155438e03,
            1.19207665e04,
            1.49813502e04,
            1.83158132e04,
            2.19486686e04,
            2.59064353e04,
            3.02185290e04,
            3.49163713e04,
            4.00345035e04,
            4.56105867e04,
            5.16856248e04,
            5.83041873e04,
            6.55148549e04,
            7.33708880e04,
            8.19297812e04,
            9.12546004e04,
            1.01413760e05,
            1.12481913e05,
            1.24540399e05,
            1.37677688e05,
            1.51990497e05,
            1.67584008e05,
            1.84572766e05,
            2.03081788e05,
            2.23246787e05,
            2.45215958e05,
            2.69151087e05,
            2.95227777e05,
            3.23637678e05,
            3.54589597e05,
            3.88310839e05,
            4.25049432e05,
            4.65075467e05,
            5.08682658e05,
            5.56191680e05,
            6.07951735e05,
            6.64343225e05,
            7.25780199e05,
            7.92714367e05,
            8.65637555e05,
            9.45085710e05,
            1.03164269e06,
            1.12594450e06,
            1.22868398e06,
            1.34061613e06,
            1.46256370e06,
            1.59542276e06,
            1.74016939e06,
            1.89786767e06,
            2.06967619e06,
            2.25685739e06,
            2.46078689e06,
            2.68296356e06,
            2.92501955e06,
        ]
    )

    sizes_SCW = np.array(
        [
            33850.15878472,
            33850.15878472,
            42415.73746528,
            42415.73746528,
            122471.26215278,
            133968.84173611,
            137971.28927083,
            143358.39802083,
            153983.53072917,
            170674.89909722,
            188417.32638889,
            202479.8775,
            214028.48909722,
            222119.29177083,
            223082.21461806,
            216006.36961806,
            203076.21666667,
            187279.35704861,
            171195.58163194,
            156701.04118056,
            144929.35951389,
            136399.5478125,
            131158.73826389,
            129003.69420139,
            129664.10194444,
            132892.7125,
            138474.47829861,
            146157.35913194,
            155600.39875,
            166325.92413194,
            177722.21965278,
            189072.05152778,
            199623.53322917,
            208686.39548611,
            215704.30020833,
            220336.06829861,
            222488.32677083,
            222294.56111111,
            220072.32861111,
            216231.4453125,
            211220.05697917,
            205477.84118056,
            199411.93979167,
            193390.05072917,
            187720.14885417,
            182639.90072917,
            178297.72260417,
            174756.79166667,
            172003.29138889,
            169978.50152778,
            168579.8009375,
            167699.88868056,
            167230.23815972,
            167076.80784722,
            167462.11069444,
            167462.11069444,
            167462.11069444,
            167462.11069444,
        ]
    )

    demographic_events = []
    population_configurations = [
        msprime.PopulationConfiguration(
            # SCW
            initial_size=sizes_SCW[0],
            metadata=populations[0].asdict(),
        )
    ]

    for i, t in enumerate(times_SCW):
        curr_time = t
        demographic_events.append(
            msprime.PopulationParametersChange(
                time=curr_time, initial_size=sizes_SCW[i + 1], population_id=2
            )
        )
    citations = [
        _ZhangEtAl,
    ]
    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=3,
        mutation_rate=mutation_rate,
        recombination_rate=recombination_rate,
        population_configurations=population_configurations,
        demographic_events=demographic_events,
    )


_species.add_demographic_model(_WildBoar_4Z22_SCW())


def _WildBoar_4Z22_NCW():
    id = "WildBoar_4Z22_NCW"
    description = "Piecewise size model for North Chinese wild boar (Zhang et al. 2022)"
    long_description = """
    This demographic model is a piecewise size model for North Chinese wild boar
    from Zhang et al. (2022).
    """
    populations = [
        stdpopsim.Population(id="NCW", description="North Chinese wild boar")
    ]

    # Recombination rate used in simulations by Zhang et al. (2022).
    recombination_rate = 1.0e-8

    # Mutation rate as reported by Zhang et al. (2022, p. 1042).
    mutation_rate = 3.6e-9

    # Arrays of time intervals and effective population sizes (Ne)
    # extracted from Figure 3A of Zhang et al. (2022).

    times_NCW = np.array(
        [
            1.65039886e03,
            3.45249869e03,
            5.42012314e03,
            7.56889111e03,
            9.91496011e03,
            1.24768215e04,
            1.52742233e04,
            1.83287089e04,
            2.16639760e04,
            2.53058765e04,
            2.92825963e04,
            3.36248347e04,
            3.83663428e04,
            4.35437444e04,
            4.91968946e04,
            5.53699572e04,
            6.21103276e04,
            6.94704279e04,
            7.75071684e04,
            8.62826658e04,
            9.58647814e04,
            1.06327840e05,
            1.17752808e05,
            1.30228014e05,
            1.43850043e05,
            1.58724278e05,
            1.74965977e05,
            1.92700628e05,
            2.12065572e05,
            2.33210892e05,
            2.56299960e05,
            2.81511584e05,
            3.09040912e05,
            3.39100864e05,
            3.71924286e05,
            4.07765032e05,
            4.46900652e05,
            4.89633828e05,
            5.36295431e05,
            5.87246671e05,
            6.42881612e05,
            7.03631119e05,
            7.69965196e05,
            8.42397294e05,
            9.21487897e05,
            1.00784920e06,
            1.10214957e06,
            1.20511881e06,
            1.31755382e06,
            1.44032479e06,
            1.57438195e06,
            1.72076280e06,
            1.88060015e06,
            2.05513113e06,
            2.24570651e06,
            2.45380111e06,
            2.68102546e06,
        ]
    )
    sizes_NCW = np.array(
        [
            96756.31630556,
            96756.31630556,
            23458.80495833,
            23458.80495833,
            29425.85947222,
            33474.83920833,
            37881.70756944,
            44466.247875,
            54119.27695833,
            66767.63725,
            82821.19066667,
            102492.677625,
            124769.92051389,
            147268.88020833,
            167371.86218056,
            183211.86730556,
            193154.92351389,
            195884.55361111,
            191547.07273611,
            181982.01226389,
            169846.56283333,
            157649.26608333,
            147192.93995833,
            139407.98716667,
            134604.18288889,
            132773.80743056,
            133793.61497222,
            137496.10609722,
            143635.95609722,
            151873.229125,
            161743.12776389,
            172654.01822222,
            183895.50908333,
            194696.07972222,
            204305.66308333,
            212105.24738889,
            217677.25102778,
            220830.11786111,
            221589.25106944,
            220150.69498611,
            216828.44369444,
            212005.34018056,
            206099.23534722,
            199534.26356944,
            192720.28676389,
            186032.24869444,
            179786.74659722,
            174219.32091667,
            169471.61906944,
            165602.61593056,
            162594.35872222,
            160382.21744444,
            158874.45340278,
            157958.32315278,
            157604.38413889,
            157604.38413889,
            157604.38413889,
            157604.38413889,
        ]
    )

    demographic_events = []
    population_configurations = [
        msprime.PopulationConfiguration(
            # NCW
            initial_size=sizes_NCW[0],
            metadata=populations[0].asdict(),
        )
    ]

    for i, t in enumerate(times_NCW):
        curr_time = t
        demographic_events.append(
            msprime.PopulationParametersChange(
                time=curr_time, initial_size=sizes_NCW[i + 1], population_id=1
            )
        )
    citations = [
        _ZhangEtAl,
    ]
    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=3,
        mutation_rate=mutation_rate,
        recombination_rate=recombination_rate,
        population_configurations=population_configurations,
        demographic_events=demographic_events,
    )


_species.add_demographic_model(_WildBoar_4Z22_NCW())


def _WildBoar_4Z22_EUW():
    id = "WildBoar_4Z22_EUW"
    description = "Piecewise size model for European wild boar (Zhang et al. 2022)"
    long_description = """
    This demographic model is a piecewise size model for European wild boar
    from Zhang et al. (2022).
    """
    populations = [
        stdpopsim.Population(id="EUW", description="European wild boar"),
    ]

    # Recombination rate used in simulations by Zhang et al. (2022).
    recombination_rate = 1.0e-8

    # Mutation rate as reported by Zhang et al. (2022, p. 1042).
    mutation_rate = 3.6e-9

    # Arrays of time intervals and effective population sizes (Ne)
    # extracted from Figure 3A of Zhang et al. (2022).

    times_EUW = np.array(
        [
            4.73944067e02,
            9.96184317e02,
            1.57169447e03,
            2.20586658e03,
            2.90469702e03,
            3.67478640e03,
            4.52343262e03,
            5.45853783e03,
            6.48902685e03,
            7.62456820e03,
            8.87589953e03,
            1.02547811e04,
            1.17742748e04,
            1.34486974e04,
            1.52938534e04,
            1.73270808e04,
            1.95676705e04,
            2.20366792e04,
            2.47574417e04,
            2.77555702e04,
            3.10594196e04,
            3.47000872e04,
            3.87120172e04,
            4.31329076e04,
            4.80045934e04,
            5.33730000e04,
            5.92887014e04,
            6.58076171e04,
            7.29911514e04,
            8.09071236e04,
            8.96301860e04,
            9.92426139e04,
            1.09835096e05,
            1.21507619e05,
            1.34370208e05,
            1.48544232e05,
            1.64163422e05,
            1.81375131e05,
            2.00341726e05,
            2.21242073e05,
            2.44273403e05,
            2.69652931e05,
            2.97620093e05,
            3.28438776e05,
            3.62399639e05,
            3.99823091e05,
            4.41062175e05,
            4.86505909e05,
            5.36582961e05,
            5.91765790e05,
            6.52574915e05,
            7.19584033e05,
            7.93425272e05,
            8.74795183e05,
            9.64461440e05,
            1.06326985e06,
            1.17215260e06,
        ]
    )
    sizes_EUW = np.array(
        [
            5547.85555833,
            5547.85555833,
            3441.95138333,
            3441.95138333,
            6541.99460833,
            6554.498625,
            7333.81495,
            8576.63983333,
            10218.268475,
            12347.62349167,
            15571.03348333,
            19965.80023333,
            24556.33154167,
            27640.593675,
            27863.434775,
            25457.68985833,
            21804.402,
            18286.82223333,
            15815.00125833,
            14837.34055,
            15508.11829167,
            17894.224,
            22151.841675,
            28668.64231667,
            38124.11929167,
            51412.077375,
            69317.2482,
            91987.3558,
            118313.02894167,
            145672.11955,
            170202.141525,
            187892.04803333,
            196750.74873333,
            198505.42484167,
            197046.38273333,
            195519.707375,
            195313.25165,
            196421.64673333,
            198227.826375,
            200394.78640833,
            203122.56785833,
            206885.76555833,
            212188.14263333,
            219431.71019167,
            228687.02993333,
            239404.45968333,
            250383.241975,
            260037.1563,
            266737.61259167,
            269079.62885833,
            266139.37209167,
            257709.20124167,
            244366.97386667,
            227322.813825,
            185569.531525,
            185569.531525,
            185569.531525,
            185569.531525,
        ]
    )

    demographic_events = []
    population_configurations = [
        msprime.PopulationConfiguration(
            # EUW
            initial_size=sizes_EUW[0],
            metadata=populations[0].asdict(),
        )
    ]

    for i, t in enumerate(times_EUW):
        curr_time = t
        demographic_events.append(
            msprime.PopulationParametersChange(
                time=curr_time, initial_size=sizes_EUW[i + 1], population_id=0
            )
        )
    citations = [
        _ZhangEtAl,
    ]
    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=3,
        mutation_rate=mutation_rate,
        recombination_rate=recombination_rate,
        population_configurations=population_configurations,
        demographic_events=demographic_events,
    )


_species.add_demographic_model(_WildBoar_4Z22_EUW())

_WangEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1016/j.xgen.2025.100954",
    year=2025,
    author="Wang et al.",
    reasons={stdpopsim.CiteReason.DEM_MODEL},
    # Figure 2b Demographic history of four geographical populations using MSMC.
)


def _WildBoar_6W25():
    id = "WildBoar_6W25"
    description = "The demographic model for wild boar (Wang et al. 2025)"
    long_description = """
    This demographic model of wild boar with Sus Cebifrons (outgroup) and 5 wild boar
    populations including Central Asian wild boars (CAW), European wild boars (EUW),
    Southern Chinese wild boar (NAW), Northeast Asian wild boars (SCW),
    and Near East wild boars (NEW), from Wang et al. (2025).
    SVDquartets estimated the branching pattern among the five populations
    (SAW, SCW, NAW, CAW, and EUW).
    Fastsimcoal2 analysed more recent demographic fluctuations and respective
    divergence times based on the species tree.
    Model parameters were a generation interval of 6.5 years (Figure S11) and
    a mutation rate of 2.5e-8 (page e2).
    Estimated effective population sizes, divergence times,
    and migration rates are given in Figure S11.
    """
    # Mutation rate used by Wang et al. (2025).
    mutation_rate = 2.5e-8
    # generation time declaimed in Figure S11 by Wang et al. (2025).
    generation_time = 6.5

    # Estimated effective population sizes. from Figure S11D
    N_CAW = 87338
    N_EUW = 87275
    N_NAW = 85535
    N_SCW = 88237
    N_OUT = 1158
    N_NEW = 81722
    # Estimated divergence times, from Figure S11B
    T1 = 3_662_321  # (OUT,(SCW,NAW,CAW,EUW,NEW))
    T2 = 1_826_091  # ((SCW,NAW),(CAW,EUW,NEW))
    T3 = 1_064_310  # (SCW,NAW)
    T4 = 948_246  # (CAW,(EUW,NEW))
    T5 = 682571.5  # (CAW,EUW)

    # Migration rates, from Figure S11B
    # migration is bidirectional.
    m_OUT_SCW = 1.29e-2
    m_SCW_NAW = 4e-4
    m_NAW_CAW = 1.01e-1
    m_CAW_EUW = 5.8e-3
    m_CAW_NEW = 1e-4
    m_EUW_NEW = 2.48e-2

    model = msprime.Demography()
    model.add_population(
        initial_size=N_CAW,
        name="CAW",
        description="Central Asian wild boars",
    )
    model.add_population(
        initial_size=N_EUW,
        name="EUW",
        description="European wild boars",
    )
    model.add_population(
        initial_size=N_NAW,
        name="NAW",
        description="Northeast Asian wild boars",
    )
    model.add_population(
        initial_size=N_SCW,
        name="SCW",
        description="Southern Chinese wild boar",
    )
    model.add_population(
        initial_size=N_OUT,
        name="OUT",
        description="Sus Cebifrons",
    )
    model.add_population(
        initial_size=N_NEW,
        name="NEW",
        description="Near East wild boars",
    )

    model.set_migration_rate(rate=m_OUT_SCW, source="OUT", dest="SCW")
    model.set_migration_rate(rate=m_OUT_SCW, source="SCW", dest="OUT")
    model.set_migration_rate(rate=m_SCW_NAW, source="SCW", dest="NAW")
    model.set_migration_rate(rate=m_SCW_NAW, source="NAW", dest="SCW")
    model.set_migration_rate(rate=m_NAW_CAW, source="NAW", dest="CAW")
    model.set_migration_rate(rate=m_NAW_CAW, source="CAW", dest="NAW")
    model.set_migration_rate(rate=m_CAW_EUW, source="CAW", dest="EUW")
    model.set_migration_rate(rate=m_CAW_EUW, source="EUW", dest="CAW")
    model.set_migration_rate(rate=m_CAW_NEW, source="CAW", dest="NEW")
    model.set_migration_rate(rate=m_CAW_NEW, source="NEW", dest="CAW")
    model.set_migration_rate(rate=m_EUW_NEW, source="EUW", dest="NEW")
    model.set_migration_rate(rate=m_EUW_NEW, source="NEW", dest="EUW")

    # model.add_population_split(
    #     time=T1, derived=["OUT", "ancWB"], ancestral="root"
    # )
    # model.add_population_split(
    #     time=T2, derived=["ancSCW_NAW", "ancEUW_NEW_CAW"], ancestral="ancWB"
    # )
    # model.add_population_split(
    #     time=T3, derived=["SCW", "NAW"], ancestral="ancSCW_NAW"
    # )
    # model.add_population_split(
    #     time=T4, derived=["CAW", "ancEUW_NEW"], ancestral="ancEUW_NEW_CAW"
    # )
    # model.add_population_split(
    #     time=T5, derived=["EUW", "NEW"], ancestral="ancEUW_NEW"
    # )

    # Modeling Strategy: Budding vs. Bifurcating
    # Although Figures 2C and S11A depict a standard tree structure, Figure 2D and the
    # main text (Wang et al. 2025) suggest a "budding" or "successive divergence" model
    # where lineages emerge from a central ancestral trunk.

    # Assumption on Ancestral Ne:
    # Since the authors did not report ancestral population
    # for internal nodes, we treat the 'ancestral' labels in the split events as
    # persistent lineages. This approach assumes the ancestral trunk maintains
    # the Ne of the source population (e.g., OUT or CAW) until the next
    # divergence event.

    # From figure 2D and main text (Wang et al. 2025), Asian wild boars and Southeast
    # Asian Suids split ∼3.6 million years ago (mya),
    # with Central Asian and Southern Chinese ancestors diverging ∼1.8 mya.
    # The split between Central Asian and European-Near East ancestors occurred ∼0.9 mya,
    # followed by a European-Near East divergence ∼0.6 mya.

    model.add_population_split(time=T1, derived=["SCW"], ancestral="OUT")
    model.add_population_split(time=T2, derived=["CAW"], ancestral="SCW")
    model.add_population_split(time=T3, derived=["NAW"], ancestral="SCW")
    model.add_population_split(time=T4, derived=["EUW"], ancestral="CAW")
    model.add_population_split(time=T5, derived=["NEW"], ancestral="CAW")

    citations = [
        _WangEtAl,
    ]
    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        citations=citations,
        generation_time=generation_time,
        mutation_rate=mutation_rate,
        model=model,
    )


_species.add_demographic_model(_WildBoar_6W25())
