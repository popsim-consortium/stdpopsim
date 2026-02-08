import msprime
import stdpopsim
from . import species

_species = stdpopsim.get_species("SusScr")


def _WildBoar_4Z22():
    id = "WildBoar_4Z22"
    description = "Piecewise size model for wild boar (pig) (Zhang et al. 2022)"
    long_description = """
    This demographic model is a piecewise size model for wild boar
    from Zhang et al. (2022).
    """
    populations = [
        stdpopsim.Population(id="EUW", description="European wild boar"),
        stdpopsim.Population(id="NCW", description="North Chinese wild boar"),
        stdpopsim.Population(id="SCW", description="South Chinese wild boar"),
        stdpopsim.Population(id="SMW", description="Sumatran wild boar"),
    ]
    citations = [species._ZhangEtAl]
    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=3,
        mutation_rate=3.6e-9,
        population_configurations=[
            stdpopsim.PopulationConfiguration(
                # EUW,0.0 generations at 5547.855558333333, [0.0 - 473.94406666666663)
                initial_size=5547.855558333333,
                metadata=populations[0].asdict(),
            ),
            stdpopsim.PopulationConfiguration(
                # NCW,0.0 generations at 96756.31630555558, [0.0 - 1650.3988611111113)
                initial_size=96756.31630555558,
                metadata=populations[1].asdict(),
            ),
            stdpopsim.PopulationConfiguration(
                # SCW,0.0 generations at 33850.158784722225, [0.0 - 1994.0369444444445)
                initial_size=33850.158784722225,
                metadata=populations[2].asdict(),
            ),
            stdpopsim.PopulationConfiguration(
                # SMW,0.0 generations at 1779489.4096979166, [0.0 - 273.92735833333325)
                initial_size=1779489.4096979166,
                metadata=populations[3].asdict(),
            ),
        ],
        demographic_events=[
            # EUW, 473.94406666666663 generations at 5547.855558333333,
            # [473.94406666666663 - 996.1843166666665)
            msprime.PopulationParametersChange(
                time=473.94406666666663, initial_size=5547.855558333333, population_id=0
            ),
            # EUW, 996.1843166666665 generations at 3441.951383333333,
            # [996.1843166666665 - 1571.6944666666666)
            msprime.PopulationParametersChange(
                time=996.1843166666665, initial_size=3441.951383333333, population_id=0
            ),
            # EUW, 1571.6944666666666 generations at 3441.951383333333,
            # [1571.6944666666666 - 2205.866583333333)
            msprime.PopulationParametersChange(
                time=1571.6944666666666, initial_size=3441.951383333333, population_id=0
            ),
            # EUW, 2205.866583333333 generations at 6541.9946083333325,
            # [2205.866583333333 - 2904.6970166666665)
            msprime.PopulationParametersChange(
                time=2205.866583333333, initial_size=6541.9946083333325, population_id=0
            ),
            # EUW, 2904.6970166666665 generations at 6554.498624999999,
            # [2904.6970166666665 - 3674.7863999999995)
            msprime.PopulationParametersChange(
                time=2904.6970166666665, initial_size=6554.498624999999, population_id=0
            ),
            # EUW, 3674.7863999999995 generations at 7333.814949999999,
            # [3674.7863999999995 - 4523.432616666666)
            msprime.PopulationParametersChange(
                time=3674.7863999999995, initial_size=7333.814949999999, population_id=0
            ),
            # EUW, 4523.432616666666 generations at 8576.639833333333,
            # [4523.432616666666 - 5458.537833333333)
            msprime.PopulationParametersChange(
                time=4523.432616666666, initial_size=8576.639833333333, population_id=0
            ),
            # EUW, 5458.537833333333 generations at 10218.268474999999,
            # [5458.537833333333 - 6489.026849999999)
            msprime.PopulationParametersChange(
                time=5458.537833333333, initial_size=10218.268474999999, population_id=0
            ),
            # EUW, 6489.026849999999 generations at 12347.623491666667,
            # [6489.026849999999 - 7624.5682)
            msprime.PopulationParametersChange(
                time=6489.026849999999, initial_size=12347.623491666667, population_id=0
            ),
            # EUW, 7624.5682 generations at 15571.03348333333,
            # [7624.5682 - 8875.899533333333)
            msprime.PopulationParametersChange(
                time=7624.5682, initial_size=15571.03348333333, population_id=0
            ),
            # EUW, 8875.899533333333 generations at 19965.80023333333,
            # [8875.899533333333 - 10254.781133333332)
            msprime.PopulationParametersChange(
                time=8875.899533333333, initial_size=19965.80023333333, population_id=0
            ),
            # EUW, 10254.781133333332 generations at 24556.331541666663,
            # [10254.781133333332 - 11774.274816666666)
            msprime.PopulationParametersChange(
                time=10254.781133333332,
                initial_size=24556.331541666663,
                population_id=0,
            ),
            # EUW, 11774.274816666666 generations at 27640.593674999996,
            # [11774.274816666666 - 13448.697449999998)
            msprime.PopulationParametersChange(
                time=11774.274816666666,
                initial_size=27640.593674999996,
                population_id=0,
            ),
            # EUW, 13448.697449999998 generations at 27863.434774999998,
            # [13448.697449999998 - 15293.853366666664)
            msprime.PopulationParametersChange(
                time=13448.697449999998,
                initial_size=27863.434774999998,
                population_id=0,
            ),
            # EUW, 15293.853366666664 generations at 25457.68985833333,
            # [15293.853366666664 - 17327.08085)
            msprime.PopulationParametersChange(
                time=15293.853366666664, initial_size=25457.68985833333, population_id=0
            ),
            # EUW, 17327.08085 generations at 21804.402,
            # [17327.08085 - 19567.67048333333)
            msprime.PopulationParametersChange(
                time=17327.08085, initial_size=21804.402, population_id=0
            ),
            # EUW, 19567.67048333333 generations at 18286.82223333333,
            # [19567.67048333333 - 22036.679216666667)
            msprime.PopulationParametersChange(
                time=19567.67048333333, initial_size=18286.82223333333, population_id=0
            ),
            # EUW, 22036.679216666667 generations at 15815.001258333332,
            # [22036.679216666667 - 24757.44168333333)
            msprime.PopulationParametersChange(
                time=22036.679216666667,
                initial_size=15815.001258333332,
                population_id=0,
            ),
            # EUW, 24757.44168333333 generations at 14837.340549999999,
            # [24757.44168333333 - 27755.570199999995)
            msprime.PopulationParametersChange(
                time=24757.44168333333, initial_size=14837.340549999999, population_id=0
            ),
            # EUW, 27755.570199999995 generations at 15508.118291666666,
            # [27755.570199999995 - 31059.419599999997)
            msprime.PopulationParametersChange(
                time=27755.570199999995,
                initial_size=15508.118291666666,
                population_id=0,
            ),
            # EUW, 31059.419599999997 generations at 17894.224,
            # [31059.419599999997 - 34700.08723333333)
            msprime.PopulationParametersChange(
                time=31059.419599999997, initial_size=17894.224, population_id=0
            ),
            # EUW, 34700.08723333333 generations at 22151.841674999996,
            # [34700.08723333333 - 38712.01725)
            msprime.PopulationParametersChange(
                time=34700.08723333333, initial_size=22151.841674999996, population_id=0
            ),
            # EUW, 38712.01725 generations at 28668.642316666665,
            # [38712.01725 - 43132.90763333333)
            msprime.PopulationParametersChange(
                time=38712.01725, initial_size=28668.642316666665, population_id=0
            ),
            # EUW, 43132.90763333333 generations at 38124.11929166666,
            # [43132.90763333333 - 48004.593383333326)
            msprime.PopulationParametersChange(
                time=43132.90763333333, initial_size=38124.11929166666, population_id=0
            ),
            # EUW, 48004.593383333326 generations at 51412.07737499999,
            # [48004.593383333326 - 53373.00003333332)
            msprime.PopulationParametersChange(
                time=48004.593383333326, initial_size=51412.07737499999, population_id=0
            ),
            # EUW, 53373.00003333332 generations at 69317.24819999999,
            # [53373.00003333332 - 59288.70144999999)
            msprime.PopulationParametersChange(
                time=53373.00003333332, initial_size=69317.24819999999, population_id=0
            ),
            # EUW, 59288.70144999999 generations at 91987.35579999999,
            # [59288.70144999999 - 65807.61708333332)
            msprime.PopulationParametersChange(
                time=59288.70144999999, initial_size=91987.35579999999, population_id=0
            ),
            # EUW, 65807.61708333332 generations at 118313.02894166666,
            # [65807.61708333332 - 72991.15141666666)
            msprime.PopulationParametersChange(
                time=65807.61708333332, initial_size=118313.02894166666, population_id=0
            ),
            # EUW, 72991.15141666666 generations at 145672.11954999997,
            # [72991.15141666666 - 80907.12363333332)
            msprime.PopulationParametersChange(
                time=72991.15141666666, initial_size=145672.11954999997, population_id=0
            ),
            # EUW, 80907.12363333332 generations at 170202.14152499998,
            # [80907.12363333332 - 89630.18596666666)
            msprime.PopulationParametersChange(
                time=80907.12363333332, initial_size=170202.14152499998, population_id=0
            ),
            # EUW, 89630.18596666666 generations at 187892.0480333333,
            # [89630.18596666666 - 99242.61391666665)
            msprime.PopulationParametersChange(
                time=89630.18596666666, initial_size=187892.0480333333, population_id=0
            ),
            # EUW, 99242.61391666665 generations at 196750.7487333333,
            # [99242.61391666665 - 109835.09646666666)
            msprime.PopulationParametersChange(
                time=99242.61391666665, initial_size=196750.7487333333, population_id=0
            ),
            # EUW, 109835.09646666666 generations at 198505.42484166665,
            # [109835.09646666666 - 121507.61926666665)
            msprime.PopulationParametersChange(
                time=109835.09646666666,
                initial_size=198505.42484166665,
                population_id=0,
            ),
            # EUW, 121507.61926666665 generations at 197046.3827333333,
            # [121507.61926666665 - 134370.20836666666)
            msprime.PopulationParametersChange(
                time=121507.61926666665, initial_size=197046.3827333333, population_id=0
            ),
            # EUW, 134370.20836666666 generations at 195519.70737499997,
            # [134370.20836666666 - 148544.23174999998)
            msprime.PopulationParametersChange(
                time=134370.20836666666,
                initial_size=195519.70737499997,
                population_id=0,
            ),
            # EUW, 148544.23174999998 generations at 195313.25165,
            # [148544.23174999998 - 164163.42196666665)
            msprime.PopulationParametersChange(
                time=148544.23174999998, initial_size=195313.25165, population_id=0
            ),
            # EUW, 164163.42196666665 generations at 196421.6467333333,
            # [164163.42196666665 - 181375.13118333332)
            msprime.PopulationParametersChange(
                time=164163.42196666665, initial_size=196421.6467333333, population_id=0
            ),
            # EUW, 181375.13118333332 generations at 198227.826375,
            # [181375.13118333332 - 200341.7256833333)
            msprime.PopulationParametersChange(
                time=181375.13118333332, initial_size=198227.826375, population_id=0
            ),
            # EUW, 200341.7256833333 generations at 200394.7864083333,
            # [200341.7256833333 - 221242.0733333333)
            msprime.PopulationParametersChange(
                time=200341.7256833333, initial_size=200394.7864083333, population_id=0
            ),
            # EUW, 221242.0733333333 generations at 203122.56785833332,
            # [221242.0733333333 - 244273.40291666662)
            msprime.PopulationParametersChange(
                time=221242.0733333333, initial_size=203122.56785833332, population_id=0
            ),
            # EUW, 244273.40291666662 generations at 206885.7655583333,
            # [244273.40291666662 - 269652.93104999996)
            msprime.PopulationParametersChange(
                time=244273.40291666662, initial_size=206885.7655583333, population_id=0
            ),
            # EUW, 269652.93104999996 generations at 212188.14263333334,
            # [269652.93104999996 - 297620.0933833333)
            msprime.PopulationParametersChange(
                time=269652.93104999996,
                initial_size=212188.14263333334,
                population_id=0,
            ),
            # EUW, 297620.0933833333 generations at 219431.71019166664,
            # [297620.0933833333 - 328438.77579999994)
            msprime.PopulationParametersChange(
                time=297620.0933833333, initial_size=219431.71019166664, population_id=0
            ),
            # EUW, 328438.77579999994 generations at 228687.0299333333,
            # [328438.77579999994 - 362399.6385833333)
            msprime.PopulationParametersChange(
                time=328438.77579999994, initial_size=228687.0299333333, population_id=0
            ),
            # EUW, 362399.6385833333 generations at 239404.45968333332,
            # [362399.6385833333 - 399823.09134999994)
            msprime.PopulationParametersChange(
                time=362399.6385833333, initial_size=239404.45968333332, population_id=0
            ),
            # EUW, 399823.09134999994 generations at 250383.24197499998,
            # [399823.09134999994 - 441062.17501666665)
            msprime.PopulationParametersChange(
                time=399823.09134999994,
                initial_size=250383.24197499998,
                population_id=0,
            ),
            # EUW, 441062.17501666665 generations at 260037.15629999997,
            # [441062.17501666665 - 486505.9085999999)
            msprime.PopulationParametersChange(
                time=441062.17501666665,
                initial_size=260037.15629999997,
                population_id=0,
            ),
            # EUW, 486505.9085999999 generations at 266737.61259166664,
            # [486505.9085999999 - 536582.9614)
            msprime.PopulationParametersChange(
                time=486505.9085999999, initial_size=266737.61259166664, population_id=0
            ),
            # EUW, 536582.9614 generations at 269079.62885833334,
            # [536582.9614 - 591765.7900166665)
            msprime.PopulationParametersChange(
                time=536582.9614, initial_size=269079.62885833334, population_id=0
            ),
            # EUW, 591765.7900166665 generations at 266139.3720916666,
            # [591765.7900166665 - 652574.9148166665)
            msprime.PopulationParametersChange(
                time=591765.7900166665, initial_size=266139.3720916666, population_id=0
            ),
            # EUW, 652574.9148166665 generations at 257709.20124166663,
            # [652574.9148166665 - 719584.0330999999)
            msprime.PopulationParametersChange(
                time=652574.9148166665, initial_size=257709.20124166663, population_id=0
            ),
            # EUW, 719584.0330999999 generations at 244366.97386666667,
            # [719584.0330999999 - 793425.2717166665)
            msprime.PopulationParametersChange(
                time=719584.0330999999, initial_size=244366.97386666667, population_id=0
            ),
            # EUW, 793425.2717166665 generations at 227322.81382499996,
            # [793425.2717166665 - 874795.1834166667)
            msprime.PopulationParametersChange(
                time=793425.2717166665, initial_size=227322.81382499996, population_id=0
            ),
            # EUW, 874795.1834166667 generations at 185569.53152499997,
            # [874795.1834166667 - 964461.44045)
            msprime.PopulationParametersChange(
                time=874795.1834166667, initial_size=185569.53152499997, population_id=0
            ),
            # EUW, 964461.44045 generations at 185569.53152499997,
            # [964461.44045 - 1063269.8535499999)
            msprime.PopulationParametersChange(
                time=964461.44045, initial_size=185569.53152499997, population_id=0
            ),
            # EUW, 1063269.8535499999 generations at 185569.53152499997,
            # [1063269.8535499999 - 1172152.5994833333)
            msprime.PopulationParametersChange(
                time=1063269.8535499999,
                initial_size=185569.53152499997,
                population_id=0,
            ),
            # EUW, 1172152.5994833333 generations at 185569.53152499997,
            # [1172152.5994833333 - inf)
            msprime.PopulationParametersChange(
                time=1172152.5994833333,
                initial_size=185569.53152499997,
                population_id=0,
            ),
            # NCW, 1650.3988611111113 generations at 96756.31630555558,
            # [1650.3988611111113 - 3452.4986944444454)
            msprime.PopulationParametersChange(
                time=1650.3988611111113, initial_size=96756.31630555558, population_id=1
            ),
            # NCW, 3452.4986944444454 generations at 23458.804958333338,
            # [3452.4986944444454 - 5420.12313888889)
            msprime.PopulationParametersChange(
                time=3452.4986944444454,
                initial_size=23458.804958333338,
                population_id=1,
            ),
            # NCW, 5420.12313888889 generations at 23458.804958333338,
            # [5420.12313888889 - 7568.891111111113)
            msprime.PopulationParametersChange(
                time=5420.12313888889, initial_size=23458.804958333338, population_id=1
            ),
            # NCW, 7568.891111111113 generations at 29425.859472222226,
            # [7568.891111111113 - 9914.960111111113)
            msprime.PopulationParametersChange(
                time=7568.891111111113, initial_size=29425.859472222226, population_id=1
            ),
            # NCW, 9914.960111111113 generations at 33474.83920833334,
            # [9914.960111111113 - 12476.821500000004)
            msprime.PopulationParametersChange(
                time=9914.960111111113, initial_size=33474.83920833334, population_id=1
            ),
            # NCW, 12476.821500000004 generations at 37881.70756944445,
            # [12476.821500000004 - 15274.223333333337)
            msprime.PopulationParametersChange(
                time=12476.821500000004, initial_size=37881.70756944445, population_id=1
            ),
            # NCW, 15274.223333333337 generations at 44466.24787500001,
            # [15274.223333333337 - 18328.70894444445)
            msprime.PopulationParametersChange(
                time=15274.223333333337, initial_size=44466.24787500001, population_id=1
            ),
            # NCW, 18328.70894444445 generations at 54119.27695833334,
            # [18328.70894444445 - 21663.976000000002)
            msprime.PopulationParametersChange(
                time=18328.70894444445, initial_size=54119.27695833334, population_id=1
            ),
            # NCW, 21663.976000000002 generations at 66767.63725000001,
            # [21663.976000000002 - 25305.876500000006)
            msprime.PopulationParametersChange(
                time=21663.976000000002, initial_size=66767.63725000001, population_id=1
            ),
            # NCW, 25305.876500000006 generations at 82821.19066666669,
            # [25305.876500000006 - 29282.596305555562)
            msprime.PopulationParametersChange(
                time=25305.876500000006, initial_size=82821.19066666669, population_id=1
            ),
            # NCW, 29282.596305555562 generations at 102492.67762500001,
            # [29282.596305555562 - 33624.83466666667)
            msprime.PopulationParametersChange(
                time=29282.596305555562,
                initial_size=102492.67762500001,
                population_id=1,
            ),
            # NCW, 33624.83466666667 generations at 124769.92051388892,
            # [33624.83466666667 - 38366.34280555556)
            msprime.PopulationParametersChange(
                time=33624.83466666667, initial_size=124769.92051388892, population_id=1
            ),
            # NCW, 38366.34280555556 generations at 147268.88020833337,
            # [38366.34280555556 - 43543.744388888896)
            msprime.PopulationParametersChange(
                time=38366.34280555556, initial_size=147268.88020833337, population_id=1
            ),
            # NCW, 43543.744388888896 generations at 167371.8621805556,
            # [43543.744388888896 - 49196.894583333335)
            msprime.PopulationParametersChange(
                time=43543.744388888896, initial_size=167371.8621805556, population_id=1
            ),
            # NCW, 49196.894583333335 generations at 183211.86730555558,
            # [49196.894583333335 - 55369.957222222234)
            msprime.PopulationParametersChange(
                time=49196.894583333335,
                initial_size=183211.86730555558,
                population_id=1,
            ),
            # NCW, 55369.957222222234 generations at 193154.9235138889,
            # [55369.957222222234 - 62110.327638888906)
            msprime.PopulationParametersChange(
                time=55369.957222222234, initial_size=193154.9235138889, population_id=1
            ),
            # NCW, 62110.327638888906 generations at 195884.55361111116,
            # [62110.327638888906 - 69470.42794444445)
            msprime.PopulationParametersChange(
                time=62110.327638888906,
                initial_size=195884.55361111116,
                population_id=1,
            ),
            # NCW, 69470.42794444445 generations at 191547.07273611115,
            # [69470.42794444445 - 77507.16844444445)
            msprime.PopulationParametersChange(
                time=69470.42794444445, initial_size=191547.07273611115, population_id=1
            ),
            # NCW, 77507.16844444445 generations at 181982.0122638889,
            # [77507.16844444445 - 86282.66575000001)
            msprime.PopulationParametersChange(
                time=77507.16844444445, initial_size=181982.0122638889, population_id=1
            ),
            # NCW, 86282.66575000001 generations at 169846.56283333336,
            # [86282.66575000001 - 95864.78136111112)
            msprime.PopulationParametersChange(
                time=86282.66575000001, initial_size=169846.56283333336, population_id=1
            ),
            # NCW, 95864.78136111112 generations at 157649.26608333338,
            # [95864.78136111112 - 106327.8397777778)
            msprime.PopulationParametersChange(
                time=95864.78136111112, initial_size=157649.26608333338, population_id=1
            ),
            # NCW, 106327.8397777778 generations at 147192.93995833336,
            # [106327.8397777778 - 117752.8080277778)
            msprime.PopulationParametersChange(
                time=106327.8397777778, initial_size=147192.93995833336, population_id=1
            ),
            # NCW, 117752.8080277778 generations at 139407.9871666667,
            # [117752.8080277778 - 130228.01377777781)
            msprime.PopulationParametersChange(
                time=117752.8080277778, initial_size=139407.9871666667, population_id=1
            ),
            # NCW, 130228.01377777781 generations at 134604.1828888889,
            # [130228.01377777781 - 143850.04297222226)
            msprime.PopulationParametersChange(
                time=130228.01377777781, initial_size=134604.1828888889, population_id=1
            ),
            # NCW, 143850.04297222226 generations at 132773.80743055558,
            # [143850.04297222226 - 158724.2784166667)
            msprime.PopulationParametersChange(
                time=143850.04297222226,
                initial_size=132773.80743055558,
                population_id=1,
            ),
            # NCW, 158724.2784166667 generations at 133793.61497222225,
            # [158724.2784166667 - 174965.97694444447)
            msprime.PopulationParametersChange(
                time=158724.2784166667, initial_size=133793.61497222225, population_id=1
            ),
            # NCW, 174965.97694444447 generations at 137496.10609722225,
            # [174965.97694444447 - 192700.62847222225)
            msprime.PopulationParametersChange(
                time=174965.97694444447,
                initial_size=137496.10609722225,
                population_id=1,
            ),
            # NCW, 192700.62847222225 generations at 143635.95609722225,
            # [192700.62847222225 - 212065.57175000003)
            msprime.PopulationParametersChange(
                time=192700.62847222225,
                initial_size=143635.95609722225,
                population_id=1,
            ),
            # NCW, 212065.57175000003 generations at 151873.229125,
            # [212065.57175000003 - 233210.89200000002)
            msprime.PopulationParametersChange(
                time=212065.57175000003, initial_size=151873.229125, population_id=1
            ),
            # NCW, 233210.89200000002 generations at 161743.12776388892,
            # [233210.89200000002 - 256299.95950000006)
            msprime.PopulationParametersChange(
                time=233210.89200000002,
                initial_size=161743.12776388892,
                population_id=1,
            ),
            # NCW, 256299.95950000006 generations at 172654.01822222225,
            # [256299.95950000006 - 281511.58391666674)
            msprime.PopulationParametersChange(
                time=256299.95950000006,
                initial_size=172654.01822222225,
                population_id=1,
            ),
            # NCW, 281511.58391666674 generations at 183895.50908333337,
            # [281511.58391666674 - 309040.9119444445)
            msprime.PopulationParametersChange(
                time=281511.58391666674,
                initial_size=183895.50908333337,
                population_id=1,
            ),
            # NCW, 309040.9119444445 generations at 194696.07972222226,
            # [309040.9119444445 - 339100.86352777784)
            msprime.PopulationParametersChange(
                time=309040.9119444445, initial_size=194696.07972222226, population_id=1
            ),
            # NCW, 339100.86352777784 generations at 204305.66308333338,
            # [339100.86352777784 - 371924.2861944445)
            msprime.PopulationParametersChange(
                time=339100.86352777784,
                initial_size=204305.66308333338,
                population_id=1,
            ),
            # NCW, 371924.2861944445 generations at 212105.24738888894,
            # [371924.2861944445 - 407765.03222222225)
            msprime.PopulationParametersChange(
                time=371924.2861944445, initial_size=212105.24738888894, population_id=1
            ),
            # NCW, 407765.03222222225 generations at 217677.25102777782,
            # [407765.03222222225 - 446900.65155555564)
            msprime.PopulationParametersChange(
                time=407765.03222222225,
                initial_size=217677.25102777782,
                population_id=1,
            ),
            # NCW, 446900.65155555564 generations at 220830.11786111115,
            # [446900.65155555564 - 489633.82802777784)
            msprime.PopulationParametersChange(
                time=446900.65155555564,
                initial_size=220830.11786111115,
                population_id=1,
            ),
            # NCW, 489633.82802777784 generations at 221589.2510694445,
            # [489633.82802777784 - 536295.4313333334)
            msprime.PopulationParametersChange(
                time=489633.82802777784, initial_size=221589.2510694445, population_id=1
            ),
            # NCW, 536295.4313333334 generations at 220150.69498611114,
            # [536295.4313333334 - 587246.6713611112)
            msprime.PopulationParametersChange(
                time=536295.4313333334, initial_size=220150.69498611114, population_id=1
            ),
            # NCW, 587246.6713611112 generations at 216828.44369444446,
            # [587246.6713611112 - 642881.6115833335)
            msprime.PopulationParametersChange(
                time=587246.6713611112, initial_size=216828.44369444446, population_id=1
            ),
            # NCW, 642881.6115833335 generations at 212005.3401805556,
            # [642881.6115833335 - 703631.1186666668)
            msprime.PopulationParametersChange(
                time=642881.6115833335, initial_size=212005.3401805556, population_id=1
            ),
            # NCW, 703631.1186666668 generations at 206099.23534722227,
            # [703631.1186666668 - 769965.1963333335)
            msprime.PopulationParametersChange(
                time=703631.1186666668, initial_size=206099.23534722227, population_id=1
            ),
            # NCW, 769965.1963333335 generations at 199534.26356944448,
            # [769965.1963333335 - 842397.2940277779)
            msprime.PopulationParametersChange(
                time=769965.1963333335, initial_size=199534.26356944448, population_id=1
            ),
            # NCW, 842397.2940277779 generations at 192720.28676388893,
            # [842397.2940277779 - 921487.8974722225)
            msprime.PopulationParametersChange(
                time=842397.2940277779, initial_size=192720.28676388893, population_id=1
            ),
            # NCW, 921487.8974722225 generations at 186032.24869444445,
            # [921487.8974722225 - 1007849.196388889)
            msprime.PopulationParametersChange(
                time=921487.8974722225, initial_size=186032.24869444445, population_id=1
            ),
            # NCW, 1007849.196388889 generations at 179786.74659722226,
            # [1007849.196388889 - 1102149.5726944448)
            msprime.PopulationParametersChange(
                time=1007849.196388889, initial_size=179786.74659722226, population_id=1
            ),
            # NCW, 1102149.5726944448 generations at 174219.3209166667,
            # [1102149.5726944448 - 1205118.8068055557)
            msprime.PopulationParametersChange(
                time=1102149.5726944448, initial_size=174219.3209166667, population_id=1
            ),
            # NCW, 1205118.8068055557 generations at 169471.61906944448,
            # [1205118.8068055557 - 1317553.822527778)
            msprime.PopulationParametersChange(
                time=1205118.8068055557,
                initial_size=169471.61906944448,
                population_id=1,
            ),
            # NCW, 1317553.822527778 generations at 165602.6159305556,
            # [1317553.822527778 - 1440324.7910000002)
            msprime.PopulationParametersChange(
                time=1317553.822527778, initial_size=165602.6159305556, population_id=1
            ),
            # NCW, 1440324.7910000002 generations at 162594.35872222224,
            # [1440324.7910000002 - 1574381.9527500002)
            msprime.PopulationParametersChange(
                time=1440324.7910000002,
                initial_size=162594.35872222224,
                population_id=1,
            ),
            # NCW, 1574381.9527500002 generations at 160382.21744444448,
            # [1574381.9527500002 - 1720762.798805556)
            msprime.PopulationParametersChange(
                time=1574381.9527500002,
                initial_size=160382.21744444448,
                population_id=1,
            ),
            # NCW, 1720762.798805556 generations at 158874.45340277778,
            # [1720762.798805556 - 1880600.1494444448)
            msprime.PopulationParametersChange(
                time=1720762.798805556, initial_size=158874.45340277778, population_id=1
            ),
            # NCW, 1880600.1494444448 generations at 157958.3231527778,
            # [1880600.1494444448 - 2055131.1305833336)
            msprime.PopulationParametersChange(
                time=1880600.1494444448, initial_size=157958.3231527778, population_id=1
            ),
            # NCW, 2055131.1305833336 generations at 157604.3841388889,
            # [2055131.1305833336 - 2245706.509222223)
            msprime.PopulationParametersChange(
                time=2055131.1305833336, initial_size=157604.3841388889, population_id=1
            ),
            # NCW, 2245706.509222223 generations at 157604.3841388889,
            # [2245706.509222223 - 2453801.1060555563)
            msprime.PopulationParametersChange(
                time=2245706.509222223, initial_size=157604.3841388889, population_id=1
            ),
            # NCW, 2453801.1060555563 generations at 157604.3841388889,
            # [2453801.1060555563 - 2681025.4647777784)
            msprime.PopulationParametersChange(
                time=2453801.1060555563, initial_size=157604.3841388889, population_id=1
            ),
            # NCW, 2681025.4647777784 generations at 157604.3841388889,
            # [2681025.4647777784 - inf)
            msprime.PopulationParametersChange(
                time=2681025.4647777784, initial_size=157604.3841388889, population_id=1
            ),
            # SCW, 1994.0369444444445 generations at 33850.158784722225,
            # [1994.0369444444445 - 4166.351666666667)
            msprime.PopulationParametersChange(
                time=1994.0369444444445,
                initial_size=33850.158784722225,
                population_id=2,
            ),
            # SCW, 4166.351666666667 generations at 42415.73746527778,
            # [4166.351666666667 - 6532.989166666667)
            msprime.PopulationParametersChange(
                time=4166.351666666667, initial_size=42415.73746527778, population_id=2
            ),
            # SCW, 6532.989166666667 generations at 42415.73746527778,
            # [6532.989166666667 - 9111.554375000002)
            msprime.PopulationParametersChange(
                time=6532.989166666667, initial_size=42415.73746527778, population_id=2
            ),
            # SCW, 9111.554375000002 generations at 122471.2621527778,
            # [9111.554375000002 - 11920.766458333335)
            msprime.PopulationParametersChange(
                time=9111.554375000002, initial_size=122471.2621527778, population_id=2
            ),
            # SCW, 11920.766458333335 generations at 133968.84173611112,
            # [11920.766458333335 - 14981.350208333333)
            msprime.PopulationParametersChange(
                time=11920.766458333335,
                initial_size=133968.84173611112,
                population_id=2,
            ),
            # SCW, 14981.350208333333 generations at 137971.28927083334,
            # [14981.350208333333 - 18315.813194444447)
            msprime.PopulationParametersChange(
                time=14981.350208333333,
                initial_size=137971.28927083334,
                population_id=2,
            ),
            # SCW, 18315.813194444447 generations at 143358.39802083335,
            # [18315.813194444447 - 21948.668611111112)
            msprime.PopulationParametersChange(
                time=18315.813194444447,
                initial_size=143358.39802083335,
                population_id=2,
            ),
            # SCW, 21948.668611111112 generations at 153983.5307291667,
            # [21948.668611111112 - 25906.43527777778)
            msprime.PopulationParametersChange(
                time=21948.668611111112, initial_size=153983.5307291667, population_id=2
            ),
            # SCW, 25906.43527777778 generations at 170674.89909722222,
            # [25906.43527777778 - 30218.529027777782)
            msprime.PopulationParametersChange(
                time=25906.43527777778, initial_size=170674.89909722222, population_id=2
            ),
            # SCW, 30218.529027777782 generations at 188417.3263888889,
            # [30218.529027777782 - 34916.37131944444)
            msprime.PopulationParametersChange(
                time=30218.529027777782, initial_size=188417.3263888889, population_id=2
            ),
            # SCW, 34916.37131944444 generations at 202479.8775,
            # [34916.37131944444 - 40034.503472222226)
            msprime.PopulationParametersChange(
                time=34916.37131944444, initial_size=202479.8775, population_id=2
            ),
            # SCW, 40034.503472222226 generations at 214028.48909722225,
            # [40034.503472222226 - 45610.58666666667)
            msprime.PopulationParametersChange(
                time=40034.503472222226,
                initial_size=214028.48909722225,
                population_id=2,
            ),
            # SCW, 45610.58666666667 generations at 222119.29177083337,
            # [45610.58666666667 - 51685.624791666676)
            msprime.PopulationParametersChange(
                time=45610.58666666667, initial_size=222119.29177083337, population_id=2
            ),
            # SCW, 51685.624791666676 generations at 223082.2146180556,
            # [51685.624791666676 - 58304.187291666676)
            msprime.PopulationParametersChange(
                time=51685.624791666676, initial_size=223082.2146180556, population_id=2
            ),
            # SCW, 58304.187291666676 generations at 216006.36961805559,
            # [58304.187291666676 - 65514.85486111111)
            msprime.PopulationParametersChange(
                time=58304.187291666676,
                initial_size=216006.36961805559,
                population_id=2,
            ),
            # SCW, 65514.85486111111 generations at 203076.21666666667,
            # [65514.85486111111 - 73370.88798611112)
            msprime.PopulationParametersChange(
                time=65514.85486111111, initial_size=203076.21666666667, population_id=2
            ),
            # SCW, 73370.88798611112 generations at 187279.35704861113,
            # [73370.88798611112 - 81929.78125)
            msprime.PopulationParametersChange(
                time=73370.88798611112, initial_size=187279.35704861113, population_id=2
            ),
            # SCW, 81929.78125 generations at 171195.5816319445,
            # [81929.78125 - 91254.60041666668)
            msprime.PopulationParametersChange(
                time=81929.78125, initial_size=171195.5816319445, population_id=2
            ),
            # SCW, 91254.60041666668 generations at 156701.0411805556,
            # [91254.60041666668 - 101413.75958333335)
            msprime.PopulationParametersChange(
                time=91254.60041666668, initial_size=156701.0411805556, population_id=2
            ),
            # SCW, 101413.75958333335 generations at 144929.3595138889,
            # [101413.75958333335 - 112481.91256944445)
            msprime.PopulationParametersChange(
                time=101413.75958333335, initial_size=144929.3595138889, population_id=2
            ),
            # SCW, 112481.91256944445 generations at 136399.5478125,
            # [112481.91256944445 - 124540.39861111113)
            msprime.PopulationParametersChange(
                time=112481.91256944445, initial_size=136399.5478125, population_id=2
            ),
            # SCW, 124540.39861111113 generations at 131158.73826388892,
            # [124540.39861111113 - 137677.68805555557)
            msprime.PopulationParametersChange(
                time=124540.39861111113,
                initial_size=131158.73826388892,
                population_id=2,
            ),
            # SCW, 137677.68805555557 generations at 129003.6942013889,
            # [137677.68805555557 - 151990.49659722223)
            msprime.PopulationParametersChange(
                time=137677.68805555557, initial_size=129003.6942013889, population_id=2
            ),
            # SCW, 151990.49659722223 generations at 129664.10194444447,
            # [151990.49659722223 - 167584.00812500002)
            msprime.PopulationParametersChange(
                time=151990.49659722223,
                initial_size=129664.10194444447,
                population_id=2,
            ),
            # SCW, 167584.00812500002 generations at 132892.71250000002,
            # [167584.00812500002 - 184572.76611111112)
            msprime.PopulationParametersChange(
                time=167584.00812500002,
                initial_size=132892.71250000002,
                population_id=2,
            ),
            # SCW, 184572.76611111112 generations at 138474.47829861112,
            # [184572.76611111112 - 203081.78784722227)
            msprime.PopulationParametersChange(
                time=184572.76611111112,
                initial_size=138474.47829861112,
                population_id=2,
            ),
            # SCW, 203081.78784722227 generations at 146157.35913194448,
            # [203081.78784722227 - 223246.78729166667)
            msprime.PopulationParametersChange(
                time=203081.78784722227,
                initial_size=146157.35913194448,
                population_id=2,
            ),
            # SCW, 223246.78729166667 generations at 155600.39875000002,
            # [223246.78729166667 - 245215.95784722225)
            msprime.PopulationParametersChange(
                time=223246.78729166667,
                initial_size=155600.39875000002,
                population_id=2,
            ),
            # SCW, 245215.95784722225 generations at 166325.92413194446,
            # [245215.95784722225 - 269151.08659722225)
            msprime.PopulationParametersChange(
                time=245215.95784722225,
                initial_size=166325.92413194446,
                population_id=2,
            ),
            # SCW, 269151.08659722225 generations at 177722.2196527778,
            # [269151.08659722225 - 295227.7771527778)
            msprime.PopulationParametersChange(
                time=269151.08659722225, initial_size=177722.2196527778, population_id=2
            ),
            # SCW, 295227.7771527778 generations at 189072.0515277778,
            # [295227.7771527778 - 323637.67812500003)
            msprime.PopulationParametersChange(
                time=295227.7771527778, initial_size=189072.0515277778, population_id=2
            ),
            # SCW, 323637.67812500003 generations at 199623.5332291667,
            # [323637.67812500003 - 354589.59736111114)
            msprime.PopulationParametersChange(
                time=323637.67812500003, initial_size=199623.5332291667, population_id=2
            ),
            # SCW, 354589.59736111114 generations at 208686.39548611114,
            # [354589.59736111114 - 388310.83902777784)
            msprime.PopulationParametersChange(
                time=354589.59736111114,
                initial_size=208686.39548611114,
                population_id=2,
            ),
            # SCW, 388310.83902777784 generations at 215704.30020833336,
            # [388310.83902777784 - 425049.43208333344)
            msprime.PopulationParametersChange(
                time=388310.83902777784,
                initial_size=215704.30020833336,
                population_id=2,
            ),
            # SCW, 425049.43208333344 generations at 220336.06829861112,
            # [425049.43208333344 - 465075.46736111114)
            msprime.PopulationParametersChange(
                time=425049.43208333344,
                initial_size=220336.06829861112,
                population_id=2,
            ),
            # SCW, 465075.46736111114 generations at 222488.32677083337,
            # [465075.46736111114 - 508682.6575000001)
            msprime.PopulationParametersChange(
                time=465075.46736111114,
                initial_size=222488.32677083337,
                population_id=2,
            ),
            # SCW, 508682.6575000001 generations at 222294.56111111114,
            # [508682.6575000001 - 556191.6796527778)
            msprime.PopulationParametersChange(
                time=508682.6575000001, initial_size=222294.56111111114, population_id=2
            ),
            # SCW, 556191.6796527778 generations at 220072.32861111112,
            # [556191.6796527778 - 607951.7354166667)
            msprime.PopulationParametersChange(
                time=556191.6796527778, initial_size=220072.32861111112, population_id=2
            ),
            # SCW, 607951.7354166667 generations at 216231.44531250003,
            # [607951.7354166667 - 664343.2250000001)
            msprime.PopulationParametersChange(
                time=607951.7354166667, initial_size=216231.44531250003, population_id=2
            ),
            # SCW, 664343.2250000001 generations at 211220.0569791667,
            # [664343.2250000001 - 725780.1985416667)
            msprime.PopulationParametersChange(
                time=664343.2250000001, initial_size=211220.0569791667, population_id=2
            ),
            # SCW, 725780.1985416667 generations at 205477.84118055558,
            # [725780.1985416667 - 792714.3673611112)
            msprime.PopulationParametersChange(
                time=725780.1985416667, initial_size=205477.84118055558, population_id=2
            ),
            # SCW, 792714.3673611112 generations at 199411.9397916667,
            # [792714.3673611112 - 865637.5552777778)
            msprime.PopulationParametersChange(
                time=792714.3673611112, initial_size=199411.9397916667, population_id=2
            ),
            # SCW, 865637.5552777778 generations at 193390.0507291667,
            # [865637.5552777778 - 945085.7098611112)
            msprime.PopulationParametersChange(
                time=865637.5552777778, initial_size=193390.0507291667, population_id=2
            ),
            # SCW, 945085.7098611112 generations at 187720.1488541667,
            # [945085.7098611112 - 1031642.6908333334)
            msprime.PopulationParametersChange(
                time=945085.7098611112, initial_size=187720.1488541667, population_id=2
            ),
            # SCW, 1031642.6908333334 generations at 182639.9007291667,
            # [1031642.6908333334 - 1125944.5041666667)
            msprime.PopulationParametersChange(
                time=1031642.6908333334, initial_size=182639.9007291667, population_id=2
            ),
            # SCW, 1125944.5041666667 generations at 178297.7226041667,
            # [1125944.5041666667 - 1228683.981875)
            msprime.PopulationParametersChange(
                time=1125944.5041666667, initial_size=178297.7226041667, population_id=2
            ),
            # SCW, 1228683.981875 generations at 174756.7916666667,
            # [1228683.981875 - 1340616.1303472223)
            msprime.PopulationParametersChange(
                time=1228683.981875, initial_size=174756.7916666667, population_id=2
            ),
            # SCW, 1340616.1303472223 generations at 172003.2913888889,
            # [1340616.1303472223 - 1462563.7015277778)
            msprime.PopulationParametersChange(
                time=1340616.1303472223, initial_size=172003.2913888889, population_id=2
            ),
            # SCW, 1462563.7015277778 generations at 169978.50152777781,
            # [1462563.7015277778 - 1595422.7640972224)
            msprime.PopulationParametersChange(
                time=1462563.7015277778,
                initial_size=169978.50152777781,
                population_id=2,
            ),
            # SCW, 1595422.7640972224 generations at 168579.80093750003,
            # [1595422.7640972224 - 1740169.388888889)
            msprime.PopulationParametersChange(
                time=1595422.7640972224,
                initial_size=168579.80093750003,
                population_id=2,
            ),
            # SCW, 1740169.388888889 generations at 167699.88868055557,
            # [1740169.388888889 - 1897867.671388889)
            msprime.PopulationParametersChange(
                time=1740169.388888889, initial_size=167699.88868055557, population_id=2
            ),
            # SCW, 1897867.671388889 generations at 167230.23815972224,
            # [1897867.671388889 - 2069676.1943055557)
            msprime.PopulationParametersChange(
                time=1897867.671388889, initial_size=167230.23815972224, population_id=2
            ),
            # SCW, 2069676.1943055557 generations at 167076.80784722223,
            # [2069676.1943055557 - 2256857.387152778)
            msprime.PopulationParametersChange(
                time=2069676.1943055557,
                initial_size=167076.80784722223,
                population_id=2,
            ),
            # SCW, 2256857.387152778 generations at 167462.11069444445,
            # [2256857.387152778 - 2460786.8858333337)
            msprime.PopulationParametersChange(
                time=2256857.387152778, initial_size=167462.11069444445, population_id=2
            ),
            # SCW, 2460786.8858333337 generations at 167462.11069444445,
            # [2460786.8858333337 - 2682963.560763889)
            msprime.PopulationParametersChange(
                time=2460786.8858333337,
                initial_size=167462.11069444445,
                population_id=2,
            ),
            # SCW, 2682963.560763889 generations at 167462.11069444445,
            # [2682963.560763889 - 2925019.5450000004)
            msprime.PopulationParametersChange(
                time=2682963.560763889, initial_size=167462.11069444445, population_id=2
            ),
            # SCW, 2925019.5450000004 generations at 167462.11069444445,
            # [2925019.5450000004 - inf)
            msprime.PopulationParametersChange(
                time=2925019.5450000004,
                initial_size=167462.11069444445,
                population_id=2,
            ),
            # SMW, 273.92735833333325 generations at 1779489.4096979166,
            # [273.92735833333325 - 578.0230374999999)
            msprime.PopulationParametersChange(
                time=273.92735833333325,
                initial_size=1779489.4096979166,
                population_id=3,
            ),
            # SMW, 578.0230374999999 generations at 13416.260156249999,
            # [578.0230374999999 - 915.6445999999999)
            msprime.PopulationParametersChange(
                time=578.0230374999999, initial_size=13416.260156249999, population_id=3
            ),
            # SMW, 915.6445999999999 generations at 13416.260156249999,
            # [915.6445999999999 - 1290.4231874999998)
            msprime.PopulationParametersChange(
                time=915.6445999999999, initial_size=13416.260156249999, population_id=3
            ),
            # SMW, 1290.4231874999998 generations at 4897.875385416666,
            # [1290.4231874999998 - 1706.4624874999995)
            msprime.PopulationParametersChange(
                time=1290.4231874999998, initial_size=4897.875385416666, population_id=3
            ),
            # SMW, 1706.4624874999995 generations at 4407.559341666665,
            # [1706.4624874999995 - 2168.3636041666664)
            msprime.PopulationParametersChange(
                time=1706.4624874999995, initial_size=4407.559341666665, population_id=3
            ),
            # SMW, 2168.3636041666664 generations at 4602.385014583333,
            # [2168.3636041666664 - 2681.1255749999996)
            msprime.PopulationParametersChange(
                time=2168.3636041666664, initial_size=4602.385014583333, population_id=3
            ),
            # SMW, 2681.1255749999996 generations at 4882.169454166666,
            # [2681.1255749999996 - 3250.3443374999993)
            msprime.PopulationParametersChange(
                time=2681.1255749999996, initial_size=4882.169454166666, population_id=3
            ),
            # SMW, 3250.3443374999993 generations at 5175.035952083333,
            # [3250.3443374999993 - 3882.2873416666657)
            msprime.PopulationParametersChange(
                time=3250.3443374999993, initial_size=5175.035952083333, population_id=3
            ),
            # SMW, 3882.2873416666657 generations at 5521.001679166666,
            # [3882.2873416666657 - 4583.818937499999)
            msprime.PopulationParametersChange(
                time=3882.2873416666657, initial_size=5521.001679166666, population_id=3
            ),
            # SMW, 4583.818937499999 generations at 5931.532089583332,
            # [4583.818937499999 - 5362.624212499999)
            msprime.PopulationParametersChange(
                time=4583.818937499999, initial_size=5931.532089583332, population_id=3
            ),
            # SMW, 5362.624212499999 generations at 6379.940779166665,
            # [5362.624212499999 - 6227.2089916666655)
            msprime.PopulationParametersChange(
                time=5362.624212499999, initial_size=6379.940779166665, population_id=3
            ),
            # SMW, 6227.2089916666655 generations at 6835.213818749998,
            # [6227.2089916666655 - 7186.999320833332)
            msprime.PopulationParametersChange(
                time=6227.2089916666655, initial_size=6835.213818749998, population_id=3
            ),
            # SMW, 7186.999320833332 generations at 7291.755270833331,
            # [7186.999320833332 - 8252.540433333332)
            msprime.PopulationParametersChange(
                time=7186.999320833332, initial_size=7291.755270833331, population_id=3
            ),
            # SMW, 8252.540433333332 generations at 7769.623462499998,
            # [8252.540433333332 - 9435.422137499998)
            msprime.PopulationParametersChange(
                time=8252.540433333332, initial_size=7769.623462499998, population_id=3
            ),
            # SMW, 9435.422137499998 generations at 8299.036456249998,
            # [9435.422137499998 - 10748.577266666665)
            msprime.PopulationParametersChange(
                time=9435.422137499998, initial_size=8299.036456249998, population_id=3
            ),
            # SMW, 10748.577266666665 generations at 8915.174045833332,
            # [10748.577266666665 - 12206.381162499996)
            msprime.PopulationParametersChange(
                time=10748.577266666665, initial_size=8915.174045833332, population_id=3
            ),
            # SMW, 12206.381162499996 generations at 9670.289852083331,
            # [12206.381162499996 - 13824.751158333329)
            msprime.PopulationParametersChange(
                time=12206.381162499996, initial_size=9670.289852083331, population_id=3
            ),
            # SMW, 13824.751158333329 generations at 10648.260760416664,
            # [13824.751158333329 - 15621.34554583333)
            msprime.PopulationParametersChange(
                time=13824.751158333329,
                initial_size=10648.260760416664,
                population_id=3,
            ),
            # SMW, 15621.34554583333 generations at 11973.826435416664,
            # [15621.34554583333 - 17615.862024999995)
            msprime.PopulationParametersChange(
                time=15621.34554583333, initial_size=11973.826435416664, population_id=3
            ),
            # SMW, 17615.862024999995 generations at 13810.972716666665,
            # [17615.862024999995 - 19830.03770416666)
            msprime.PopulationParametersChange(
                time=17615.862024999995,
                initial_size=13810.972716666665,
                population_id=3,
            ),
            # SMW, 19830.03770416666 generations at 16361.98652708333,
            # [19830.03770416666 - 22288.071904166663)
            msprime.PopulationParametersChange(
                time=19830.03770416666, initial_size=16361.98652708333, population_id=3
            ),
            # SMW, 22288.071904166663 generations at 19850.04628958333,
            # [22288.071904166663 - 25016.849995833327)
            msprime.PopulationParametersChange(
                time=22288.071904166663, initial_size=19850.04628958333, population_id=3
            ),
            # SMW, 25016.849995833327 generations at 24476.406787499996,
            # [25016.849995833327 - 28046.19210833333)
            msprime.PopulationParametersChange(
                time=25016.849995833327,
                initial_size=24476.406787499996,
                population_id=3,
            ),
            # SMW, 28046.19210833333 generations at 30345.97319583333,
            # [28046.19210833333 - 31409.17644999999)
            msprime.PopulationParametersChange(
                time=28046.19210833333, initial_size=30345.97319583333, population_id=3
            ),
            # SMW, 31409.17644999999 generations at 37379.656264583326,
            # [31409.17644999999 - 35142.562112499996)
            msprime.PopulationParametersChange(
                time=31409.17644999999, initial_size=37379.656264583326, population_id=3
            ),
            # SMW, 35142.562112499996 generations at 45251.64807708333,
            # [35142.562112499996 - 39287.162133333324)
            msprime.PopulationParametersChange(
                time=35142.562112499996, initial_size=45251.64807708333, population_id=3
            ),
            # SMW, 39287.162133333324 generations at 53429.50139791666,
            # [39287.162133333324 - 43888.26629999999)
            msprime.PopulationParametersChange(
                time=39287.162133333324, initial_size=53429.50139791666, population_id=3
            ),
            # SMW, 43888.26629999999 generations at 61330.815922916656,
            # [43888.26629999999 - 48996.11369583332)
            msprime.PopulationParametersChange(
                time=43888.26629999999, initial_size=61330.815922916656, population_id=3
            ),
            # SMW, 48996.11369583332 generations at 68519.94169999998,
            # [48996.11369583332 - 54666.589083333325)
            msprime.PopulationParametersChange(
                time=48996.11369583332, initial_size=68519.94169999998, population_id=3
            ),
            # SMW, 54666.589083333325 generations at 74824.37462916666,
            # [54666.589083333325 - 60961.59596666666)
            msprime.PopulationParametersChange(
                time=54666.589083333325, initial_size=74824.37462916666, population_id=3
            ),
            # SMW, 60961.59596666666 generations at 80318.72721041665,
            # [60961.59596666666 - 67949.92707083332)
            msprime.PopulationParametersChange(
                time=60961.59596666666, initial_size=80318.72721041665, population_id=3
            ),
            # SMW, 67949.92707083332 generations at 85281.27919791665,
            # [67949.92707083332 - 75707.98559583332)
            msprime.PopulationParametersChange(
                time=67949.92707083332, initial_size=85281.27919791665, population_id=3
            ),
            # SMW, 75707.98559583332 generations at 90121.39207291666,
            # [75707.98559583332 - 84320.53134166665)
            msprime.PopulationParametersChange(
                time=75707.98559583332, initial_size=90121.39207291666, population_id=3
            ),
            # SMW, 84320.53134166665 generations at 95327.01914999998,
            # [84320.53134166665 - 93881.70041249998)
            msprime.PopulationParametersChange(
                time=84320.53134166665, initial_size=95327.01914999998, population_id=3
            ),
            # SMW, 93881.70041249998 generations at 101382.34581249999,
            # [93881.70041249998 - 104495.92543749999)
            msprime.PopulationParametersChange(
                time=93881.70041249998, initial_size=101382.34581249999, population_id=3
            ),
            # SMW, 104495.92543749999 generations at 108697.1563458333,
            # [104495.92543749999 - 116279.22885416665)
            msprime.PopulationParametersChange(
                time=104495.92543749999, initial_size=108697.1563458333, population_id=3
            ),
            # SMW, 116279.22885416665 generations at 117580.7444333333,
            # [116279.22885416665 - 129360.36696666664)
            msprime.PopulationParametersChange(
                time=116279.22885416665, initial_size=117580.7444333333, population_id=3
            ),
            # SMW, 129360.36696666664 generations at 128223.79475416664,
            # [129360.36696666664 - 143882.27245416664)
            msprime.PopulationParametersChange(
                time=129360.36696666664,
                initial_size=128223.79475416664,
                population_id=3,
            ),
            # SMW, 143882.27245416664 generations at 140717.28481666665,
            # [143882.27245416664 - 160003.67097499996)
            msprime.PopulationParametersChange(
                time=143882.27245416664,
                initial_size=140717.28481666665,
                population_id=3,
            ),
            # SMW, 160003.67097499996 generations at 155055.94200416663,
            # [160003.67097499996 - 177900.67289999998)
            msprime.PopulationParametersChange(
                time=160003.67097499996,
                initial_size=155055.94200416663,
                population_id=3,
            ),
            # SMW, 177900.67289999998 generations at 171127.08900624997,
            # [177900.67289999998 - 197768.9122041666)
            msprime.PopulationParametersChange(
                time=177900.67289999998,
                initial_size=171127.08900624997,
                population_id=3,
            ),
            # SMW, 197768.9122041666 generations at 188647.94444791664,
            # [197768.9122041666 - 219825.4863916666)
            msprime.PopulationParametersChange(
                time=197768.9122041666, initial_size=188647.94444791664, population_id=3
            ),
            # SMW, 219825.4863916666 generations at 207099.76492291663,
            # [219825.4863916666 - 244311.3689666666)
            msprime.PopulationParametersChange(
                time=219825.4863916666, initial_size=207099.76492291663, population_id=3
            ),
            # SMW, 244311.3689666666 generations at 225653.4065895833,
            # [244311.3689666666 - 271494.21983749996)
            msprime.PopulationParametersChange(
                time=244311.3689666666, initial_size=225653.4065895833, population_id=3
            ),
            # SMW, 271494.21983749996 generations at 243141.84289999996,
            # [271494.21983749996 - 301670.9967541666)
            msprime.PopulationParametersChange(
                time=271494.21983749996,
                initial_size=243141.84289999996,
                population_id=3,
            ),
            # SMW, 301670.9967541666 generations at 258150.88096458328,
            # [301670.9967541666 - 335171.5118375)
            msprime.PopulationParametersChange(
                time=301670.9967541666, initial_size=258150.88096458328, population_id=3
            ),
            # SMW, 335171.5118375 generations at 269176.8426354166,
            # [335171.5118375 - 372361.8637541666)
            msprime.PopulationParametersChange(
                time=335171.5118375, initial_size=269176.8426354166, population_id=3
            ),
            # SMW, 372361.8637541666 generations at 274910.6764708333,
            # [372361.8637541666 - 413648.3921791666)
            msprime.PopulationParametersChange(
                time=372361.8637541666, initial_size=274910.6764708333, population_id=3
            ),
            # SMW, 413648.3921791666 generations at 274517.54320833326,
            # [413648.3921791666 - 459482.3037708332)
            msprime.PopulationParametersChange(
                time=413648.3921791666, initial_size=274517.54320833326, population_id=3
            ),
            # SMW, 459482.3037708332 generations at 267869.51971666666,
            # [459482.3037708332 - 510364.4473708332)
            msprime.PopulationParametersChange(
                time=459482.3037708332, initial_size=267869.51971666666, population_id=3
            ),
            # SMW, 510364.4473708332 generations at 255602.0060458333,
            # [510364.4473708332 - 566850.8353291665)
            msprime.PopulationParametersChange(
                time=510364.4473708332, initial_size=255602.0060458333, population_id=3
            ),
            # SMW, 566850.8353291665 generations at 238998.64808124994,
            # [566850.8353291665 - 629558.7368583332)
            msprime.PopulationParametersChange(
                time=566850.8353291665, initial_size=238998.64808124994, population_id=3
            ),
            # SMW, 629558.7368583332 generations at 219719.26306249996,
            # [629558.7368583332 - 699173.4180291665)
            msprime.PopulationParametersChange(
                time=629558.7368583332, initial_size=219719.26306249996, population_id=3
            ),
            # SMW, 699173.4180291665 generations at 177796.99883124998,
            # [699173.4180291665 - 776455.5532791665)
            msprime.PopulationParametersChange(
                time=699173.4180291665, initial_size=177796.99883124998, population_id=3
            ),
            # SMW, 776455.5532791665 generations at 177796.99883124998,
            # [776455.5532791665 - 862249.7063666665)
            msprime.PopulationParametersChange(
                time=776455.5532791665, initial_size=177796.99883124998, population_id=3
            ),
            # SMW, 862249.7063666665 generations at 177796.99883124998,
            # [862249.7063666665 - 957493.3833541665)
            msprime.PopulationParametersChange(
                time=862249.7063666665, initial_size=177796.99883124998, population_id=3
            ),
            # SMW, 957493.3833541665 generations at 177796.99883124998,
            # [957493.3833541665 - inf)
            msprime.PopulationParametersChange(
                time=957493.3833541665, initial_size=177796.99883124998, population_id=3
            ),
        ],
    )


_species.add_demographic_model(_WildBoar_4Z22())

_WangEtAl = stdpopsim.Citation(
    doi="https://www.sciencedirect.com/science/article/pii/S2666979X25002101?via%3Dihub",
    year=2025,
    author="Wang et al.",
    reasons={stdpopsim.CiteReason.DEM_MODEL},
    # Figure 2b Demographic history of four geographical populations using MSMC.
)


def _WildBoar_13W25():
    id = "WildBoar_13W25"
    description = "Piecewise size model for wild boar (Wang et al. 2025)"
    long_description = """
    This demographic model is a piecewise size model for wild boar
    from Wang et al. (2025). Effective population sizes were estimated by
    MSMC with parameters n = 64 (4 + 50 + 4 + 6), generation time (g) = 5 years, and
    mutation rate (u) = 2.5e-8 per generation.
    """
    populations = [
        stdpopsim.Population(id="FJW", description="FuJian wild boar"),
        stdpopsim.Population(id="HBW", description="Hubei wild boar"),
        stdpopsim.Population(id="HW", description="Hungary wild boar"),
        stdpopsim.Population(id="JLW", description="Jilin wild boar"),
        stdpopsim.Population(id="JXW", description="Jiangxi wild boar"),
        stdpopsim.Population(id="HLJW", description="HeiLongJiang wild boar"),
        stdpopsim.Population(id="LNW", description="Liaoning wild boar"),
        stdpopsim.Population(id="IMW", description="Inner Mongolia wild boar"),
        stdpopsim.Population(id="NW", description="Netherland wild boar"),
        stdpopsim.Population(id="SXW", description="Shaanxi wild boar"),
        stdpopsim.Population(id="XJW", description="XinJiang wild boar"),
        stdpopsim.Population(id="YNW", description="Yunnan wild boar"),
        stdpopsim.Population(id="SCW", description="Sichuan wild boar"),
    ]
    citations = [
        _WangEtAl,
    ]
    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=5,
        mutation_rate=2.5e-8,
        population_configurations=[
            stdpopsim.PopulationConfiguration(
                # FJW,0.0 generations at 29181.701322368797,
                # [0.0 - 654.0)
                initial_size=29181.701322368797,
                metadata=populations[0].asdict(),
            ),
            stdpopsim.PopulationConfiguration(
                # HBW,0.0 generations at 122440.23386084668,
                # [0.0 - 519.0)
                initial_size=122440.23386084668,
                metadata=populations[1].asdict(),
            ),
            stdpopsim.PopulationConfiguration(
                # HungryW,0.0 generations at 26924.713117181735,
                # [0.0 - 243.0)
                initial_size=26924.713117181735,
                metadata=populations[2].asdict(),
            ),
            stdpopsim.PopulationConfiguration(
                # JLW,0.0 generations at 3268913.936033892,
                # [0.0 - 216.0)
                initial_size=3268913.936033892,
                metadata=populations[3].asdict(),
            ),
            stdpopsim.PopulationConfiguration(
                # JXW,0.0 generations at 65893.73317650624,
                # [0.0 - 618.0)
                initial_size=65893.73317650624,
                metadata=populations[4].asdict(),
            ),
            stdpopsim.PopulationConfiguration(
                # HLJW,0.0 generations at 24816.94401669684,
                # [0.0 - 631.0)
                initial_size=24816.94401669684,
                metadata=populations[5].asdict(),
            ),
            stdpopsim.PopulationConfiguration(
                # LNW,0.0 generations at 42760.35174665347,
                # [0.0 - 576.0)
                initial_size=42760.35174665347,
                metadata=populations[6].asdict(),
            ),
            stdpopsim.PopulationConfiguration(
                # MW,0.0 generations at 475403.7960993119,
                # [0.0 - 513.0)
                initial_size=475403.7960993119,
                metadata=populations[7].asdict(),
            ),
            stdpopsim.PopulationConfiguration(
                # NetherlandW,0.0 generations at 14409.948628533139,
                # [0.0 - 187.0)
                initial_size=14409.948628533139,
                metadata=populations[8].asdict(),
            ),
            stdpopsim.PopulationConfiguration(
                # SXW,0.0 generations at 56807.36905191342,
                # [0.0 - 628.0)
                initial_size=56807.36905191342,
                metadata=populations[9].asdict(),
            ),
            stdpopsim.PopulationConfiguration(
                # XJW,0.0 generations at 32182.901867895627,
                # [0.0 - 414.0)
                initial_size=32182.901867895627,
                metadata=populations[10].asdict(),
            ),
            stdpopsim.PopulationConfiguration(
                # YNW,0.0 generations at 65330.22796983051,
                # [0.0 - 621.0)
                initial_size=65330.22796983051,
                metadata=populations[11].asdict(),
            ),
            stdpopsim.PopulationConfiguration(
                # SCW,0.0 generations at 51760.3701901676,
                # [0.0 - 584.0)
                initial_size=51760.3701901676,
                metadata=populations[12].asdict(),
            ),
        ],
        demographic_events=[
            # FJW, 654.0 generations at 29181.701322368797,
            # [654.0 - 1462.0)
            msprime.PopulationParametersChange(
                time=654.0, initial_size=29181.701322368797, population_id=0
            ),
            # FJW, 1462.0 generations at 10149.39915556999,
            # [1462.0 - 2462.0)
            msprime.PopulationParametersChange(
                time=1462.0, initial_size=10149.39915556999, population_id=0
            ),
            # FJW, 2462.0 generations at 9230.669823555747,
            # [2462.0 - 3698.0)
            msprime.PopulationParametersChange(
                time=2462.0, initial_size=9230.669823555747, population_id=0
            ),
            # FJW, 3698.0 generations at 10026.520145785604,
            # [3698.0 - 5228.0)
            msprime.PopulationParametersChange(
                time=3698.0, initial_size=10026.520145785604, population_id=0
            ),
            # FJW, 5228.0 generations at 11542.545823906921,
            # [5228.0 - 7120.0)
            msprime.PopulationParametersChange(
                time=5228.0, initial_size=11542.545823906921, population_id=0
            ),
            # FJW, 7120.0 generations at 13979.366455112255,
            # [7120.0 - 9461.0)
            msprime.PopulationParametersChange(
                time=7120.0, initial_size=13979.366455112255, population_id=0
            ),
            # FJW, 9461.0 generations at 17792.48623306378,
            # [9461.0 - 12355.0)
            msprime.PopulationParametersChange(
                time=9461.0, initial_size=17792.48623306378, population_id=0
            ),
            # FJW, 12355.0 generations at 22898.753162890283,
            # [12355.0 - 15936.0)
            msprime.PopulationParametersChange(
                time=12355.0, initial_size=22898.753162890283, population_id=0
            ),
            # FJW, 15936.0 generations at 28211.410669273402,
            # [15936.0 - 20364.0)
            msprime.PopulationParametersChange(
                time=15936.0, initial_size=28211.410669273402, population_id=0
            ),
            # FJW, 20364.0 generations at 32962.12487041765,
            # [20364.0 - 25842.0)
            msprime.PopulationParametersChange(
                time=20364.0, initial_size=32962.12487041765, population_id=0
            ),
            # FJW, 25842.0 generations at 37325.89801444885,
            # [25842.0 - 32618.0)
            msprime.PopulationParametersChange(
                time=25842.0, initial_size=37325.89801444885, population_id=0
            ),
            # FJW, 32618.0 generations at 41888.93985376571,
            # [32618.0 - 40999.0)
            msprime.PopulationParametersChange(
                time=32618.0, initial_size=41888.93985376571, population_id=0
            ),
            # FJW, 40999.0 generations at 47822.74976028847,
            # [40999.0 - 51365.0)
            msprime.PopulationParametersChange(
                time=40999.0, initial_size=47822.74976028847, population_id=0
            ),
            # FJW, 51365.0 generations at 55348.91958908962,
            # [51365.0 - 64187.0)
            msprime.PopulationParametersChange(
                time=51365.0, initial_size=55348.91958908962, population_id=0
            ),
            # FJW, 64187.0 generations at 62616.81950388694,
            # [64187.0 - 80046.0)
            msprime.PopulationParametersChange(
                time=64187.0, initial_size=62616.81950388694, population_id=0
            ),
            # FJW, 80046.0 generations at 66300.24729992243,
            # [80046.0 - 99663.0)
            msprime.PopulationParametersChange(
                time=80046.0, initial_size=66300.24729992243, population_id=0
            ),
            # FJW, 99663.0 generations at 63675.00278578138,
            # [99663.0 - 123927.0)
            msprime.PopulationParametersChange(
                time=99663.0, initial_size=63675.00278578138, population_id=0
            ),
            # FJW, 123927.0 generations at 54931.92561118634,
            # [123927.0 - 153939.0)
            msprime.PopulationParametersChange(
                time=123927.0, initial_size=54931.92561118634, population_id=0
            ),
            # FJW, 153939.0 generations at 43695.59897927081,
            # [153939.0 - 191061.0)
            msprime.PopulationParametersChange(
                time=153939.0, initial_size=43695.59897927081, population_id=0
            ),
            # FJW, 191061.0 generations at 34221.143875955204,
            # [191061.0 - 236977.0)
            msprime.PopulationParametersChange(
                time=191061.0, initial_size=34221.143875955204, population_id=0
            ),
            # FJW, 236977.0 generations at 28649.477934888335,
            # [236977.0 - 293771.0)
            msprime.PopulationParametersChange(
                time=236977.0, initial_size=28649.477934888335, population_id=0
            ),
            # FJW, 293771.0 generations at 26770.60801404922,
            # [293771.0 - 364019.0)
            msprime.PopulationParametersChange(
                time=293771.0, initial_size=26770.60801404922, population_id=0
            ),
            # FJW, 364019.0 generations at 27672.889838376486,
            # [364019.0 - 450908.0)
            msprime.PopulationParametersChange(
                time=364019.0, initial_size=27672.889838376486, population_id=0
            ),
            # FJW, 450908.0 generations at 30313.410349604565,
            # [450908.0 - 558384.0)
            msprime.PopulationParametersChange(
                time=450908.0, initial_size=30313.410349604565, population_id=0
            ),
            # FJW, 558384.0 generations at 34391.44340887987,
            # [558384.0 - 691320.0)
            msprime.PopulationParametersChange(
                time=558384.0, initial_size=34391.44340887987, population_id=0
            ),
            # FJW, 691320.0 generations at 39719.02759876633,
            # [691320.0 - 855748.0)
            msprime.PopulationParametersChange(
                time=691320.0, initial_size=39719.02759876633, population_id=0
            ),
            # FJW, 855748.0 generations at 46250.04625004625,
            # [855748.0 - 1059128.0)
            msprime.PopulationParametersChange(
                time=855748.0, initial_size=46250.04625004625, population_id=0
            ),
            # FJW, 1059128.0 generations at 46250.04625004625,
            # [1059128.0 - 1310692.0)
            msprime.PopulationParametersChange(
                time=1059128.0, initial_size=46250.04625004625, population_id=0
            ),
            # FJW, 1310692.0 generations at 50311.8074265259,
            # [1310692.0 - 1621852.0)
            msprime.PopulationParametersChange(
                time=1310692.0, initial_size=50311.8074265259, population_id=0
            ),
            # FJW, 1621852.0 generations at 50311.8074265259,
            # [1621852.0 - 2006724.0)
            msprime.PopulationParametersChange(
                time=1621852.0, initial_size=50311.8074265259, population_id=0
            ),
            # FJW, 2006724.0 generations at 50311.8074265259,
            # [2006724.0 - inf)
            msprime.PopulationParametersChange(
                time=2006724.0, initial_size=50311.8074265259, population_id=0
            ),
            # HBW, 519.0 generations at 122440.23386084668,
            # [519.0 - 1162.0)
            msprime.PopulationParametersChange(
                time=519.0, initial_size=122440.23386084668, population_id=1
            ),
            # HBW, 1162.0 generations at 15583.606046439147,
            # [1162.0 - 1956.0)
            msprime.PopulationParametersChange(
                time=1162.0, initial_size=15583.606046439147, population_id=1
            ),
            # HBW, 1956.0 generations at 10614.584439019212,
            # [1956.0 - 2939.0)
            msprime.PopulationParametersChange(
                time=1956.0, initial_size=10614.584439019212, population_id=1
            ),
            # HBW, 2939.0 generations at 11810.140186364013,
            # [2939.0 - 4155.0)
            msprime.PopulationParametersChange(
                time=2939.0, initial_size=11810.140186364013, population_id=1
            ),
            # HBW, 4155.0 generations at 14705.449839710596,
            # [4155.0 - 5659.0)
            msprime.PopulationParametersChange(
                time=4155.0, initial_size=14705.449839710596, population_id=1
            ),
            # HBW, 5659.0 generations at 17962.834894603067,
            # [5659.0 - 7519.0)
            msprime.PopulationParametersChange(
                time=5659.0, initial_size=17962.834894603067, population_id=1
            ),
            # HBW, 7519.0 generations at 22598.02741818667,
            # [7519.0 - 9819.0)
            msprime.PopulationParametersChange(
                time=7519.0, initial_size=22598.02741818667, population_id=1
            ),
            # HBW, 9819.0 generations at 29704.04375999727,
            # [9819.0 - 12665.0)
            msprime.PopulationParametersChange(
                time=9819.0, initial_size=29704.04375999727, population_id=1
            ),
            # HBW, 12665.0 generations at 39308.63964590778,
            # [12665.0 - 16184.0)
            msprime.PopulationParametersChange(
                time=12665.0, initial_size=39308.63964590778, population_id=1
            ),
            # HBW, 16184.0 generations at 46987.07150727478,
            # [16184.0 - 20538.0)
            msprime.PopulationParametersChange(
                time=16184.0, initial_size=46987.07150727478, population_id=1
            ),
            # HBW, 20538.0 generations at 48346.19745070501,
            # [20538.0 - 25923.0)
            msprime.PopulationParametersChange(
                time=20538.0, initial_size=48346.19745070501, population_id=1
            ),
            # HBW, 25923.0 generations at 45711.987054365265,
            # [25923.0 - 32583.0)
            msprime.PopulationParametersChange(
                time=25923.0, initial_size=45711.987054365265, population_id=1
            ),
            # HBW, 32583.0 generations at 44650.53145295062,
            # [32583.0 - 40821.0)
            msprime.PopulationParametersChange(
                time=32583.0, initial_size=44650.53145295062, population_id=1
            ),
            # HBW, 40821.0 generations at 48088.94531325139,
            # [40821.0 - 51011.0)
            msprime.PopulationParametersChange(
                time=40821.0, initial_size=48088.94531325139, population_id=1
            ),
            # HBW, 51011.0 generations at 55654.96150068038,
            # [51011.0 - 63615.0)
            msprime.PopulationParametersChange(
                time=51011.0, initial_size=55654.96150068038, population_id=1
            ),
            # HBW, 63615.0 generations at 64661.51314406909,
            # [63615.0 - 79205.0)
            msprime.PopulationParametersChange(
                time=63615.0, initial_size=64661.51314406909, population_id=1
            ),
            # HBW, 79205.0 generations at 70504.4594070575,
            # [79205.0 - 98488.0)
            msprime.PopulationParametersChange(
                time=79205.0, initial_size=70504.4594070575, population_id=1
            ),
            # HBW, 98488.0 generations at 68972.65234334586,
            # [98488.0 - 122340.0)
            msprime.PopulationParametersChange(
                time=98488.0, initial_size=68972.65234334586, population_id=1
            ),
            # HBW, 122340.0 generations at 59468.76551763101,
            # [122340.0 - 151842.0)
            msprime.PopulationParametersChange(
                time=122340.0, initial_size=59468.76551763101, population_id=1
            ),
            # HBW, 151842.0 generations at 45961.37406123894,
            # [151842.0 - 188333.0)
            msprime.PopulationParametersChange(
                time=151842.0, initial_size=45961.37406123894, population_id=1
            ),
            # HBW, 188333.0 generations at 34138.6061548493,
            # [188333.0 - 233469.0)
            msprime.PopulationParametersChange(
                time=188333.0, initial_size=34138.6061548493, population_id=1
            ),
            # HBW, 233469.0 generations at 27039.414001819754,
            # [233469.0 - 289297.0)
            msprime.PopulationParametersChange(
                time=233469.0, initial_size=27039.414001819754, population_id=1
            ),
            # HBW, 289297.0 generations at 23650.558271428003,
            # [289297.0 - 358352.0)
            msprime.PopulationParametersChange(
                time=289297.0, initial_size=23650.558271428003, population_id=1
            ),
            # HBW, 358352.0 generations at 22930.730848253595,
            # [358352.0 - 443764.0)
            msprime.PopulationParametersChange(
                time=358352.0, initial_size=22930.730848253595, population_id=1
            ),
            # HBW, 443764.0 generations at 24326.253987985263,
            # [443764.0 - 549412.0)
            msprime.PopulationParametersChange(
                time=443764.0, initial_size=24326.253987985263, population_id=1
            ),
            # HBW, 549412.0 generations at 27050.824441501918,
            # [549412.0 - 680088.0)
            msprime.PopulationParametersChange(
                time=549412.0, initial_size=27050.824441501918, population_id=1
            ),
            # HBW, 680088.0 generations at 32548.692844495363,
            # [680088.0 - 841724.0)
            msprime.PopulationParametersChange(
                time=680088.0, initial_size=32548.692844495363, population_id=1
            ),
            # HBW, 841724.0 generations at 32548.692844495363,
            # [841724.0 - 1041648.0)
            msprime.PopulationParametersChange(
                time=841724.0, initial_size=32548.692844495363, population_id=1
            ),
            # HBW, 1041648.0 generations at 47873.00249897073,
            # [1041648.0 - 1288936.0)
            msprime.PopulationParametersChange(
                time=1041648.0, initial_size=47873.00249897073, population_id=1
            ),
            # HBW, 1288936.0 generations at 47873.00249897073,
            # [1288936.0 - 1594804.0)
            msprime.PopulationParametersChange(
                time=1288936.0, initial_size=47873.00249897073, population_id=1
            ),
            # HBW, 1594804.0 generations at 47873.00249897073,
            # [1594804.0 - inf)
            msprime.PopulationParametersChange(
                time=1594804.0, initial_size=47873.00249897073, population_id=1
            ),
            # HungryW, 243.0 generations at 26924.713117181735,
            # [243.0 - 544.0)
            msprime.PopulationParametersChange(
                time=243.0, initial_size=26924.713117181735, population_id=2
            ),
            # HungryW, 544.0 generations at 4435.022485564002,
            # [544.0 - 916.0)
            msprime.PopulationParametersChange(
                time=544.0, initial_size=4435.022485564002, population_id=2
            ),
            # HungryW, 916.0 generations at 3300.607806927646,
            # [916.0 - 1377.0)
            msprime.PopulationParametersChange(
                time=916.0, initial_size=3300.607806927646, population_id=2
            ),
            # HungryW, 1377.0 generations at 3314.721173941651,
            # [1377.0 - 1946.0)
            msprime.PopulationParametersChange(
                time=1377.0, initial_size=3314.721173941651, population_id=2
            ),
            # HungryW, 1946.0 generations at 3976.0087631233137,
            # [1946.0 - 2651.0)
            msprime.PopulationParametersChange(
                time=1946.0, initial_size=3976.0087631233137, population_id=2
            ),
            # HungryW, 2651.0 generations at 5197.20805983026,
            # [2651.0 - 3522.0)
            msprime.PopulationParametersChange(
                time=2651.0, initial_size=5197.20805983026, population_id=2
            ),
            # HungryW, 3522.0 generations at 6140.657910088487,
            # [3522.0 - 4600.0)
            msprime.PopulationParametersChange(
                time=3522.0, initial_size=6140.657910088487, population_id=2
            ),
            # HungryW, 4600.0 generations at 6407.6481688543445,
            # [4600.0 - 5933.0)
            msprime.PopulationParametersChange(
                time=4600.0, initial_size=6407.6481688543445, population_id=2
            ),
            # HungryW, 5933.0 generations at 6407.873995565752,
            # [5933.0 - 7581.0)
            msprime.PopulationParametersChange(
                time=5933.0, initial_size=6407.873995565752, population_id=2
            ),
            # HungryW, 7581.0 generations at 7112.400826460977,
            # [7581.0 - 9621.0)
            msprime.PopulationParametersChange(
                time=7581.0, initial_size=7112.400826460977, population_id=2
            ),
            # HungryW, 9621.0 generations at 9937.246289680666,
            # [9621.0 - 12143.0)
            msprime.PopulationParametersChange(
                time=9621.0, initial_size=9937.246289680666, population_id=2
            ),
            # HungryW, 12143.0 generations at 17226.380479065643,
            # [12143.0 - 15263.0)
            msprime.PopulationParametersChange(
                time=12143.0, initial_size=17226.380479065643, population_id=2
            ),
            # HungryW, 15263.0 generations at 31222.728897728393,
            # [15263.0 - 19122.0)
            msprime.PopulationParametersChange(
                time=15263.0, initial_size=31222.728897728393, population_id=2
            ),
            # HungryW, 19122.0 generations at 48107.10566004152,
            # [19122.0 - 23896.0)
            msprime.PopulationParametersChange(
                time=19122.0, initial_size=48107.10566004152, population_id=2
            ),
            # HungryW, 23896.0 generations at 59596.3538950687,
            # [23896.0 - 29800.0)
            msprime.PopulationParametersChange(
                time=23896.0, initial_size=59596.3538950687, population_id=2
            ),
            # HungryW, 29800.0 generations at 64567.57480960637,
            # [29800.0 - 37103.0)
            msprime.PopulationParametersChange(
                time=29800.0, initial_size=64567.57480960637, population_id=2
            ),
            # HungryW, 37103.0 generations at 69362.79865019994,
            # [37103.0 - 46136.0)
            msprime.PopulationParametersChange(
                time=37103.0, initial_size=69362.79865019994, population_id=2
            ),
            # HungryW, 46136.0 generations at 80290.0075071157,
            # [46136.0 - 57308.0)
            msprime.PopulationParametersChange(
                time=46136.0, initial_size=80290.0075071157, population_id=2
            ),
            # HungryW, 57308.0 generations at 95366.61310241897,
            # [57308.0 - 71128.0)
            msprime.PopulationParametersChange(
                time=57308.0, initial_size=95366.61310241897, population_id=2
            ),
            # HungryW, 71128.0 generations at 108946.10980678408,
            # [71128.0 - 88222.0)
            msprime.PopulationParametersChange(
                time=71128.0, initial_size=108946.10980678408, population_id=2
            ),
            # HungryW, 88222.0 generations at 109678.58690108638,
            # [88222.0 - 109365.0)
            msprime.PopulationParametersChange(
                time=88222.0, initial_size=109678.58690108638, population_id=2
            ),
            # HungryW, 109365.0 generations at 89326.38969530768,
            # [109365.0 - 135517.0)
            msprime.PopulationParametersChange(
                time=109365.0, initial_size=89326.38969530768, population_id=2
            ),
            # HungryW, 135517.0 generations at 58610.978422368295,
            # [135517.0 - 167864.0)
            msprime.PopulationParametersChange(
                time=135517.0, initial_size=58610.978422368295, population_id=2
            ),
            # HungryW, 167864.0 generations at 32522.176058800094,
            # [167864.0 - 207875.0)
            msprime.PopulationParametersChange(
                time=167864.0, initial_size=32522.176058800094, population_id=2
            ),
            # HungryW, 207875.0 generations at 20192.556216076508,
            # [207875.0 - 257364.0)
            msprime.PopulationParametersChange(
                time=207875.0, initial_size=20192.556216076508, population_id=2
            ),
            # HungryW, 257364.0 generations at 15234.845137799177,
            # [257364.0 - 318578.0)
            msprime.PopulationParametersChange(
                time=257364.0, initial_size=15234.845137799177, population_id=2
            ),
            # HungryW, 318578.0 generations at 11981.069909542923,
            # [318578.0 - 394293.0)
            msprime.PopulationParametersChange(
                time=318578.0, initial_size=11981.069909542923, population_id=2
            ),
            # HungryW, 394293.0 generations at 11981.069909542923,
            # [394293.0 - 487944.0)
            msprime.PopulationParametersChange(
                time=394293.0, initial_size=11981.069909542923, population_id=2
            ),
            # HungryW, 487944.0 generations at 181.48820326678768,
            # [487944.0 - 603784.0)
            msprime.PopulationParametersChange(
                time=487944.0, initial_size=181.48820326678768, population_id=2
            ),
            # HungryW, 603784.0 generations at 181.48820326678768,
            # [603784.0 - 747064.0)
            msprime.PopulationParametersChange(
                time=603784.0, initial_size=181.48820326678768, population_id=2
            ),
            # HungryW, 747064.0 generations at 181.48820326678768,
            # [747064.0 - inf)
            msprime.PopulationParametersChange(
                time=747064.0, initial_size=181.48820326678768, population_id=2
            ),
            # JLW, 216.0 generations at 3268913.936033892,
            # [216.0 - 483.0)
            msprime.PopulationParametersChange(
                time=216.0, initial_size=3268913.936033892, population_id=3
            ),
            # JLW, 483.0 generations at 2247.1733368389237,
            # [483.0 - 814.0)
            msprime.PopulationParametersChange(
                time=483.0, initial_size=2247.1733368389237, population_id=3
            ),
            # JLW, 814.0 generations at 8541.71560357898,
            # [814.0 - 1223.0)
            msprime.PopulationParametersChange(
                time=814.0, initial_size=8541.71560357898, population_id=3
            ),
            # JLW, 1223.0 generations at 8787.384830337569,
            # [1223.0 - 1729.0)
            msprime.PopulationParametersChange(
                time=1223.0, initial_size=8787.384830337569, population_id=3
            ),
            # JLW, 1729.0 generations at 8428.611765499165,
            # [1729.0 - 2354.0)
            msprime.PopulationParametersChange(
                time=1729.0, initial_size=8428.611765499165, population_id=3
            ),
            # JLW, 2354.0 generations at 10075.312964408959,
            # [2354.0 - 3128.0)
            msprime.PopulationParametersChange(
                time=2354.0, initial_size=10075.312964408959, population_id=3
            ),
            # JLW, 3128.0 generations at 14386.936661511347,
            # [3128.0 - 4085.0)
            msprime.PopulationParametersChange(
                time=3128.0, initial_size=14386.936661511347, population_id=3
            ),
            # JLW, 4085.0 generations at 21387.469295614395,
            # [4085.0 - 5269.0)
            msprime.PopulationParametersChange(
                time=4085.0, initial_size=21387.469295614395, population_id=3
            ),
            # JLW, 5269.0 generations at 30301.698864746857,
            # [5269.0 - 6733.0)
            msprime.PopulationParametersChange(
                time=5269.0, initial_size=30301.698864746857, population_id=3
            ),
            # JLW, 6733.0 generations at 38644.804012876455,
            # [6733.0 - 8544.0)
            msprime.PopulationParametersChange(
                time=6733.0, initial_size=38644.804012876455, population_id=3
            ),
            # JLW, 8544.0 generations at 43800.000438,
            # [8544.0 - 10784.0)
            msprime.PopulationParametersChange(
                time=8544.0, initial_size=43800.000438, population_id=3
            ),
            # JLW, 10784.0 generations at 44666.287003227146,
            # [10784.0 - 13555.0)
            msprime.PopulationParametersChange(
                time=10784.0, initial_size=44666.287003227146, population_id=3
            ),
            # JLW, 13555.0 generations at 42603.590630618346,
            # [13555.0 - 16982.0)
            msprime.PopulationParametersChange(
                time=13555.0, initial_size=42603.590630618346, population_id=3
            ),
            # JLW, 16982.0 generations at 40082.65042517671,
            # [16982.0 - 21222.0)
            msprime.PopulationParametersChange(
                time=16982.0, initial_size=40082.65042517671, population_id=3
            ),
            # JLW, 21222.0 generations at 38658.84727049209,
            # [21222.0 - 26465.0)
            msprime.PopulationParametersChange(
                time=21222.0, initial_size=38658.84727049209, population_id=3
            ),
            # JLW, 26465.0 generations at 38853.35965000894,
            # [26465.0 - 32951.0)
            msprime.PopulationParametersChange(
                time=26465.0, initial_size=38853.35965000894, population_id=3
            ),
            # JLW, 32951.0 generations at 40625.139648917546,
            # [32951.0 - 40973.0)
            msprime.PopulationParametersChange(
                time=32951.0, initial_size=40625.139648917546, population_id=3
            ),
            # JLW, 40973.0 generations at 43120.1759303178,
            # [40973.0 - 50896.0)
            msprime.PopulationParametersChange(
                time=40973.0, initial_size=43120.1759303178, population_id=3
            ),
            # JLW, 50896.0 generations at 44570.72817427155,
            # [50896.0 - 63169.0)
            msprime.PopulationParametersChange(
                time=50896.0, initial_size=44570.72817427155, population_id=3
            ),
            # JLW, 63169.0 generations at 43081.72171792674,
            # [63169.0 - 78350.0)
            msprime.PopulationParametersChange(
                time=63169.0, initial_size=43081.72171792674, population_id=3
            ),
            # JLW, 78350.0 generations at 38245.597931678065,
            # [78350.0 - 97127.0)
            msprime.PopulationParametersChange(
                time=78350.0, initial_size=38245.597931678065, population_id=3
            ),
            # JLW, 97127.0 generations at 31833.447403186532,
            # [97127.0 - 120352.0)
            msprime.PopulationParametersChange(
                time=97127.0, initial_size=31833.447403186532, population_id=3
            ),
            # JLW, 120352.0 generations at 26336.06131035073,
            # [120352.0 - 149080.0)
            msprime.PopulationParametersChange(
                time=120352.0, initial_size=26336.06131035073, population_id=3
            ),
            # JLW, 149080.0 generations at 22917.724224149093,
            # [149080.0 - 184614.0)
            msprime.PopulationParametersChange(
                time=149080.0, initial_size=22917.724224149093, population_id=3
            ),
            # JLW, 184614.0 generations at 21028.482027482125,
            # [184614.0 - 228565.0)
            msprime.PopulationParametersChange(
                time=184614.0, initial_size=21028.482027482125, population_id=3
            ),
            # JLW, 228565.0 generations at 20039.31714022913,
            # [228565.0 - 282929.0)
            msprime.PopulationParametersChange(
                time=228565.0, initial_size=20039.31714022913, population_id=3
            ),
            # JLW, 282929.0 generations at 19999.20003199872,
            # [282929.0 - 350171.0)
            msprime.PopulationParametersChange(
                time=282929.0, initial_size=19999.20003199872, population_id=3
            ),
            # JLW, 350171.0 generations at 19999.20003199872,
            # [350171.0 - 433344.0)
            msprime.PopulationParametersChange(
                time=350171.0, initial_size=19999.20003199872, population_id=3
            ),
            # JLW, 433344.0 generations at 20855.427051365878,
            # [433344.0 - 536220.0)
            msprime.PopulationParametersChange(
                time=433344.0, initial_size=20855.427051365878, population_id=3
            ),
            # JLW, 536220.0 generations at 20855.427051365878,
            # [536220.0 - 663464.0)
            msprime.PopulationParametersChange(
                time=536220.0, initial_size=20855.427051365878, population_id=3
            ),
            # JLW, 663464.0 generations at 20855.427051365878,
            # [663464.0 - inf)
            msprime.PopulationParametersChange(
                time=663464.0, initial_size=20855.427051365878, population_id=3
            ),
            # JXW, 618.0 generations at 65893.73317650624,
            # [618.0 - 1383.0)
            msprime.PopulationParametersChange(
                time=618.0, initial_size=65893.73317650624, population_id=4
            ),
            # JXW, 1383.0 generations at 10934.697983641692,
            # [1383.0 - 2328.0)
            msprime.PopulationParametersChange(
                time=1383.0, initial_size=10934.697983641692, population_id=4
            ),
            # JXW, 2328.0 generations at 9776.317847645863,
            # [2328.0 - 3498.0)
            msprime.PopulationParametersChange(
                time=2328.0, initial_size=9776.317847645863, population_id=4
            ),
            # JXW, 3498.0 generations at 10284.309742840835,
            # [3498.0 - 4945.0)
            msprime.PopulationParametersChange(
                time=3498.0, initial_size=10284.309742840835, population_id=4
            ),
            # JXW, 4945.0 generations at 11635.212808042259,
            # [4945.0 - 6735.0)
            msprime.PopulationParametersChange(
                time=4945.0, initial_size=11635.212808042259, population_id=4
            ),
            # JXW, 6735.0 generations at 14079.54945441746,
            # [6735.0 - 8948.0)
            msprime.PopulationParametersChange(
                time=6735.0, initial_size=14079.54945441746, population_id=4
            ),
            # JXW, 8948.0 generations at 18153.103273004523,
            # [8948.0 - 11686.0)
            msprime.PopulationParametersChange(
                time=8948.0, initial_size=18153.103273004523, population_id=4
            ),
            # JXW, 11686.0 generations at 23430.754505441208,
            # [11686.0 - 15073.0)
            msprime.PopulationParametersChange(
                time=11686.0, initial_size=23430.754505441208, population_id=4
            ),
            # JXW, 15073.0 generations at 28981.601030585734,
            # [15073.0 - 19261.0)
            msprime.PopulationParametersChange(
                time=15073.0, initial_size=28981.601030585734, population_id=4
            ),
            # JXW, 19261.0 generations at 34042.49524681661,
            # [19261.0 - 24443.0)
            msprime.PopulationParametersChange(
                time=19261.0, initial_size=34042.49524681661, population_id=4
            ),
            # JXW, 24443.0 generations at 37705.75559506281,
            # [24443.0 - 30851.0)
            msprime.PopulationParametersChange(
                time=24443.0, initial_size=37705.75559506281, population_id=4
            ),
            # JXW, 30851.0 generations at 40586.392194425054,
            # [30851.0 - 38778.0)
            msprime.PopulationParametersChange(
                time=30851.0, initial_size=40586.392194425054, population_id=4
            ),
            # JXW, 38778.0 generations at 44889.44851068032,
            # [38778.0 - 48583.0)
            msprime.PopulationParametersChange(
                time=38778.0, initial_size=44889.44851068032, population_id=4
            ),
            # JXW, 48583.0 generations at 51770.28489187776,
            # [48583.0 - 60710.0)
            msprime.PopulationParametersChange(
                time=48583.0, initial_size=51770.28489187776, population_id=4
            ),
            # JXW, 60710.0 generations at 59890.10166344758,
            # [60710.0 - 75711.0)
            msprime.PopulationParametersChange(
                time=60710.0, initial_size=59890.10166344758, population_id=4
            ),
            # JXW, 75711.0 generations at 65818.91899007451,
            # [75711.0 - 94265.0)
            msprime.PopulationParametersChange(
                time=75711.0, initial_size=65818.91899007451, population_id=4
            ),
            # JXW, 94265.0 generations at 65086.157801389585,
            # [94265.0 - 117214.0)
            msprime.PopulationParametersChange(
                time=94265.0, initial_size=65086.157801389585, population_id=4
            ),
            # JXW, 117214.0 generations at 56381.38624914371,
            # [117214.0 - 145601.0)
            msprime.PopulationParametersChange(
                time=117214.0, initial_size=56381.38624914371, population_id=4
            ),
            # JXW, 145601.0 generations at 43681.8563041655,
            # [145601.0 - 180712.0)
            msprime.PopulationParametersChange(
                time=145601.0, initial_size=43681.8563041655, population_id=4
            ),
            # JXW, 180712.0 generations at 32467.32164076857,
            # [180712.0 - 224141.0)
            msprime.PopulationParametersChange(
                time=180712.0, initial_size=32467.32164076857, population_id=4
            ),
            # JXW, 224141.0 generations at 25639.809393656968,
            # [224141.0 - 277859.0)
            msprime.PopulationParametersChange(
                time=224141.0, initial_size=25639.809393656968, population_id=4
            ),
            # JXW, 277859.0 generations at 22371.86514239692,
            # [277859.0 - 344302.0)
            msprime.PopulationParametersChange(
                time=277859.0, initial_size=22371.86514239692, population_id=4
            ),
            # JXW, 344302.0 generations at 21600.06825621569,
            # [344302.0 - 426484.0)
            msprime.PopulationParametersChange(
                time=344302.0, initial_size=21600.06825621569, population_id=4
            ),
            # JXW, 426484.0 generations at 23261.007108563776,
            # [426484.0 - 528140.0)
            msprime.PopulationParametersChange(
                time=426484.0, initial_size=23261.007108563776, population_id=4
            ),
            # JXW, 528140.0 generations at 27237.825713324608,
            # [528140.0 - 653876.0)
            msprime.PopulationParametersChange(
                time=528140.0, initial_size=27237.825713324608, population_id=4
            ),
            # JXW, 653876.0 generations at 32287.12292676312,
            # [653876.0 - 809396.0)
            msprime.PopulationParametersChange(
                time=653876.0, initial_size=32287.12292676312, population_id=4
            ),
            # JXW, 809396.0 generations at 37983.96317074931,
            # [809396.0 - 1001760.0)
            msprime.PopulationParametersChange(
                time=809396.0, initial_size=37983.96317074931, population_id=4
            ),
            # JXW, 1001760.0 generations at 37983.96317074931,
            # [1001760.0 - 1239700.0)
            msprime.PopulationParametersChange(
                time=1001760.0, initial_size=37983.96317074931, population_id=4
            ),
            # JXW, 1239700.0 generations at 42030.93476798924,
            # [1239700.0 - 1534004.0)
            msprime.PopulationParametersChange(
                time=1239700.0, initial_size=42030.93476798924, population_id=4
            ),
            # JXW, 1534004.0 generations at 42030.93476798924,
            # [1534004.0 - 1898028.0)
            msprime.PopulationParametersChange(
                time=1534004.0, initial_size=42030.93476798924, population_id=4
            ),
            # JXW, 1898028.0 generations at 42030.93476798924,
            # [1898028.0 - inf)
            msprime.PopulationParametersChange(
                time=1898028.0, initial_size=42030.93476798924, population_id=4
            ),
            # HLJW, 631.0 generations at 24816.94401669684,
            # [631.0 - 1413.0)
            msprime.PopulationParametersChange(
                time=631.0, initial_size=24816.94401669684, population_id=5
            ),
            # HLJW, 1413.0 generations at 9042.327133311028,
            # [1413.0 - 2379.0)
            msprime.PopulationParametersChange(
                time=1413.0, initial_size=9042.327133311028, population_id=5
            ),
            # HLJW, 2379.0 generations at 9811.280028648938,
            # [2379.0 - 3574.0)
            msprime.PopulationParametersChange(
                time=2379.0, initial_size=9811.280028648938, population_id=5
            ),
            # HLJW, 3574.0 generations at 11282.861333634211,
            # [3574.0 - 5052.0)
            msprime.PopulationParametersChange(
                time=3574.0, initial_size=11282.861333634211, population_id=5
            ),
            # HLJW, 5052.0 generations at 13187.740676267344,
            # [5052.0 - 6880.0)
            msprime.PopulationParametersChange(
                time=5052.0, initial_size=13187.740676267344, population_id=5
            ),
            # HLJW, 6880.0 generations at 15720.922189295625,
            # [6880.0 - 9141.0)
            msprime.PopulationParametersChange(
                time=6880.0, initial_size=15720.922189295625, population_id=5
            ),
            # HLJW, 9141.0 generations at 19339.554223275154,
            # [9141.0 - 11938.0)
            msprime.PopulationParametersChange(
                time=9141.0, initial_size=19339.554223275154, population_id=5
            ),
            # HLJW, 11938.0 generations at 23784.295467545737,
            # [11938.0 - 15398.0)
            msprime.PopulationParametersChange(
                time=11938.0, initial_size=23784.295467545737, population_id=5
            ),
            # HLJW, 15398.0 generations at 28008.96847170464,
            # [15398.0 - 19677.0)
            msprime.PopulationParametersChange(
                time=15398.0, initial_size=28008.96847170464, population_id=5
            ),
            # HLJW, 19677.0 generations at 31718.894469652158,
            # [19677.0 - 24970.0)
            msprime.PopulationParametersChange(
                time=19677.0, initial_size=31718.894469652158, population_id=5
            ),
            # HLJW, 24970.0 generations at 35526.187240769854,
            # [24970.0 - 31517.0)
            msprime.PopulationParametersChange(
                time=24970.0, initial_size=35526.187240769854, population_id=5
            ),
            # HLJW, 31517.0 generations at 40272.40253071778,
            # [31517.0 - 39615.0)
            msprime.PopulationParametersChange(
                time=31517.0, initial_size=40272.40253071778, population_id=5
            ),
            # HLJW, 39615.0 generations at 46178.178501749004,
            # [39615.0 - 49631.0)
            msprime.PopulationParametersChange(
                time=39615.0, initial_size=46178.178501749004, population_id=5
            ),
            # HLJW, 49631.0 generations at 53306.32479543699,
            # [49631.0 - 62020.0)
            msprime.PopulationParametersChange(
                time=49631.0, initial_size=53306.32479543699, population_id=5
            ),
            # HLJW, 62020.0 generations at 61254.42946093039,
            # [62020.0 - 77344.0)
            msprime.PopulationParametersChange(
                time=62020.0, initial_size=61254.42946093039, population_id=5
            ),
            # HLJW, 77344.0 generations at 66642.00912329105,
            # [77344.0 - 96299.0)
            msprime.PopulationParametersChange(
                time=77344.0, initial_size=66642.00912329105, population_id=5
            ),
            # HLJW, 96299.0 generations at 64785.09165470842,
            # [96299.0 - 119744.0)
            msprime.PopulationParametersChange(
                time=96299.0, initial_size=64785.09165470842, population_id=5
            ),
            # HLJW, 119744.0 generations at 55031.98734264291,
            # [119744.0 - 148742.0)
            msprime.PopulationParametersChange(
                time=119744.0, initial_size=55031.98734264291, population_id=5
            ),
            # HLJW, 148742.0 generations at 41885.60622085024,
            # [148742.0 - 184611.0)
            msprime.PopulationParametersChange(
                time=148742.0, initial_size=41885.60622085024, population_id=5
            ),
            # HLJW, 184611.0 generations at 31061.68851338759,
            # [184611.0 - 228978.0)
            msprime.PopulationParametersChange(
                time=184611.0, initial_size=31061.68851338759, population_id=5
            ),
            # HLJW, 228978.0 generations at 24907.56181122808,
            # [228978.0 - 283854.0)
            msprime.PopulationParametersChange(
                time=228978.0, initial_size=24907.56181122808, population_id=5
            ),
            # HLJW, 283854.0 generations at 21939.615595995143,
            # [283854.0 - 351731.0)
            msprime.PopulationParametersChange(
                time=283854.0, initial_size=21939.615595995143, population_id=5
            ),
            # HLJW, 351731.0 generations at 21006.10857637401,
            # [351731.0 - 435688.0)
            msprime.PopulationParametersChange(
                time=351731.0, initial_size=21006.10857637401, population_id=5
            ),
            # HLJW, 435688.0 generations at 21762.738274780688,
            # [435688.0 - 539536.0)
            msprime.PopulationParametersChange(
                time=435688.0, initial_size=21762.738274780688, population_id=5
            ),
            # HLJW, 539536.0 generations at 24110.764853738074,
            # [539536.0 - 667984.0)
            msprime.PopulationParametersChange(
                time=539536.0, initial_size=24110.764853738074, population_id=5
            ),
            # HLJW, 667984.0 generations at 28036.17788394144,
            # [667984.0 - 826860.0)
            msprime.PopulationParametersChange(
                time=667984.0, initial_size=28036.17788394144, population_id=5
            ),
            # HLJW, 826860.0 generations at 34467.18895947004,
            # [826860.0 - 1023376.0)
            msprime.PopulationParametersChange(
                time=826860.0, initial_size=34467.18895947004, population_id=5
            ),
            # HLJW, 1023376.0 generations at 34467.18895947004,
            # [1023376.0 - 1266448.0)
            msprime.PopulationParametersChange(
                time=1023376.0, initial_size=34467.18895947004, population_id=5
            ),
            # HLJW, 1266448.0 generations at 38471.97028425015,
            # [1266448.0 - 1567104.0)
            msprime.PopulationParametersChange(
                time=1266448.0, initial_size=38471.97028425015, population_id=5
            ),
            # HLJW, 1567104.0 generations at 38471.97028425015,
            # [1567104.0 - 1938984.0)
            msprime.PopulationParametersChange(
                time=1567104.0, initial_size=38471.97028425015, population_id=5
            ),
            # HLJW, 1938984.0 generations at 38471.97028425015,
            # [1938984.0 - inf)
            msprime.PopulationParametersChange(
                time=1938984.0, initial_size=38471.97028425015, population_id=5
            ),
            # LNW, 576.0 generations at 42760.35174665347,
            # [576.0 - 1288.0)
            msprime.PopulationParametersChange(
                time=576.0, initial_size=42760.35174665347, population_id=6
            ),
            # LNW, 1288.0 generations at 10543.407207473167,
            # [1288.0 - 2169.0)
            msprime.PopulationParametersChange(
                time=1288.0, initial_size=10543.407207473167, population_id=6
            ),
            # LNW, 2169.0 generations at 9537.93475098837,
            # [2169.0 - 3259.0)
            msprime.PopulationParametersChange(
                time=2169.0, initial_size=9537.93475098837, population_id=6
            ),
            # LNW, 3259.0 generations at 10704.116803322557,
            # [3259.0 - 4607.0)
            msprime.PopulationParametersChange(
                time=3259.0, initial_size=10704.116803322557, population_id=6
            ),
            # LNW, 4607.0 generations at 12946.743570323475,
            # [4607.0 - 6274.0)
            msprime.PopulationParametersChange(
                time=4607.0, initial_size=12946.743570323475, population_id=6
            ),
            # LNW, 6274.0 generations at 16007.811812164335,
            # [6274.0 - 8337.0)
            msprime.PopulationParametersChange(
                time=6274.0, initial_size=16007.811812164335, population_id=6
            ),
            # LNW, 8337.0 generations at 21565.64420353224,
            # [8337.0 - 10887.0)
            msprime.PopulationParametersChange(
                time=8337.0, initial_size=21565.64420353224, population_id=6
            ),
            # LNW, 10887.0 generations at 29731.37700872616,
            # [10887.0 - 14043.0)
            msprime.PopulationParametersChange(
                time=10887.0, initial_size=29731.37700872616, population_id=6
            ),
            # LNW, 14043.0 generations at 37805.39672038184,
            # [14043.0 - 17945.0)
            msprime.PopulationParametersChange(
                time=14043.0, initial_size=37805.39672038184, population_id=6
            ),
            # LNW, 17945.0 generations at 42720.80246755355,
            # [17945.0 - 22772.0)
            msprime.PopulationParametersChange(
                time=17945.0, initial_size=42720.80246755355, population_id=6
            ),
            # LNW, 22772.0 generations at 45062.6032215255,
            # [22772.0 - 28743.0)
            msprime.PopulationParametersChange(
                time=22772.0, initial_size=45062.6032215255, population_id=6
            ),
            # LNW, 28743.0 generations at 47061.3703800441,
            # [28743.0 - 36128.0)
            msprime.PopulationParametersChange(
                time=28743.0, initial_size=47061.3703800441, population_id=6
            ),
            # LNW, 36128.0 generations at 50237.62396133713,
            # [36128.0 - 45263.0)
            msprime.PopulationParametersChange(
                time=36128.0, initial_size=50237.62396133713, population_id=6
            ),
            # LNW, 45263.0 generations at 56282.221571849885,
            # [45263.0 - 56561.0)
            msprime.PopulationParametersChange(
                time=45263.0, initial_size=56282.221571849885, population_id=6
            ),
            # LNW, 56561.0 generations at 65128.97164610219,
            # [56561.0 - 70537.0)
            msprime.PopulationParametersChange(
                time=56561.0, initial_size=65128.97164610219, population_id=6
            ),
            # LNW, 70537.0 generations at 73958.47970949109,
            # [70537.0 - 87823.0)
            msprime.PopulationParametersChange(
                time=70537.0, initial_size=73958.47970949109, population_id=6
            ),
            # LNW, 87823.0 generations at 76746.55983545538,
            # [87823.0 - 109204.0)
            msprime.PopulationParametersChange(
                time=87823.0, initial_size=76746.55983545538, population_id=6
            ),
            # LNW, 109204.0 generations at 68771.53408661087,
            # [109204.0 - 135650.0)
            msprime.PopulationParametersChange(
                time=109204.0, initial_size=68771.53408661087, population_id=6
            ),
            # LNW, 135650.0 generations at 52668.30816227105,
            # [135650.0 - 168362.0)
            msprime.PopulationParametersChange(
                time=135650.0, initial_size=52668.30816227105, population_id=6
            ),
            # LNW, 168362.0 generations at 36044.60880786061,
            # [168362.0 - 208823.0)
            msprime.PopulationParametersChange(
                time=168362.0, initial_size=36044.60880786061, population_id=6
            ),
            # LNW, 208823.0 generations at 25076.388949838445,
            # [208823.0 - 258870.0)
            msprime.PopulationParametersChange(
                time=208823.0, initial_size=25076.388949838445, population_id=6
            ),
            # LNW, 258870.0 generations at 19793.749134023477,
            # [258870.0 - 320772.0)
            msprime.PopulationParametersChange(
                time=258870.0, initial_size=19793.749134023477, population_id=6
            ),
            # LNW, 320772.0 generations at 17975.750712289126,
            # [320772.0 - 397340.0)
            msprime.PopulationParametersChange(
                time=320772.0, initial_size=17975.750712289126, population_id=6
            ),
            # LNW, 397340.0 generations at 18016.07033473859,
            # [397340.0 - 492044.0)
            msprime.PopulationParametersChange(
                time=397340.0, initial_size=18016.07033473859, population_id=6
            ),
            # LNW, 492044.0 generations at 19149.93441147464,
            # [492044.0 - 609188.0)
            msprime.PopulationParametersChange(
                time=492044.0, initial_size=19149.93441147464, population_id=6
            ),
            # LNW, 609188.0 generations at 21743.857360295722,
            # [609188.0 - 754080.0)
            msprime.PopulationParametersChange(
                time=609188.0, initial_size=21743.857360295722, population_id=6
            ),
            # LNW, 754080.0 generations at 28649.0675444741,
            # [754080.0 - 933300.0)
            msprime.PopulationParametersChange(
                time=754080.0, initial_size=28649.0675444741, population_id=6
            ),
            # LNW, 933300.0 generations at 28649.0675444741,
            # [933300.0 - 1154976.0)
            msprime.PopulationParametersChange(
                time=933300.0, initial_size=28649.0675444741, population_id=6
            ),
            # LNW, 1154976.0 generations at 44198.89502762431,
            # [1154976.0 - 1429168.0)
            msprime.PopulationParametersChange(
                time=1154976.0, initial_size=44198.89502762431, population_id=6
            ),
            # LNW, 1429168.0 generations at 44198.89502762431,
            # [1429168.0 - 1768316.0)
            msprime.PopulationParametersChange(
                time=1429168.0, initial_size=44198.89502762431, population_id=6
            ),
            # LNW, 1768316.0 generations at 44198.89502762431,
            # [1768316.0 - inf)
            msprime.PopulationParametersChange(
                time=1768316.0, initial_size=44198.89502762431, population_id=6
            ),
            # MW, 513.0 generations at 475403.7960993119,
            # [513.0 - 1147.0)
            msprime.PopulationParametersChange(
                time=513.0, initial_size=475403.7960993119, population_id=7
            ),
            # MW, 1147.0 generations at 17332.3742752901,
            # [1147.0 - 1931.0)
            msprime.PopulationParametersChange(
                time=1147.0, initial_size=17332.3742752901, population_id=7
            ),
            # MW, 1931.0 generations at 10851.16541516559,
            # [1931.0 - 2901.0)
            msprime.PopulationParametersChange(
                time=1931.0, initial_size=10851.16541516559, population_id=7
            ),
            # MW, 2901.0 generations at 12282.220871177928,
            # [2901.0 - 4101.0)
            msprime.PopulationParametersChange(
                time=2901.0, initial_size=12282.220871177928, population_id=7
            ),
            # MW, 4101.0 generations at 17598.817359473444,
            # [4101.0 - 5586.0)
            msprime.PopulationParametersChange(
                time=4101.0, initial_size=17598.817359473444, population_id=7
            ),
            # MW, 5586.0 generations at 28472.588015887708,
            # [5586.0 - 7422.0)
            msprime.PopulationParametersChange(
                time=5586.0, initial_size=28472.588015887708, population_id=7
            ),
            # MW, 7422.0 generations at 42425.28377211684,
            # [7422.0 - 9693.0)
            msprime.PopulationParametersChange(
                time=7422.0, initial_size=42425.28377211684, population_id=7
            ),
            # MW, 9693.0 generations at 53361.081095502996,
            # [9693.0 - 12502.0)
            msprime.PopulationParametersChange(
                time=9693.0, initial_size=53361.081095502996, population_id=7
            ),
            # MW, 12502.0 generations at 59229.24977270776,
            # [12502.0 - 15976.0)
            msprime.PopulationParametersChange(
                time=12502.0, initial_size=59229.24977270776, population_id=7
            ),
            # MW, 15976.0 generations at 62775.03311382997,
            # [15976.0 - 20273.0)
            msprime.PopulationParametersChange(
                time=15976.0, initial_size=62775.03311382997, population_id=7
            ),
            # MW, 20273.0 generations at 61628.34410802216,
            # [20273.0 - 25589.0)
            msprime.PopulationParametersChange(
                time=20273.0, initial_size=61628.34410802216, population_id=7
            ),
            # MW, 25589.0 generations at 57202.34758434486,
            # [25589.0 - 32163.0)
            msprime.PopulationParametersChange(
                time=25589.0, initial_size=57202.34758434486, population_id=7
            ),
            # MW, 32163.0 generations at 55381.565138412385,
            # [32163.0 - 40296.0)
            msprime.PopulationParametersChange(
                time=32163.0, initial_size=55381.565138412385, population_id=7
            ),
            # MW, 40296.0 generations at 58253.33209059559,
            # [40296.0 - 50354.0)
            msprime.PopulationParametersChange(
                time=40296.0, initial_size=58253.33209059559, population_id=7
            ),
            # MW, 50354.0 generations at 66091.66914510426,
            # [50354.0 - 62796.0)
            msprime.PopulationParametersChange(
                time=50354.0, initial_size=66091.66914510426, population_id=7
            ),
            # MW, 62796.0 generations at 76658.017087072,
            # [62796.0 - 78185.0)
            msprime.PopulationParametersChange(
                time=62796.0, initial_size=76658.017087072, population_id=7
            ),
            # MW, 78185.0 generations at 83892.61744966444,
            # [78185.0 - 97220.0)
            msprime.PopulationParametersChange(
                time=78185.0, initial_size=83892.61744966444, population_id=7
            ),
            # MW, 97220.0 generations at 81435.54580138685,
            # [97220.0 - 120765.0)
            msprime.PopulationParametersChange(
                time=97220.0, initial_size=81435.54580138685, population_id=7
            ),
            # MW, 120765.0 generations at 67792.24388937661,
            # [120765.0 - 149887.0)
            msprime.PopulationParametersChange(
                time=120765.0, initial_size=67792.24388937661, population_id=7
            ),
            # MW, 149887.0 generations at 48709.44329977252,
            # [149887.0 - 185908.0)
            msprime.PopulationParametersChange(
                time=149887.0, initial_size=48709.44329977252, population_id=7
            ),
            # MW, 185908.0 generations at 32856.04452651154,
            # [185908.0 - 230462.0)
            msprime.PopulationParametersChange(
                time=185908.0, initial_size=32856.04452651154, population_id=7
            ),
            # MW, 230462.0 generations at 24603.03010918825,
            # [230462.0 - 285572.0)
            msprime.PopulationParametersChange(
                time=230462.0, initial_size=24603.03010918825, population_id=7
            ),
            # MW, 285572.0 generations at 21448.411208939695,
            # [285572.0 - 353737.0)
            msprime.PopulationParametersChange(
                time=285572.0, initial_size=21448.411208939695, population_id=7
            ),
            # MW, 353737.0 generations at 21353.445431804033,
            # [353737.0 - 438052.0)
            msprime.PopulationParametersChange(
                time=353737.0, initial_size=21353.445431804033, population_id=7
            ),
            # MW, 438052.0 generations at 23923.67391075513,
            # [438052.0 - 542336.0)
            msprime.PopulationParametersChange(
                time=438052.0, initial_size=23923.67391075513, population_id=7
            ),
            # MW, 542336.0 generations at 28908.501701265326,
            # [542336.0 - 671332.0)
            msprime.PopulationParametersChange(
                time=542336.0, initial_size=28908.501701265326, population_id=7
            ),
            # MW, 671332.0 generations at 35324.08077910793,
            # [671332.0 - 830884.0)
            msprime.PopulationParametersChange(
                time=671332.0, initial_size=35324.08077910793, population_id=7
            ),
            # MW, 830884.0 generations at 35324.08077910793,
            # [830884.0 - 1028232.0)
            msprime.PopulationParametersChange(
                time=830884.0, initial_size=35324.08077910793, population_id=7
            ),
            # MW, 1028232.0 generations at 33548.32216453775,
            # [1028232.0 - 1272336.0)
            msprime.PopulationParametersChange(
                time=1028232.0, initial_size=33548.32216453775, population_id=7
            ),
            # MW, 1272336.0 generations at 33548.32216453775,
            # [1272336.0 - 1574268.0)
            msprime.PopulationParametersChange(
                time=1272336.0, initial_size=33548.32216453775, population_id=7
            ),
            # MW, 1574268.0 generations at 33548.32216453775,
            # [1574268.0 - inf)
            msprime.PopulationParametersChange(
                time=1574268.0, initial_size=33548.32216453775, population_id=7
            ),
            # NetherlandW, 187.0 generations at 14409.948628533139,
            # [187.0 - 419.0)
            msprime.PopulationParametersChange(
                time=187.0, initial_size=14409.948628533139, population_id=8
            ),
            # NetherlandW, 419.0 generations at 3577.3951554914806,
            # [419.0 - 705.0)
            msprime.PopulationParametersChange(
                time=419.0, initial_size=3577.3951554914806, population_id=8
            ),
            # NetherlandW, 705.0 generations at 2638.508503912908,
            # [705.0 - 1059.0)
            msprime.PopulationParametersChange(
                time=705.0, initial_size=2638.508503912908, population_id=8
            ),
            # NetherlandW, 1059.0 generations at 2688.941861043552,
            # [1059.0 - 1497.0)
            msprime.PopulationParametersChange(
                time=1059.0, initial_size=2688.941861043552, population_id=8
            ),
            # NetherlandW, 1497.0 generations at 3014.381614683656,
            # [1497.0 - 2039.0)
            msprime.PopulationParametersChange(
                time=1497.0, initial_size=3014.381614683656, population_id=8
            ),
            # NetherlandW, 2039.0 generations at 3851.4792568955922,
            # [2039.0 - 2709.0)
            msprime.PopulationParametersChange(
                time=2039.0, initial_size=3851.4792568955922, population_id=8
            ),
            # NetherlandW, 2709.0 generations at 4771.9713490840195,
            # [2709.0 - 3538.0)
            msprime.PopulationParametersChange(
                time=2709.0, initial_size=4771.9713490840195, population_id=8
            ),
            # NetherlandW, 3538.0 generations at 5151.61193937583,
            # [3538.0 - 4563.0)
            msprime.PopulationParametersChange(
                time=3538.0, initial_size=5151.61193937583, population_id=8
            ),
            # NetherlandW, 4563.0 generations at 5124.499720714765,
            # [4563.0 - 5831.0)
            msprime.PopulationParametersChange(
                time=4563.0, initial_size=5124.499720714765, population_id=8
            ),
            # NetherlandW, 5831.0 generations at 5386.799110100787,
            # [5831.0 - 7399.0)
            msprime.PopulationParametersChange(
                time=5831.0, initial_size=5386.799110100787, population_id=8
            ),
            # NetherlandW, 7399.0 generations at 6947.677044180279,
            # [7399.0 - 9339.0)
            msprime.PopulationParametersChange(
                time=7399.0, initial_size=6947.677044180279, population_id=8
            ),
            # NetherlandW, 9339.0 generations at 11343.591721446763,
            # [9339.0 - 11739.0)
            msprime.PopulationParametersChange(
                time=9339.0, initial_size=11343.591721446763, population_id=8
            ),
            # NetherlandW, 11739.0 generations at 20032.994341680747,
            # [11739.0 - 14707.0)
            msprime.PopulationParametersChange(
                time=11739.0, initial_size=20032.994341680747, population_id=8
            ),
            # NetherlandW, 14707.0 generations at 32805.92343753589,
            # [14707.0 - 18378.0)
            msprime.PopulationParametersChange(
                time=14707.0, initial_size=32805.92343753589, population_id=8
            ),
            # NetherlandW, 18378.0 generations at 47107.03896929799,
            # [18378.0 - 22919.0)
            msprime.PopulationParametersChange(
                time=18378.0, initial_size=47107.03896929799, population_id=8
            ),
            # NetherlandW, 22919.0 generations at 58610.11962325415,
            # [22919.0 - 28536.0)
            msprime.PopulationParametersChange(
                time=22919.0, initial_size=58610.11962325415, population_id=8
            ),
            # NetherlandW, 28536.0 generations at 66264.6610562587,
            # [28536.0 - 35483.0)
            msprime.PopulationParametersChange(
                time=28536.0, initial_size=66264.6610562587, population_id=8
            ),
            # NetherlandW, 35483.0 generations at 71951.50468584175,
            # [35483.0 - 44076.0)
            msprime.PopulationParametersChange(
                time=35483.0, initial_size=71951.50468584175, population_id=8
            ),
            # NetherlandW, 44076.0 generations at 78095.10421791658,
            # [44076.0 - 54705.0)
            msprime.PopulationParametersChange(
                time=44076.0, initial_size=78095.10421791658, population_id=8
            ),
            # NetherlandW, 54705.0 generations at 85281.66401582828,
            # [54705.0 - 67852.0)
            msprime.PopulationParametersChange(
                time=54705.0, initial_size=85281.66401582828, population_id=8
            ),
            # NetherlandW, 67852.0 generations at 92375.76614151044,
            # [67852.0 - 84114.0)
            msprime.PopulationParametersChange(
                time=67852.0, initial_size=92375.76614151044, population_id=8
            ),
            # NetherlandW, 84114.0 generations at 94208.08682217283,
            # [84114.0 - 104227.0)
            msprime.PopulationParametersChange(
                time=84114.0, initial_size=94208.08682217283, population_id=8
            ),
            # NetherlandW, 104227.0 generations at 84847.14786312458,
            # [104227.0 - 129106.0)
            msprime.PopulationParametersChange(
                time=104227.0, initial_size=84847.14786312458, population_id=8
            ),
            # NetherlandW, 129106.0 generations at 66174.3296540406,
            # [129106.0 - 159878.0)
            msprime.PopulationParametersChange(
                time=129106.0, initial_size=66174.3296540406, population_id=8
            ),
            # NetherlandW, 159878.0 generations at 45770.98943147854,
            # [159878.0 - 197941.0)
            msprime.PopulationParametersChange(
                time=159878.0, initial_size=45770.98943147854, population_id=8
            ),
            # NetherlandW, 197941.0 generations at 30535.19028003823,
            # [197941.0 - 245021.0)
            msprime.PopulationParametersChange(
                time=197941.0, initial_size=30535.19028003823, population_id=8
            ),
            # NetherlandW, 245021.0 generations at 22397.26932492391,
            # [245021.0 - 303254.0)
            msprime.PopulationParametersChange(
                time=245021.0, initial_size=22397.26932492391, population_id=8
            ),
            # NetherlandW, 303254.0 generations at 22397.26932492391,
            # [303254.0 - 375282.0)
            msprime.PopulationParametersChange(
                time=303254.0, initial_size=22397.26932492391, population_id=8
            ),
            # NetherlandW, 375282.0 generations at 23136.9813411814,
            # [375282.0 - 464376.0)
            msprime.PopulationParametersChange(
                time=375282.0, initial_size=23136.9813411814, population_id=8
            ),
            # NetherlandW, 464376.0 generations at 23136.9813411814,
            # [464376.0 - 574572.0)
            msprime.PopulationParametersChange(
                time=464376.0, initial_size=23136.9813411814, population_id=8
            ),
            # NetherlandW, 574572.0 generations at 23136.9813411814,
            # [574572.0 - inf)
            msprime.PopulationParametersChange(
                time=574572.0, initial_size=23136.9813411814, population_id=8
            ),
            # SXW, 628.0 generations at 56807.36905191342,
            # [628.0 - 1404.0)
            msprime.PopulationParametersChange(
                time=628.0, initial_size=56807.36905191342, population_id=9
            ),
            # SXW, 1404.0 generations at 11527.44396221304,
            # [1404.0 - 2364.0)
            msprime.PopulationParametersChange(
                time=1404.0, initial_size=11527.44396221304, population_id=9
            ),
            # SXW, 2364.0 generations at 9426.669934578911,
            # [2364.0 - 3552.0)
            msprime.PopulationParametersChange(
                time=2364.0, initial_size=9426.669934578911, population_id=9
            ),
            # SXW, 3552.0 generations at 10144.765808081322,
            # [3552.0 - 5021.0)
            msprime.PopulationParametersChange(
                time=3552.0, initial_size=10144.765808081322, population_id=9
            ),
            # SXW, 5021.0 generations at 11731.30616362826,
            # [5021.0 - 6839.0)
            msprime.PopulationParametersChange(
                time=5021.0, initial_size=11731.30616362826, population_id=9
            ),
            # SXW, 6839.0 generations at 14329.213684399068,
            # [6839.0 - 9086.0)
            msprime.PopulationParametersChange(
                time=6839.0, initial_size=14329.213684399068, population_id=9
            ),
            # SXW, 9086.0 generations at 18143.058012427995,
            # [9086.0 - 11867.0)
            msprime.PopulationParametersChange(
                time=9086.0, initial_size=18143.058012427995, population_id=9
            ),
            # SXW, 11867.0 generations at 23561.03899469759,
            # [11867.0 - 15305.0)
            msprime.PopulationParametersChange(
                time=11867.0, initial_size=23561.03899469759, population_id=9
            ),
            # SXW, 15305.0 generations at 30528.292094698765,
            # [15305.0 - 19559.0)
            msprime.PopulationParametersChange(
                time=15305.0, initial_size=30528.292094698765, population_id=9
            ),
            # SXW, 19559.0 generations at 37101.758994393924,
            # [19559.0 - 24820.0)
            msprime.PopulationParametersChange(
                time=19559.0, initial_size=37101.758994393924, population_id=9
            ),
            # SXW, 24820.0 generations at 41195.91748457728,
            # [24820.0 - 31328.0)
            msprime.PopulationParametersChange(
                time=24820.0, initial_size=41195.91748457728, population_id=9
            ),
            # SXW, 31328.0 generations at 43624.30746411901,
            # [31328.0 - 39377.0)
            msprime.PopulationParametersChange(
                time=31328.0, initial_size=43624.30746411901, population_id=9
            ),
            # SXW, 39377.0 generations at 47277.41203446523,
            # [39377.0 - 49333.0)
            msprime.PopulationParametersChange(
                time=39377.0, initial_size=47277.41203446523, population_id=9
            ),
            # SXW, 49333.0 generations at 53750.2922672142,
            # [49333.0 - 61648.0)
            msprime.PopulationParametersChange(
                time=49333.0, initial_size=53750.2922672142, population_id=9
            ),
            # SXW, 61648.0 generations at 62296.95087573939,
            # [61648.0 - 76880.0)
            msprime.PopulationParametersChange(
                time=61648.0, initial_size=62296.95087573939, population_id=9
            ),
            # SXW, 76880.0 generations at 69008.35001035125,
            # [76880.0 - 95721.0)
            msprime.PopulationParametersChange(
                time=76880.0, initial_size=69008.35001035125, population_id=9
            ),
            # SXW, 95721.0 generations at 68371.15967742486,
            # [95721.0 - 119025.0)
            msprime.PopulationParametersChange(
                time=95721.0, initial_size=68371.15967742486, population_id=9
            ),
            # SXW, 119025.0 generations at 58981.73925352711,
            # [119025.0 - 147850.0)
            msprime.PopulationParametersChange(
                time=119025.0, initial_size=58981.73925352711, population_id=9
            ),
            # SXW, 147850.0 generations at 45438.74516361357,
            # [147850.0 - 183503.0)
            msprime.PopulationParametersChange(
                time=147850.0, initial_size=45438.74516361357, population_id=9
            ),
            # SXW, 183503.0 generations at 33828.81267633269,
            # [183503.0 - 227603.0)
            msprime.PopulationParametersChange(
                time=183503.0, initial_size=33828.81267633269, population_id=9
            ),
            # SXW, 227603.0 generations at 27119.416281683953,
            # [227603.0 - 282150.0)
            msprime.PopulationParametersChange(
                time=227603.0, initial_size=27119.416281683953, population_id=9
            ),
            # SXW, 282150.0 generations at 24178.826601545028,
            # [282150.0 - 349620.0)
            msprime.PopulationParametersChange(
                time=282150.0, initial_size=24178.826601545028, population_id=9
            ),
            # SXW, 349620.0 generations at 24094.55667822804,
            # [349620.0 - 433072.0)
            msprime.PopulationParametersChange(
                time=349620.0, initial_size=24094.55667822804, population_id=9
            ),
            # SXW, 433072.0 generations at 26648.65217779448,
            # [433072.0 - 536296.0)
            msprime.PopulationParametersChange(
                time=433072.0, initial_size=26648.65217779448, population_id=9
            ),
            # SXW, 536296.0 generations at 31281.672693602275,
            # [536296.0 - 663972.0)
            msprime.PopulationParametersChange(
                time=536296.0, initial_size=31281.672693602275, population_id=9
            ),
            # SXW, 663972.0 generations at 37780.258681431194,
            # [663972.0 - 821896.0)
            msprime.PopulationParametersChange(
                time=663972.0, initial_size=37780.258681431194, population_id=9
            ),
            # SXW, 821896.0 generations at 48919.012129469054,
            # [821896.0 - 1017236.0)
            msprime.PopulationParametersChange(
                time=821896.0, initial_size=48919.012129469054, population_id=9
            ),
            # SXW, 1017236.0 generations at 48919.012129469054,
            # [1017236.0 - 1258848.0)
            msprime.PopulationParametersChange(
                time=1017236.0, initial_size=48919.012129469054, population_id=9
            ),
            # SXW, 1258848.0 generations at 71862.57002109165,
            # [1258848.0 - 1557696.0)
            msprime.PopulationParametersChange(
                time=1258848.0, initial_size=71862.57002109165, population_id=9
            ),
            # SXW, 1557696.0 generations at 71862.57002109165,
            # [1557696.0 - 1927344.0)
            msprime.PopulationParametersChange(
                time=1557696.0, initial_size=71862.57002109165, population_id=9
            ),
            # SXW, 1927344.0 generations at 71862.57002109165,
            # [1927344.0 - inf)
            msprime.PopulationParametersChange(
                time=1927344.0, initial_size=71862.57002109165, population_id=9
            ),
            # XJW, 414.0 generations at 32182.901867895627,
            # [414.0 - 927.0)
            msprime.PopulationParametersChange(
                time=414.0, initial_size=32182.901867895627, population_id=10
            ),
            # XJW, 927.0 generations at 7911.830560236722,
            # [927.0 - 1561.0)
            msprime.PopulationParametersChange(
                time=927.0, initial_size=7911.830560236722, population_id=10
            ),
            # XJW, 1561.0 generations at 6585.055874199093,
            # [1561.0 - 2345.0)
            msprime.PopulationParametersChange(
                time=1561.0, initial_size=6585.055874199093, population_id=10
            ),
            # XJW, 2345.0 generations at 6870.774473698675,
            # [2345.0 - 3315.0)
            msprime.PopulationParametersChange(
                time=2345.0, initial_size=6870.774473698675, population_id=10
            ),
            # XJW, 3315.0 generations at 8063.833304437931,
            # [3315.0 - 4515.0)
            msprime.PopulationParametersChange(
                time=3315.0, initial_size=8063.833304437931, population_id=10
            ),
            # XJW, 4515.0 generations at 9721.243347024085,
            # [4515.0 - 5999.0)
            msprime.PopulationParametersChange(
                time=4515.0, initial_size=9721.243347024085, population_id=10
            ),
            # XJW, 5999.0 generations at 10870.15598673841,
            # [5999.0 - 7835.0)
            msprime.PopulationParametersChange(
                time=5999.0, initial_size=10870.15598673841, population_id=10
            ),
            # XJW, 7835.0 generations at 13108.994736738614,
            # [7835.0 - 10105.0)
            msprime.PopulationParametersChange(
                time=7835.0, initial_size=13108.994736738614, population_id=10
            ),
            # XJW, 10105.0 generations at 19161.30948389013,
            # [10105.0 - 12913.0)
            msprime.PopulationParametersChange(
                time=10105.0, initial_size=19161.30948389013, population_id=10
            ),
            # XJW, 12913.0 generations at 29815.679469519433,
            # [12913.0 - 16387.0)
            msprime.PopulationParametersChange(
                time=12913.0, initial_size=29815.679469519433, population_id=10
            ),
            # XJW, 16387.0 generations at 39991.20193557418,
            # [16387.0 - 20683.0)
            msprime.PopulationParametersChange(
                time=16387.0, initial_size=39991.20193557418, population_id=10
            ),
            # XJW, 20683.0 generations at 44415.722277371744,
            # [20683.0 - 25998.0)
            msprime.PopulationParametersChange(
                time=20683.0, initial_size=44415.722277371744, population_id=10
            ),
            # XJW, 25998.0 generations at 45482.144846986695,
            # [25998.0 - 32571.0)
            msprime.PopulationParametersChange(
                time=25998.0, initial_size=45482.144846986695, population_id=10
            ),
            # XJW, 32571.0 generations at 48548.87415160842,
            # [32571.0 - 40701.0)
            msprime.PopulationParametersChange(
                time=32571.0, initial_size=48548.87415160842, population_id=10
            ),
            # XJW, 40701.0 generations at 55380.184969817805,
            # [40701.0 - 50758.0)
            msprime.PopulationParametersChange(
                time=40701.0, initial_size=55380.184969817805, population_id=10
            ),
            # XJW, 50758.0 generations at 64930.42704741869,
            # [50758.0 - 63197.0)
            msprime.PopulationParametersChange(
                time=50758.0, initial_size=64930.42704741869, population_id=10
            ),
            # XJW, 63197.0 generations at 74475.1364756876,
            # [63197.0 - 78583.0)
            msprime.PopulationParametersChange(
                time=63197.0, initial_size=74475.1364756876, population_id=10
            ),
            # XJW, 78583.0 generations at 78971.78732897673,
            # [78583.0 - 97614.0)
            msprime.PopulationParametersChange(
                time=78583.0, initial_size=78971.78732897673, population_id=10
            ),
            # XJW, 97614.0 generations at 73258.73152506365,
            # [97614.0 - 121153.0)
            msprime.PopulationParametersChange(
                time=97614.0, initial_size=73258.73152506365, population_id=10
            ),
            # XJW, 121153.0 generations at 57620.94636642313,
            # [121153.0 - 150269.0)
            msprime.PopulationParametersChange(
                time=121153.0, initial_size=57620.94636642313, population_id=10
            ),
            # XJW, 150269.0 generations at 39265.653743685594,
            # [150269.0 - 186282.0)
            msprime.PopulationParametersChange(
                time=150269.0, initial_size=39265.653743685594, population_id=10
            ),
            # XJW, 186282.0 generations at 25963.640517818847,
            # [186282.0 - 230827.0)
            msprime.PopulationParametersChange(
                time=186282.0, initial_size=25963.640517818847, population_id=10
            ),
            # XJW, 230827.0 generations at 19153.05203884239,
            # [230827.0 - 285924.0)
            msprime.PopulationParametersChange(
                time=230827.0, initial_size=19153.05203884239, population_id=10
            ),
            # XJW, 285924.0 generations at 15756.840438355302,
            # [285924.0 - 354075.0)
            msprime.PopulationParametersChange(
                time=285924.0, initial_size=15756.840438355302, population_id=10
            ),
            # XJW, 354075.0 generations at 14418.88296913638,
            # [354075.0 - 438372.0)
            msprime.PopulationParametersChange(
                time=354075.0, initial_size=14418.88296913638, population_id=10
            ),
            # XJW, 438372.0 generations at 14936.519790888724,
            # [438372.0 - 542636.0)
            msprime.PopulationParametersChange(
                time=438372.0, initial_size=14936.519790888724, population_id=10
            ),
            # XJW, 542636.0 generations at 18123.000806473538,
            # [542636.0 - 671600.0)
            msprime.PopulationParametersChange(
                time=542636.0, initial_size=18123.000806473538, population_id=10
            ),
            # XJW, 671600.0 generations at 18123.000806473538,
            # [671600.0 - 831120.0)
            msprime.PopulationParametersChange(
                time=671600.0, initial_size=18123.000806473538, population_id=10
            ),
            # XJW, 831120.0 generations at 23769.03156646237,
            # [831120.0 - 1028428.0)
            msprime.PopulationParametersChange(
                time=831120.0, initial_size=23769.03156646237, population_id=10
            ),
            # XJW, 1028428.0 generations at 23769.03156646237,
            # [1028428.0 - 1272476.0)
            msprime.PopulationParametersChange(
                time=1028428.0, initial_size=23769.03156646237, population_id=10
            ),
            # XJW, 1272476.0 generations at 23769.03156646237,
            # [1272476.0 - inf)
            msprime.PopulationParametersChange(
                time=1272476.0, initial_size=23769.03156646237, population_id=10
            ),
            # YNW, 621.0 generations at 65330.22796983051,
            # [621.0 - 1388.0)
            msprime.PopulationParametersChange(
                time=621.0, initial_size=65330.22796983051, population_id=11
            ),
            # YNW, 1388.0 generations at 10536.131027325457,
            # [1388.0 - 2338.0)
            msprime.PopulationParametersChange(
                time=1388.0, initial_size=10536.131027325457, population_id=11
            ),
            # YNW, 2338.0 generations at 9339.553475948316,
            # [2338.0 - 3513.0)
            msprime.PopulationParametersChange(
                time=2338.0, initial_size=9339.553475948316, population_id=11
            ),
            # YNW, 3513.0 generations at 10421.38885849317,
            # [3513.0 - 4965.0)
            msprime.PopulationParametersChange(
                time=3513.0, initial_size=10421.38885849317, population_id=11
            ),
            # YNW, 4965.0 generations at 12366.365957867793,
            # [4965.0 - 6762.0)
            msprime.PopulationParametersChange(
                time=4965.0, initial_size=12366.365957867793, population_id=11
            ),
            # YNW, 6762.0 generations at 14871.76817887763,
            # [6762.0 - 8985.0)
            msprime.PopulationParametersChange(
                time=6762.0, initial_size=14871.76817887763, population_id=11
            ),
            # YNW, 8985.0 generations at 18485.484273474252,
            # [8985.0 - 11734.0)
            msprime.PopulationParametersChange(
                time=8985.0, initial_size=18485.484273474252, population_id=11
            ),
            # YNW, 11734.0 generations at 23740.02770461233,
            # [11734.0 - 15135.0)
            msprime.PopulationParametersChange(
                time=11734.0, initial_size=23740.02770461233, population_id=11
            ),
            # YNW, 15135.0 generations at 29474.39780120993,
            # [15135.0 - 19341.0)
            msprime.PopulationParametersChange(
                time=15135.0, initial_size=29474.39780120993, population_id=11
            ),
            # YNW, 19341.0 generations at 34178.74458053282,
            # [19341.0 - 24544.0)
            msprime.PopulationParametersChange(
                time=19341.0, initial_size=34178.74458053282, population_id=11
            ),
            # YNW, 24544.0 generations at 37964.20734531484,
            # [24544.0 - 30979.0)
            msprime.PopulationParametersChange(
                time=24544.0, initial_size=37964.20734531484, population_id=11
            ),
            # YNW, 30979.0 generations at 41980.558803218235,
            # [30979.0 - 38939.0)
            msprime.PopulationParametersChange(
                time=30979.0, initial_size=41980.558803218235, population_id=11
            ),
            # YNW, 38939.0 generations at 47134.57155852818,
            # [38939.0 - 48784.0)
            msprime.PopulationParametersChange(
                time=38939.0, initial_size=47134.57155852818, population_id=11
            ),
            # YNW, 48784.0 generations at 53964.64775925291,
            # [48784.0 - 60962.0)
            msprime.PopulationParametersChange(
                time=48784.0, initial_size=53964.64775925291, population_id=11
            ),
            # YNW, 60962.0 generations at 61821.12675185618,
            # [60962.0 - 76024.0)
            msprime.PopulationParametersChange(
                time=60962.0, initial_size=61821.12675185618, population_id=11
            ),
            # YNW, 76024.0 generations at 67232.99257411597,
            # [76024.0 - 94656.0)
            msprime.PopulationParametersChange(
                time=76024.0, initial_size=67232.99257411597, population_id=11
            ),
            # YNW, 94656.0 generations at 65379.98849312203,
            # [94656.0 - 117700.0)
            msprime.PopulationParametersChange(
                time=94656.0, initial_size=65379.98849312203, population_id=11
            ),
            # YNW, 117700.0 generations at 55363.62831074498,
            # [117700.0 - 146204.0)
            msprime.PopulationParametersChange(
                time=117700.0, initial_size=55363.62831074498, population_id=11
            ),
            # YNW, 146204.0 generations at 41940.23935294599,
            # [146204.0 - 181461.0)
            msprime.PopulationParametersChange(
                time=146204.0, initial_size=41940.23935294599, population_id=11
            ),
            # YNW, 181461.0 generations at 30985.797659642707,
            # [181461.0 - 225070.0)
            msprime.PopulationParametersChange(
                time=181461.0, initial_size=30985.797659642707, population_id=11
            ),
            # YNW, 225070.0 generations at 24828.37386565367,
            # [225070.0 - 279011.0)
            msprime.PopulationParametersChange(
                time=225070.0, initial_size=24828.37386565367, population_id=11
            ),
            # YNW, 279011.0 generations at 22164.152143605977,
            # [279011.0 - 345729.0)
            msprime.PopulationParametersChange(
                time=279011.0, initial_size=22164.152143605977, population_id=11
            ),
            # YNW, 345729.0 generations at 22088.58627525692,
            # [345729.0 - 428252.0)
            msprime.PopulationParametersChange(
                time=345729.0, initial_size=22088.58627525692, population_id=11
            ),
            # YNW, 428252.0 generations at 24609.841721802968,
            # [428252.0 - 530328.0)
            msprime.PopulationParametersChange(
                time=428252.0, initial_size=24609.841721802968, population_id=11
            ),
            # YNW, 530328.0 generations at 30050.876133293666,
            # [530328.0 - 656584.0)
            msprime.PopulationParametersChange(
                time=530328.0, initial_size=30050.876133293666, population_id=11
            ),
            # YNW, 656584.0 generations at 37758.432874002,
            # [656584.0 - 812752.0)
            msprime.PopulationParametersChange(
                time=656584.0, initial_size=37758.432874002, population_id=11
            ),
            # YNW, 812752.0 generations at 47873.117090463435,
            # [812752.0 - 1005916.0)
            msprime.PopulationParametersChange(
                time=812752.0, initial_size=47873.117090463435, population_id=11
            ),
            # YNW, 1005916.0 generations at 47873.117090463435,
            # [1005916.0 - 1244840.0)
            msprime.PopulationParametersChange(
                time=1005916.0, initial_size=47873.117090463435, population_id=11
            ),
            # YNW, 1244840.0 generations at 66348.19532908706,
            # [1244840.0 - 1540364.0)
            msprime.PopulationParametersChange(
                time=1244840.0, initial_size=66348.19532908706, population_id=11
            ),
            # YNW, 1540364.0 generations at 66348.19532908706,
            # [1540364.0 - 1905896.0)
            msprime.PopulationParametersChange(
                time=1540364.0, initial_size=66348.19532908706, population_id=11
            ),
            # YNW, 1905896.0 generations at 66348.19532908706,
            # [1905896.0 - inf)
            msprime.PopulationParametersChange(
                time=1905896.0, initial_size=66348.19532908706, population_id=11
            ),
            # SCW, 584.0 generations at 51760.3701901676,
            # [584.0 - 1307.0)
            msprime.PopulationParametersChange(
                time=584.0, initial_size=51760.3701901676, population_id=12
            ),
            # SCW, 1307.0 generations at 11188.24786444319,
            # [1307.0 - 2201.0)
            msprime.PopulationParametersChange(
                time=1307.0, initial_size=11188.24786444319, population_id=12
            ),
            # SCW, 2201.0 generations at 9883.083126612177,
            # [2201.0 - 3307.0)
            msprime.PopulationParametersChange(
                time=2201.0, initial_size=9883.083126612177, population_id=12
            ),
            # SCW, 3307.0 generations at 11182.117557601885,
            # [3307.0 - 4674.0)
            msprime.PopulationParametersChange(
                time=3307.0, initial_size=11182.117557601885, population_id=12
            ),
            # SCW, 4674.0 generations at 13164.21702528188,
            # [4674.0 - 6366.0)
            msprime.PopulationParametersChange(
                time=4674.0, initial_size=13164.21702528188, population_id=12
            ),
            # SCW, 6366.0 generations at 16184.241404144785,
            # [6366.0 - 8458.0)
            msprime.PopulationParametersChange(
                time=6366.0, initial_size=16184.241404144785, population_id=12
            ),
            # SCW, 8458.0 generations at 20373.79807324992,
            # [8458.0 - 11046.0)
            msprime.PopulationParametersChange(
                time=8458.0, initial_size=20373.79807324992, population_id=12
            ),
            # SCW, 11046.0 generations at 25202.02573882889,
            # [11046.0 - 14248.0)
            msprime.PopulationParametersChange(
                time=11046.0, initial_size=25202.02573882889, population_id=12
            ),
            # SCW, 14248.0 generations at 30145.04287378723,
            # [14248.0 - 18207.0)
            msprime.PopulationParametersChange(
                time=14248.0, initial_size=30145.04287378723, population_id=12
            ),
            # SCW, 18207.0 generations at 34073.930206368765,
            # [18207.0 - 23105.0)
            msprime.PopulationParametersChange(
                time=18207.0, initial_size=34073.930206368765, population_id=12
            ),
            # SCW, 23105.0 generations at 36220.40114094264,
            # [23105.0 - 29163.0)
            msprime.PopulationParametersChange(
                time=23105.0, initial_size=36220.40114094264, population_id=12
            ),
            # SCW, 29163.0 generations at 38334.812418179135,
            # [29163.0 - 36656.0)
            msprime.PopulationParametersChange(
                time=29163.0, initial_size=38334.812418179135, population_id=12
            ),
            # SCW, 36656.0 generations at 42224.55864780073,
            # [36656.0 - 45924.0)
            msprime.PopulationParametersChange(
                time=36656.0, initial_size=42224.55864780073, population_id=12
            ),
            # SCW, 45924.0 generations at 48756.10975000305,
            # [45924.0 - 57388.0)
            msprime.PopulationParametersChange(
                time=45924.0, initial_size=48756.10975000305, population_id=12
            ),
            # SCW, 57388.0 generations at 56966.42398970047,
            # [57388.0 - 71567.0)
            msprime.PopulationParametersChange(
                time=57388.0, initial_size=56966.42398970047, population_id=12
            ),
            # SCW, 71567.0 generations at 64192.886786215226,
            # [71567.0 - 89106.0)
            msprime.PopulationParametersChange(
                time=71567.0, initial_size=64192.886786215226, population_id=12
            ),
            # SCW, 89106.0 generations at 66432.60246398523,
            # [89106.0 - 110799.0)
            msprime.PopulationParametersChange(
                time=89106.0, initial_size=66432.60246398523, population_id=12
            ),
            # SCW, 110799.0 generations at 61044.09825658056,
            # [110799.0 - 137632.0)
            msprime.PopulationParametersChange(
                time=110799.0, initial_size=61044.09825658056, population_id=12
            ),
            # SCW, 137632.0 generations at 49929.848562769315,
            # [137632.0 - 170822.0)
            msprime.PopulationParametersChange(
                time=137632.0, initial_size=49929.848562769315, population_id=12
            ),
            # SCW, 170822.0 generations at 37914.763819457476,
            # [170822.0 - 211874.0)
            msprime.PopulationParametersChange(
                time=170822.0, initial_size=37914.763819457476, population_id=12
            ),
            # SCW, 211874.0 generations at 28990.00279753527,
            # [211874.0 - 262651.0)
            msprime.PopulationParametersChange(
                time=211874.0, initial_size=28990.00279753527, population_id=12
            ),
            # SCW, 262651.0 generations at 23837.504499328974,
            # [262651.0 - 325458.0)
            msprime.PopulationParametersChange(
                time=262651.0, initial_size=23837.504499328974, population_id=12
            ),
            # SCW, 325458.0 generations at 21373.93811603697,
            # [325458.0 - 403144.0)
            msprime.PopulationParametersChange(
                time=325458.0, initial_size=21373.93811603697, population_id=12
            ),
            # SCW, 403144.0 generations at 20910.24382389811,
            # [403144.0 - 499232.0)
            msprime.PopulationParametersChange(
                time=403144.0, initial_size=20910.24382389811, population_id=12
            ),
            # SCW, 499232.0 generations at 22117.58148669959,
            # [499232.0 - 618088.0)
            msprime.PopulationParametersChange(
                time=499232.0, initial_size=22117.58148669959, population_id=12
            ),
            # SCW, 618088.0 generations at 24984.97778210851,
            # [618088.0 - 765096.0)
            msprime.PopulationParametersChange(
                time=618088.0, initial_size=24984.97778210851, population_id=12
            ),
            # SCW, 765096.0 generations at 31521.677457587582,
            # [765096.0 - 946936.0)
            msprime.PopulationParametersChange(
                time=765096.0, initial_size=31521.677457587582, population_id=12
            ),
            # SCW, 946936.0 generations at 31521.677457587582,
            # [946936.0 - 1171848.0)
            msprime.PopulationParametersChange(
                time=946936.0, initial_size=31521.677457587582, population_id=12
            ),
            # SCW, 1171848.0 generations at 42066.29648325761,
            # [1171848.0 - 1450044.0)
            msprime.PopulationParametersChange(
                time=1171848.0, initial_size=42066.29648325761, population_id=12
            ),
            # SCW, 1450044.0 generations at 42066.29648325761,
            # [1450044.0 - 1794148.0)
            msprime.PopulationParametersChange(
                time=1450044.0, initial_size=42066.29648325761, population_id=12
            ),
            # SCW, 1794148.0 generations at 42066.29648325761,
            # [1794148.0 - inf)
            msprime.PopulationParametersChange(
                time=1794148.0, initial_size=42066.29648325761, population_id=12
            ),
        ],
    )


_species.add_demographic_model(_WildBoar_13W25())
