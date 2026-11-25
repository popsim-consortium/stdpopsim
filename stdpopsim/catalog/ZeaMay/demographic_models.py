import msprime

import stdpopsim

_species = stdpopsim.get_species("ZeaMay")


def _maize_arg_size_history():
    id = "Domestication_1M26"
    description = "Maize size history from ARG-based pairwise coalescence rates"
    long_description = """
        Piecewise-constant size history for a single maize population,
        inferred from pairwise coalescence rates estimated from an ancestral
        recombination graph (ARG) (Magalang et al., unpublished; presented in
        the cited talk). The inferred history was truncated to the
        interval between 10 and 1,000,000 generations ago: the size at 10
        generations is held constant to the present, and the size at 1,000,000
        generations is held constant into the past. The inference assumed a
        mutation rate of 7.3e-8 per bp per generation and a generation time
        of 1 year.
    """
    citations = [
        stdpopsim.Citation(
            author="Magalang et al.",
            year=2026,
            doi="https://www.youtube.com/watch?v=9K2iWaLZQLY",
            reasons={stdpopsim.CiteReason.DEM_MODEL},
        )
    ]

    # Size-change times (generations ago) and the population size from each
    # time back to the next one. N_present applies from 0 to the first
    # change time; the last size applies indefinitely into the past.
    N_present = 37541.46690608026
    size_changes = [
        (10.693114888921986, 22962.34489733273),
        (14.594442040531693, 16028.711413520263),
        (19.919148039370977, 10391.84776497947),
        (27.186545228139586, 8656.378677071596),
        (37.105414347079645, 8290.749694997996),
        (50.64313108983731, 8725.268369126936),
        (69.12001312240564, 9372.012936071196),
        (94.33808911945931, 9628.093589091939),
        (128.7568485114504, 9191.527316930018),
        (175.73311260955938, 7989.08678650474),
        (239.8484214585115, 6469.407867680865),
        (327.3558660737588, 6855.7484087491),
        (446.7899450880375, 9529.72994699974),
        (609.798924412106, 14163.610760522828),
        (832.2808789729802, 19153.164114963092),
        (1135.9342133504822, 19325.5252394869),
        (1550.3738817746778, 19022.99414904503),
        (2116.019699943183, 16693.623534348033),
        (2888.038442328698, 12052.597268108708),
        (3941.7241931123463, 8783.121201955355),
        (5379.84169007084, 6431.112069648686),
        (7342.648849150294, 5300.7331003489935),
        (10021.575954815564, 11296.791282384609),
        (13677.895631596157, 29124.27241616258),
        (18668.20445729789, 42848.927629585625),
        (25479.201409786452, 49133.13683204264),
        (34775.155048544446, 57289.485241967704),
        (47462.68884965199, 64094.93705716451),
        (64779.20313206969, 66720.96606834889),
        (88413.5572621844, 70158.10460701899),
        (120670.7821307497, 74803.77653057255),
        (164696.88711727667, 81577.23767685854),
        (224785.68670194177, 93430.92396801861),
        (306797.5711652845, 113584.27693251942),
        (418731.0635918018, 138257.24238835784),
        (571502.9064629112, 161929.56581549303),
        (780012.7587714745, 178972.8550065785),
    ]

    model = msprime.Demography()
    model.add_population(
        initial_size=N_present,
        name="Maize",
        description="Maize",
    )
    for time, size in size_changes:
        model.add_population_parameters_change(
            time=time, initial_size=size, population="Maize"
        )

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        citations=citations,
        generation_time=1,
        mutation_rate=7.3e-8,
        model=model,
    )


_species.add_demographic_model(_maize_arg_size_history())
