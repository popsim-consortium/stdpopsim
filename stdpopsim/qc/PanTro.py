import msprime
import stdpopsim

_species = stdpopsim.get_species("PanTro")


def KuhlwilmPan():
    id = "PanTro"
    populations = [
        stdpopsim.Population(id="western", description=""),
        stdpopsim.Population(id="central", description=""),
        stdpopsim.Population(id="bonobo", description=""),
        stdpopsim.Population(id="ghost", description=""),
    ]
    citations = [
        stdpopsim.Citation(
            author="Kuhlwilm et al.",
            year="2019",
            doi="",
            reasons={stdpopsim.CiteReason.DEM_MODEL},
        )
    ]

    generation_time = 25
    mutation_rate = 1.2e-8  # per bp per generation

    # parameter value definitions based on published values
    # these are parameters from Table S3, B7 in Kuhlwilm et al., 2019
    # corresponding to Table S4, Column H: "ABC flow AT"
    Ne_west = 10100
    Ne_cent = 73000
    Ne_bonobo = 33500
    Ne_ghost = 10000

    Ne_WAnc = 11000
    Ne_BAnc = 3900
    Ne_BonAnc = 11100
    Ne_CAnc = 103100
    Ne_ChimpAnc = 10100
    Ne_PanAnc = 14200
    Ne_Anc = 10000

    m_W_C = 57.5 / 4000
    m_C_W = 89.7 / 4000
    m_C_B = 14.5 / 4000
    m_B_C = 21.7 / 4000
    m_B_gh = 68 / 4000
    m_Ch_B = 0.0004 / 4000
    m_B_Ch = 0.0004 / 4000

    T_chimp_ex = 1.00075 * 4000
    T_bon_ex = 1.5505 * 4000
    T_west_re = 2.61 * 4000
    T_bon_re = 3.08 * 4000
    T_cent_re = 3.78 * 4000
    T_ghost_ex = 5 * 4000
    T_ChimpAnc = 7 * 4000
    T_PanAnc = 19.9 * 4000
    T_BonAnc_re = 19.875 * 4000
    T_Anc = 33.01 * 4000
    T_Anc_ex = 12 * 4000
    T_Anc_ex_stop = 15 * 4000

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        mutation_rate=mutation_rate,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=Ne_west, metadata=populations[0].asdict()
            ),
            msprime.PopulationConfiguration(
                initial_size=Ne_cent, metadata=populations[1].asdict()
            ),
            msprime.PopulationConfiguration(
                initial_size=Ne_bonobo, metadata=populations[2].asdict()
            ),
            msprime.PopulationConfiguration(
                initial_size=Ne_ghost, metadata=populations[3].asdict()
            ),
        ],
        demographic_events=[
            msprime.MassMigration(
                time=T_chimp_ex,
                source=0,
                destination=1,
                proportion=m_W_C,
            ),
            msprime.MassMigration(
                time=T_chimp_ex,
                source=1,
                destination=0,
                proportion=m_C_W,
            ),
            msprime.MassMigration(
                time=T_bon_ex,
                source=1,
                destination=2,
                proportion=m_C_B,
            ),
            msprime.MassMigration(
                time=T_bon_ex,
                source=2,
                destination=1,
                proportion=m_B_C,
            ),
            msprime.PopulationParametersChange(
                time=T_west_re,
                initial_size=Ne_WAnc,
                growth_rate=0.0,
                population_id=0,
            ),
            msprime.PopulationParametersChange(
                time=T_bon_re,
                initial_size=Ne_BAnc,
                growth_rate=0.0,
                population_id=2,
            ),
            msprime.PopulationParametersChange(
                time=T_cent_re,
                initial_size=Ne_CAnc,
                growth_rate=0.0,
                population_id=1,
            ),
            msprime.MassMigration(
                time=T_ghost_ex,
                source=2,
                destination=3,
                proportion=m_B_gh,
            ),
            msprime.MassMigration(
                time=T_ChimpAnc,
                source=0,
                destination=1,
                proportion=1.0,
            ),
            msprime.PopulationParametersChange(
                time=T_ChimpAnc,
                initial_size=Ne_ChimpAnc,
                growth_rate=0.0,
                population_id=1,
            ),
            msprime.MigrationRateChange(
                time=T_Anc_ex,
                rate=m_Ch_B,
                matrix_index=(2, 1),
            ),
            msprime.MigrationRateChange(
                time=T_Anc_ex,
                rate=m_B_Ch,
                matrix_index=(1, 2),
            ),
            msprime.MigrationRateChange(
                time=T_Anc_ex_stop,
                rate=0,
                matrix_index=(2, 1),
            ),
            msprime.MigrationRateChange(
                time=T_Anc_ex_stop,
                rate=0,
                matrix_index=(1, 2),
            ),
            msprime.PopulationParametersChange(
                time=T_BonAnc_re,
                initial_size=Ne_BonAnc,
                growth_rate=0.0,
                population_id=2,
            ),
            msprime.MassMigration(
                time=T_PanAnc,
                source=1,
                destination=2,
                proportion=1.0,
            ),
            msprime.PopulationParametersChange(
                time=T_PanAnc,
                initial_size=Ne_PanAnc,
                growth_rate=0.0,
                population_id=2,
            ),
            msprime.MassMigration(
                time=T_Anc,
                source=2,
                destination=3,
                proportion=1.0,
            ),
            msprime.PopulationParametersChange(
                time=T_Anc,
                initial_size=Ne_Anc,
                growth_rate=0.0,
                population_id=3,
            ),
        ],
    )


_species.get_demographic_model("BonoboGhost_4K19").register_qc(KuhlwilmPan)
