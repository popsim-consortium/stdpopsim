import msprime
import stdpopsim

_species = stdpopsim.get_species("PanTro")

###########################################################
#
# Demographic models
#
###########################################################

_kuhlwilm2019 = stdpopsim.Citation(
    author="Kuhlwilm et al. 2019",
    year=2019,
    doi="https://doi.org/10.1038/s41559-019-0881-7",
    reasons={stdpopsim.CiteReason.DEM_MODEL},
)


def _bonobo_ghost():
    id = "BonoboGhost_4K19"
    description = "Ghost admixture into bonobos"
    long_description = """
        Demographic model of ghost admixture into bonobos from
        Kuhlwilm et al. (2019) Supplementary Table S3 row 7.
        This model simulates four populations:
        western chimpanzees, central chimpanzees, bonobos,
        and a extinct ghost lineage. The ghost admixture event is
        modelled as a 1.7% pulse from the ghost lineage to bonobos.
        Migration events among western chimpanzees, central chimpanzees,
        and bonobos are modelled as single generation pulses.
        Populatio size changes are also modelled.
    """

    populations = [
        stdpopsim.Population(
            id="western", description="Contemporary Western Chimpanzees"
        ),
        stdpopsim.Population(
            id="central", description="Contemporary Central Chimpanzees"
        ),
        stdpopsim.Population(id="bonobo", description="Contemporary Bonobos"),
        stdpopsim.Population(
            id="ghost", description="Extinct ghost lineage", sampling_time=None
        ),
    ]

    citations = [_kuhlwilm2019]
    generation_time = 25
    mutation_rate = 1.2e-8

    N_western = 10100
    N_central = 73000
    N_bonobo = 33500
    N_ghost = 10000

    N_anc_western = 11000
    N_anc_central = 103100
    N_anc_bonobo1 = 3900
    N_anc_bonobo2 = 11100
    N_anc_chimp = 10100
    N_anc_chimp_bonobo = 14200
    N_anc = 10000

    m_western_central = 0.014375
    m_central_western = 0.022425
    m_bonobo_central = 0.005425
    m_central_bonobo = 0.003625
    m_bonobo_ghost = 0.017
    m_bonobo_anc_chimp = 1e-7

    T_ghost_split = 3301 * 1000 / generation_time
    T_bonobo_split = 1990 * 1000 / generation_time
    T_bonobo_anc_chimp_migration_start = 1500 * 1000 / generation_time
    T_bonobo_anc_chimp_migration_stop = 1200 * 1000 / generation_time
    T_western_split = 700 * 1000 / generation_time
    T_bonobo_ghost_admixture = 500 * 1000 / generation_time
    T_bonobo_central_admixture = 155.05 * 1000 / generation_time
    T_western_central_admixture = 100.075 * 1000 / generation_time
    T_bonobo_resize1 = 308 * 1000 / generation_time
    T_bonobo_resize2 = 1987.5 * 1000 / generation_time
    T_central_resize = 378 * 1000 / generation_time
    T_western_resize = 261 * 1000 / generation_time

    population_configurations = [
        msprime.PopulationConfiguration(
            initial_size=N_western, metadata=populations[0].asdict()
        ),
        msprime.PopulationConfiguration(
            initial_size=N_central, metadata=populations[1].asdict()
        ),
        msprime.PopulationConfiguration(
            initial_size=N_bonobo, metadata=populations[2].asdict()
        ),
        msprime.PopulationConfiguration(
            initial_size=N_ghost, metadata=populations[3].asdict()
        ),
    ]
    migration_matrix = [
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
    ]
    demographic_events = [
        msprime.MassMigration(
            time=T_western_central_admixture,
            source=0,
            destination=1,
            proportion=m_western_central,
        ),
        msprime.MassMigration(
            time=T_western_central_admixture,
            source=1,
            destination=0,
            proportion=m_central_western,
        ),
        msprime.MassMigration(
            time=T_bonobo_central_admixture,
            source=1,
            destination=2,
            proportion=m_central_bonobo,
        ),
        msprime.MassMigration(
            time=T_bonobo_central_admixture,
            source=2,
            destination=1,
            proportion=m_bonobo_central,
        ),
        msprime.PopulationParametersChange(
            time=T_western_resize,
            initial_size=N_anc_western,
            growth_rate=0.0,
            population_id=0,
        ),
        msprime.PopulationParametersChange(
            time=T_bonobo_resize1,
            initial_size=N_anc_bonobo1,
            growth_rate=0.0,
            population_id=2,
        ),
        msprime.PopulationParametersChange(
            time=T_central_resize,
            initial_size=N_anc_central,
            growth_rate=0.0,
            population_id=1,
        ),
        msprime.MassMigration(
            time=T_bonobo_ghost_admixture,
            source=2,
            destination=3,
            proportion=m_bonobo_ghost,
        ),
        msprime.MassMigration(
            time=T_western_split,
            source=0,
            destination=1,
            proportion=1.0,
        ),
        msprime.PopulationParametersChange(
            time=T_western_split,
            initial_size=N_anc_chimp,
            growth_rate=0.0,
            population_id=1,
        ),
        msprime.MigrationRateChange(
            time=T_bonobo_anc_chimp_migration_stop,
            rate=m_bonobo_anc_chimp,
            matrix_index=(1, 2),
        ),
        msprime.MigrationRateChange(
            time=T_bonobo_anc_chimp_migration_stop,
            rate=m_bonobo_anc_chimp,
            matrix_index=(2, 1),
        ),
        msprime.MigrationRateChange(
            time=T_bonobo_anc_chimp_migration_start,
            rate=0,
            matrix_index=(1, 2),
        ),
        msprime.MigrationRateChange(
            time=T_bonobo_anc_chimp_migration_start,
            rate=0,
            matrix_index=(2, 1),
        ),
        msprime.PopulationParametersChange(
            time=T_bonobo_resize2,
            initial_size=N_anc_bonobo2,
            growth_rate=0.0,
            population_id=2,
        ),
        msprime.MassMigration(
            time=T_bonobo_split,
            source=1,
            destination=2,
            proportion=1.0,
        ),
        msprime.PopulationParametersChange(
            time=T_bonobo_split,
            initial_size=N_anc_chimp_bonobo,
            growth_rate=0.0,
            population_id=2,
        ),
        msprime.MassMigration(
            time=T_ghost_split,
            source=2,
            destination=3,
            proportion=1.0,
        ),
        msprime.PopulationParametersChange(
            time=T_ghost_split,
            initial_size=N_anc,
            growth_rate=0.0,
            population_id=3,
        ),
    ]

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        mutation_rate=mutation_rate,
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events,
    )


_species.add_demographic_model(_bonobo_ghost())
