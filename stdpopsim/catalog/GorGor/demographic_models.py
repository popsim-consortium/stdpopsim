import msprime
import stdpopsim

_species = stdpopsim.get_species("GorGor")

###########################################################
#
# Demographic models
#
###########################################################

_pawar2023 = stdpopsim.Citation(
    author="Pawar et al. 2023",
    year=2023,
    doi="https://doi.org/10.1038/s41559-023-02145-2",
    reasons={stdpopsim.CiteReason.DEM_MODEL},
)


def _gorilla_ghost():
    id = "GorillaGhost_5P23"
    description = "Ghost admixture in eastern gorillas"
    long_description = """
        Demographic model of ghost admixture into eastern gorillas
        from Pawar et al. (2023) Fig. 2A and Supp. Table 1.
        This model simulates five populations:
        mountain gorillas, eastern lowland gorillas, western lowland gorillas,
        cross river gorillas, and a extinct ghost lineage. The ghost admixture
        event is modelled as a 2.47% pulse from the ghost lineage to a
        population ancestral to eastern gorillas. Migration events among
        eastern and western lowland gorillas are modelled as single generation
        pulses. Population size changes are also modelled.
    """

    populations = [
        stdpopsim.Population(
            id="cross_river", description="Contemporary Cross River Gorillas"
        ),
        stdpopsim.Population(
            id="western_lowland", description="Contemporary Western Lowland Gorillas"
        ),
        stdpopsim.Population(
            id="eastern_lowland", description="Contemporary Eastern Lowland Gorillas"
        ),
        stdpopsim.Population(
            id="mountain", description="Contemporary Mountain Gorillas"
        ),
        stdpopsim.Population(
            id="ghost", description="Extinct ghost lineage", sampling_time=None
        ),
    ]

    citations = [_pawar2023]
    generation_time = 19
    mutation_rate = 1.235e-8

    N_cross_river = 14558
    N_western_lowland = 64749
    N_eastern_lowland = 20894
    N_mountain = 2158

    # The effective population size of the extinct ghost population is
    # difficult to estmiate due to the low admixture proportion into Eastern
    # Gorillas. Here, it is somewhat arbitrarly set to 25,000 individuals.
    # See also:
    # https://github.com/h-pawar/gor_ghost_introg/blob/main/4.Demog_ABC_model_comparison/4.1.model_comparison_simulations/scripts/mc.ghoste.12apr22.R
    N_ghost = 25000

    # *********************

    N_bottleneck_ELG = 243
    N_anc_ELG = 21982
    N_bottleneck_MG = 115
    N_anc_MG = 3009
    N_anc_EG = 5325
    N_change_WLG = 48288
    N_anc_WG = 98135
    N_gor_species_split = 14364
    N_anc_ghost_gor = 35381

    m_ghost_EG = 0.02477
    m_WLG_EG = 0.002768
    m_EG_WLG = 0.00827

    T_EG_subspecies_split = 15.048 * 1000 / generation_time
    T_ghost_introg_end = 38.281 * 1000 / generation_time
    T_ne_change_WLG = 40.896 * 1000 / generation_time
    T_WG_subspecies_split = 454.313 * 1000 / generation_time
    T_gor_species_split = 965.481 * 1000 / generation_time
    T_ghost_gor_split = 3421.360 * 1000 / generation_time

    T_ne_recovery_ELG = 30 / generation_time
    T_ne_bottleneck_ELG = 7.220 * 1000 / generation_time
    T_ne_recovery_MG = 9.120 * 1000 / generation_time
    T_ne_bottleneck_MG = 9.310 * 1000 / generation_time
    T_extant_admix_end = 34.0 * 1000 / generation_time
    T_extant_admix_start = 34.0019 * 1000 / generation_time
    T_ghost_introg_start = 38.300 * 1000 / generation_time

    population_configurations = [
        msprime.PopulationConfiguration(
            initial_size=N_cross_river, metadata=populations[0].asdict()
        ),
        msprime.PopulationConfiguration(
            initial_size=N_western_lowland, metadata=populations[1].asdict()
        ),
        msprime.PopulationConfiguration(
            initial_size=N_eastern_lowland, metadata=populations[2].asdict()
        ),
        msprime.PopulationConfiguration(
            initial_size=N_mountain, metadata=populations[3].asdict()
        ),
        msprime.PopulationConfiguration(
            initial_size=N_ghost, metadata=populations[4].asdict()
        ),
    ]

    migration_matrix = [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ]

    demographic_events = [
        # ELG bottleneck
        msprime.PopulationParametersChange(
            time=T_ne_recovery_ELG,
            initial_size=N_bottleneck_ELG,
            growth_rate=0.0,
            population_id=2,
        ),
        msprime.PopulationParametersChange(
            time=T_ne_bottleneck_ELG,
            initial_size=N_anc_ELG,
            growth_rate=0.0,
            population_id=2,
        ),
        # MG bottleneck
        msprime.PopulationParametersChange(
            time=T_ne_recovery_MG,
            initial_size=N_bottleneck_MG,
            growth_rate=0.0,
            population_id=3,
        ),
        msprime.PopulationParametersChange(
            time=T_ne_bottleneck_MG,
            initial_size=N_anc_MG,
            growth_rate=0.0,
            population_id=3,
        ),
        # EG subspecies split
        msprime.MassMigration(
            time=T_EG_subspecies_split,
            source=2,
            destination=3,
            proportion=1.0,
        ),
        msprime.PopulationParametersChange(
            time=T_EG_subspecies_split,
            initial_size=N_anc_EG,
            growth_rate=0.0,
            population_id=3,
        ),
        # Admixture between EG and WLG
        msprime.MigrationRateChange(
            time=T_extant_admix_end,
            rate=m_EG_WLG,
            matrix_index=(1, 3),
        ),
        msprime.MigrationRateChange(
            time=T_extant_admix_end,
            rate=m_WLG_EG,
            matrix_index=(3, 1),
        ),
        msprime.MigrationRateChange(
            time=T_extant_admix_start,
            rate=0,
            matrix_index=(1, 3),
        ),
        msprime.MigrationRateChange(
            time=T_extant_admix_start,
            rate=0,
            matrix_index=(3, 1),
        ),
        # Ghost admixture into EG
        msprime.MigrationRateChange(
            time=T_ghost_introg_end,
            rate=m_ghost_EG,
            matrix_index=(3, 4),
        ),
        msprime.MigrationRateChange(
            time=T_ghost_introg_start,
            rate=0,
            matrix_index=(3, 4),
        ),
        # Ancestral size change WLG
        msprime.PopulationParametersChange(
            time=T_ne_change_WLG,
            initial_size=N_change_WLG,
            growth_rate=0.0,
            population_id=1,
        ),
        # WG subspecies split
        msprime.MassMigration(
            time=T_WG_subspecies_split,
            source=0,
            destination=1,
            proportion=1.0,
        ),
        msprime.PopulationParametersChange(
            time=T_WG_subspecies_split,
            initial_size=N_anc_WG,
            growth_rate=0.0,
            population_id=1,
        ),
        # Gorilla species split
        msprime.MassMigration(
            time=T_gor_species_split,
            source=1,
            destination=3,
            proportion=1.0,
        ),
        msprime.PopulationParametersChange(
            time=T_gor_species_split,
            initial_size=N_gor_species_split,
            growth_rate=0.0,
            population_id=3,
        ),
        # Ghost Gorilla split
        msprime.MassMigration(
            time=T_ghost_gor_split,
            source=3,
            destination=4,
            proportion=1.0,
        ),
        msprime.PopulationParametersChange(
            time=T_ghost_gor_split,
            initial_size=N_anc_ghost_gor,
            growth_rate=0.0,
            population_id=4,
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


_species.add_demographic_model(_gorilla_ghost())
