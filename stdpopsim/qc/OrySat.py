import msprime
import stdpopsim

_species = stdpopsim.get_species("OrySat")

def _BottMig_3C07_qc():
    id = "QC-BottMig_3C07"
    
    #populations
    rufipogon = stdpopsim.Population(id = "RUF", description="O. rufipogon population")
    indica = stdpopsim.Population(id = "IND", description="O. sativa indica population")
    japonica = stdpopsim.Population(id = "TRJ", description="O. sativa tropical japonica population")
    populations = [rufipogon, indica, japonica]

    # generation time of 1 year and a mutation rate of 3.2e-9 per bp per gen
    generation_time = 1
    mutation_rate = 3.2e-9

    # All demographic parameters were extracted from Caicedo et al. (2007)
    # Table 2, bottleneck + migration model

    # Times
    T_1 = 12000  # time when domestication started
    T_B = 9000  # (12000-3000), bottleneck lasted for 3000 years 

    # Population sizes
    N_anc = 150000 # Ancestral population size (12000 / (0.04 * 2))
    N_rufi = N_anc  # population size for rufipogon
    N_ind = 0.27 * N_anc  # population size for indica
    N_trj = 0.12 * N_anc  # population size for tropical japonica
    N_B = (0.0055 * 150000) # bottleneck population size for both indica and tropical japonica

    # Migration rates 
    # to japonica = 0.42, to indica = 0.945, to rufipogon = 3.5
    m_ind = 0.945 / N_ind
    m_trj = 0.42 / N_trj
    m_rufi = 3.5 / N_rufi
    # they are symmetrical
    
    migration_matrix = [
        [0, m_rufi, m_rufi],
        [m_ind, 0, m_ind],
        [m_trj, m_trj, 0],
        ]

    demographic_events = [
            # There is migration between all populations:
            # (forward in time, indexes are flipped because this is backwards)
        
            # to rufipogon
            # from indica to rufipogon 
            msprime.MigrationRateChange(time = 0, rate = m_rufi, source = 0, dest = 1),
            # from japonica to rufipogon
            msprime.MigrationRateChange(time = 0, rate = m_rufi, source = 0, dest = 2),
            
            # to indica
            # from rufipogon to indica
            msprime.MigrationRateChange(time = 0, rate = m_ind, source = 1, dest = 0),
            # from japonica to indica
            msprime.MigrationRateChange(time = 0, rate = m_ind, source = 1, dest = 2),
            
            # to japonica
            # from rufipogon to japonica
            msprime.MigrationRateChange(time = 0, rate = m_trj, source = 2, dest = 0),
            # from indica to japonica
            msprime.MigrationRateChange(time = 0, rate = m_trj, source = 2, dest = 1),
            
            # Independent bottlenecks at T_B for indica and japonica
            
            # bottleneck for indica
            msprime.PopulationParametersChange(time = T_B, initial_size = N_B, population_id = 1), 
            # bottleneck for japonica
            msprime.PopulationParametersChange(time = T_B, initial_size = N_B, population_id = 2),
            
            # Migration stops before that
            msprime.MigrationRateChange(time = T_B, rate = 0),
            
            # Both merge onto rufipogon at T_1
            # indica merges onto rufipogon
            msprime.MassMigration(time = T_1, source = 1, destination = 0, proportion = 1.0),
            # japonica merges onto rufipogon
            msprime.MassMigration(time = T_1, source = 2, destination = 0, proportion = 1.0),
        ]

    return stdpopsim.DemographicModel(
        id = id,
        description = id,
        long_description = id,
        generation_time = generation_time,
        populations = populations,
        population_configurations = [
            msprime.PopulationConfiguration(initial_size = N_rufi, metadata = populations[0].asdict()
            ),
            msprime.PopulationConfiguration(initial_size = N_ind, metadata = populations[1].asdict()
            ),
            msprime.PopulationConfiguration(initial_size = N_trj, metadata = populations[2].asdict()
            ),
        ],
        migration_matrix = migration_matrix,
        demographic_events = demographic_events,
        mutation_rate = mutation_rate,
    )


_species.get_demographic_model("BottleneckMigration_3C07").register_qc(_BottMig_3C07_qc())
