import msprime
import stdpopsim

_species = stdpopsim.get_species("OrySat")

# population definitions that are reused.
_rufipogon = stdpopsim.Population(id="RUF", description="Oryza rufipogon (wild rice)")
_indica = stdpopsim.Population(id="IND", description="Oryza sativa indica")
_tropical_japonica = stdpopsim.Population(
    id="TRJ", description="Tropical Oryza sativa japonica"
)


def _BottMig_3C07():
    id = "BottleneckMigration_3C07"
    description = (
        "Bottleneck + migration model for the origin of domesticated rice varieties"
    )
    long_description = """
    The bottleneck + migration model of domesticated rice varieties
    (indica and tropical japonica) from Caicedo et al. (2007).
    Parameters were inferred from the site frequency spectrum (SFS),
    with parameter values taken from Table 2.
    """
    populations = [_rufipogon, _indica, _tropical_japonica]

    citations = [
        stdpopsim.Citation(
            author="Caicedo et al.",
            year=2007,
            doi="https://doi.org/10.1371/journal.pgen.0030163",
            reasons={stdpopsim.CiteReason.DEM_MODEL},
        ),
    ]

    generation_time = 1
    mutation_rate = 6.5e-9

    # First we set out the maximum likelihood values of the various parameters
    # given in Table 2.
    # To rescale the parameters, Caicedo et al. assumed that rice cultivation
    # began around 12,000 years (and therefore, generations) ago
    # T_1 × 2 x N_rufi = 12000

    T_1 = 12000  # time to the start of domestication
    T_B = 9000  # bottleneck lasted for 3000 years/generations

    N_rufi = 12000 / (0.04 * 2)  # ancestral population size = 150000
    N_trj = 0.12 * 150000  # population size for tropical japonica
    N_ind = 0.27 * 150000  # population size for indica
    N_B = (
        0.0055 * 150000
    )  # bottleneck population size for both indica and tropical japonica

    # to japonica = 0.42, to indica = 0.945, to rufipogon = 3.5

    # Migration rates in the form M_from_to, thinking forward in time
    m_trj_rufi = 3.5 / N_rufi
    m_trj_ind = 0.945 / N_ind
    m_rufi_trj = 0.42 / N_trj
    m_rufi_ind = 0.945 / N_ind
    m_ind_rufi = 3.5 / N_rufi
    m_ind_trj = 0.42 / N_trj

    # This actually just translates to symmetrical
    # migration between all the populations

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        mutation_rate=mutation_rate,
        # Population IDs correspond to their indexes in the population
        # configuration array. Therefore, we have 0=RUF, 1=IND and 2=TRJ
        # initially.
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_rufi, metadata=populations[0].asdict()
            ),
            msprime.PopulationConfiguration(
                initial_size=N_ind, metadata=populations[1].asdict()
            ),
            msprime.PopulationConfiguration(
                initial_size=N_trj, metadata=populations[2].asdict()
            ),
        ],
        migration_matrix=[
            [0, m_ind_rufi, m_trj_rufi],
            [m_rufi_ind, 0, m_trj_ind],
            [m_rufi_trj, m_ind_trj, 0],
        ],
        demographic_events=[
            # There is migration between all populations
            msprime.MigrationRateChange(time=0, rate=m_ind_rufi, source=0, dest=1),
            msprime.MigrationRateChange(time=0, rate=m_ind_trj, source=2, dest=1),
            msprime.MigrationRateChange(time=0, rate=m_trj_rufi, source=0, dest=2),
            msprime.MigrationRateChange(time=0, rate=m_trj_ind, source=2, dest=1),
            msprime.MigrationRateChange(time=0, rate=m_rufi_trj, source=2, dest=0),
            msprime.MigrationRateChange(time=0, rate=m_rufi_ind, source=1, dest=0),
            # IND and TRJ undergo independent bottlenecks at T_B
            msprime.PopulationParametersChange(
                time=T_B, initial_size=N_B, population_id=1
            ),
            msprime.PopulationParametersChange(
                time=T_B, initial_size=N_B, population_id=2
            ),
            # Migration stops before that point
            msprime.MigrationRateChange(time=T_B, rate=0),
            # IND and TRJ both merge onto RUF at T_1
            msprime.MassMigration(time=T_1, source=1, destination=0, proportion=1.0),
            msprime.MassMigration(time=T_1, source=2, destination=0, proportion=1.0),
        ],
    )


_species.add_demographic_model(_BottMig_3C07())
