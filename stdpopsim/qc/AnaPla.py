import msprime

import stdpopsim


_species = stdpopsim.get_species("AnaPla")


def LavertskyEtAl2019TwoPop():
    id = "QC-MallardBlackDuck_2L19"

    # Parameters are taken from Fig 6 of Lavertsky et al. (2019)
    # analysis of contemporary samples
    generation_time = 4  # 4 years per generation
    # (see page 7, section 2.6, in Lavertsky et al. (2019)
    N_Mallard = 1.37e6  # Mallard estimated Ne
    # N_Mallard       = 10.37e6  # Mallard estimated Ne
    N_Black_duck = 1.57e6  # Black duck estimated Ne
    N_anc = 819535  # ancestral population size; not reported in the paper,
    # but reported to Peter Ralph in personal communication
    T_div = 632305 / generation_time  # Divergence time

    # symmetric migration model. Reported rates correspond to
    # number of migrants per generation. scaled by the ancestral Ne
    # so m_ij = M_ij / 2N_anc
    # (m_ij is the simulated rate and M_ij is the rate reported in paper)
    m_Mallard_Black = 2.82 / (2 * N_anc)
    m_Black_Mallard = m_Mallard_Black

    model = msprime.Demography()
    model.add_population(name="Mallard", description="Mallard", initial_size=N_Mallard)
    model.add_population(
        name="Black_duck", description="Black_duck", initial_size=N_Black_duck
    )
    model.add_population(name="Ancestral", description="Ancestral", initial_size=N_anc)
    model.add_population_split(
        time=T_div, derived=["Mallard", "Black_duck"], ancestral="Ancestral"
    )
    model.set_migration_rate(source="Mallard", dest="Black_duck", rate=m_Mallard_Black)
    model.set_migration_rate(source="Black_duck", dest="Mallard", rate=m_Black_Mallard)

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=generation_time,
        mutation_rate=4.83e-9,
        model=model,
    )


_species.get_demographic_model("MallardBlackDuck_2L19").register_qc(
    LavertskyEtAl2019TwoPop
)
