import msprime
import stdpopsim

_species = stdpopsim.get_species("CanFam")


###########################################################
#
# Demographic models
#
###########################################################


def _dog_wolf_admixture():
    id = "EarlyWolfAdmixture_6F14"
    description = "Dog domestication and admixture with wolf (Freedman et al., 2014)"
    long_description = """
    Demographic model of dog domestication with the Boxer reference
    genome, two basal dog breeds, three wolves, and a golden jackal
    (outgroup) from Freedman et al. (2014).
    The topology was based on the neighbor-joining tree constructed
    from genome-wide pairwise divergence (Figure 5A).
    The parameters of the model were inferred using the
    Generalized Phylogenetic Coalescent Sampler (G-PhoCS) on neutral
    autosomal loci from whole genome sequences.
    Model parameters were calibrated with
    a generation interval of 3 years (page 3) and
    a mutation rate of 1e-8 (page 3).
    Estimated (calibrated) effective population sizes, divergence times,
    and migration rates are given in Table S12.
    Migration is constant continuous geneflow from the start and
    end times of the two populations that define it (section S9.2.3 in S9).
    While Boxer reference genome was used, no Boxer-specific parameters
    were estimated in the model, hence Boxer was removed from this model
    implementation. Thence the ancestral Basenji population with its specific
    Ne is for the Ancestral (Boxer, Basenji) population.
    """
    citations = [
        stdpopsim.Citation(
            author="Freedman et al.",
            year="2014",
            doi="https://doi.org/10.1371/journal.pgen.1004016",
            reasons={stdpopsim.CiteReason.DEM_MODEL},
        )
    ]

    generation_time = 3  # page 3
    mutation_rate = 1e-8  # page 3

    # Estimated (calibrated) effective population sizes, from Table S12
    # no N_BOX provided in the paper
    N_BSJ = 2639
    N_DNG = 1914
    N_ISW = 26092
    N_CRW = 11427
    N_CHW = 5426
    N_GLJ = 19446
    N_ancDOG1 = 793  # (Boxer, Basenji)
    N_ancDOG = 1999  # ((Boxer, Basenji), Dingo)
    N_ancWLF1 = 1393  # (Israeli, Croatian) wolf
    N_ancWLF = 12627  # ((Israeli, Croatian), Chinese) wolf
    N_ancDW = 44993  # (dog, wolf)
    N_root = 18169  # ((dog, wolf), golden jackal)root

    # Estimated (calibrated) divergence times, from Table S12
    T_ancDOG1 = 12102  # (Boxer, Basenji)
    T_ancDOG = 12795  # ((Boxer, Basenji), Dingo)
    T_ancWLF1 = 13389  # (Israeli, Croatian) wolf
    T_ancWLF = 13455  # ((Israeli, Croatian), Chinese) wolf
    T_ancDW = 14874  # (dog, wolf)
    T_ancroot = 398262  # ((dog, wolf), golden jackal)root

    # Migration rates, from Table S12
    # migration is constant continuous geneflow from the start and end times of
    # the two populations that define it. See section S9.2.3 in S9.
    m_ISW_BSJ = 0.18
    m_BSJ_ISW = 0.07
    m_CHW_DNG = 0.03
    m_DNG_CHW = 0.04
    m_GLJ_ISW = 0
    m_ISW_GLJ = 0.05
    m_GLJ_ancDW = 0.02
    m_ancDW_GLJ = 0.99

    model = msprime.Demography()
    model.add_population(
        initial_size=N_BSJ,
        name="BSJ",
        description="Basenji",
    )
    model.add_population(
        initial_size=N_DNG,
        name="DNG",
        description="Dingo",
    )
    model.add_population(
        initial_size=N_CHW,
        name="CHW",
        description="Chinese wolf",
    )
    model.add_population(
        initial_size=N_ISW,
        name="ISW",
        description="Israeli wolf",
    )
    model.add_population(
        initial_size=N_CRW,
        name="CRW",
        description="Croatian wolf",
    )
    model.add_population(
        initial_size=N_GLJ,
        name="GLJ",
        description="Golden jackal",
    )
    # model.add_population(
    #     initial_size=N_ancDOG1,
    #     name="ancDOG1",
    #     description="Ancestral (Boxer, Basenji)",
    # )
    # commenting ancDOG1 pop due to the removal of Boxer
    model.add_population(
        initial_size=N_ancDOG,
        name="ancDOG",
        description="Ancestral ((Boxer, Basenji), Dingo)",
    )
    model.add_population(
        initial_size=N_ancWLF1,
        name="ancWLF1",
        description="Ancestral (Israeli, Croatian) wolf",
    )
    model.add_population(
        initial_size=N_ancWLF,
        name="ancWLF",
        description="Ancestral ((Israeli, Croatian), Chinese) wolf",
    )
    model.add_population(
        initial_size=N_ancDW,
        name="ancDW",
        description="Ancestral (dog, wolf)",
    )
    model.add_population(
        initial_size=N_root,
        name="root",
        description="Ancestral ((dog, wolf), golden jackal) root",
    )

    model.set_migration_rate(rate=m_ISW_BSJ, source="BSJ", dest="ISW")
    model.set_migration_rate(rate=m_BSJ_ISW, source="ISW", dest="BSJ")
    model.set_migration_rate(rate=m_CHW_DNG, source="DNG", dest="CHW")
    model.set_migration_rate(rate=m_DNG_CHW, source="CHW", dest="DNG")
    model.set_migration_rate(rate=m_GLJ_ISW, source="ISW", dest="GLJ")
    model.set_migration_rate(rate=m_ISW_GLJ, source="GLJ", dest="ISW")

    # model.add_population_split(
    #     time=T_ancDOG1, derived=["BOX", "BSJ"], ancestral="ancDOG1"
    # )
    # commenting out pop split above due to the removal of Boxer and
    # instead changing Ne below for "Ancestral (Boxer, Basenji)"
    model.add_population_parameters_change(
        time=T_ancDOG1,
        initial_size=N_ancDOG1,
        population="BSJ",
    )
    # In the paper, the below split is (ancDOG1, DNG)ancDOG
    model.add_population_split(
        time=T_ancDOG, derived=["BSJ", "DNG"], ancestral="ancDOG"
    )
    model.add_population_split(
        time=T_ancWLF1, derived=["ISW", "CRW"], ancestral="ancWLF1"
    )
    model.add_population_split(
        time=T_ancWLF, derived=["ancWLF1", "CHW"], ancestral="ancWLF"
    )
    model.add_population_split(
        time=T_ancDW, derived=["ancDOG", "ancWLF"], ancestral="ancDW"
    )
    model.add_migration_rate_change(
        time=T_ancDW, rate=m_GLJ_ancDW, source="ancDW", dest="GLJ"
    )
    model.add_migration_rate_change(
        time=T_ancDW, rate=m_ancDW_GLJ, source="GLJ", dest="ancDW"
    )
    model.add_population_split(
        time=T_ancroot, derived=["ancDW", "GLJ"], ancestral="root"
    )

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        citations=citations,
        generation_time=generation_time,
        mutation_rate=mutation_rate,
        model=model,
    )


_species.add_demographic_model(_dog_wolf_admixture())
