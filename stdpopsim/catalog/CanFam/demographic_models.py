import msprime
import stdpopsim

_species = stdpopsim.get_species("CanFam")


###########################################################
#
# Demographic models
#
###########################################################


def _dog_wolf_admixture():
    id = "EarlyWolfAdmixture_7F14"
    description = "Dog domestication with wolf and jackal"
    long_description = """
    Demographic model of dog domestication with the Boxer reference
    genome, two basal dog breeds, three wolves, and a golden jackal.
    The topology was based on the neighbor-joining tree constructed
    from genome-wide pairwise divergence. The parameters of the model
    were inferred using the Generalized Phylogenetic Coalescent Sampler
    (G-PhoCS) (Gronau et al. 2011) on neutral autosomal loci from
    whole genome sequences. See Figure 5A in Freedman et al. 2014.
    Note that Freedman et al. 2014 calibrated the parameters with a
    mutation rate of 1e-8.
    """
    citations = [
        stdpopsim.Citation(
            author="Freedman et al.",
            year="2014",
            doi="https://doi.org/10.1371/journal.pgen.1004016",
            reasons={stdpopsim.CiteReason.DEM_MODEL},
        )
    ]

    generation_time = 3

    #  parameter value definitions based on published values

    # Estimated (calibrated) effective population sizes, from Table S12
    # ******Not correct - putting in for flake8
    N_BOX = 1000  # Have not found this value yet
    # ******Not correct

    N_BSJ = 2639
    N_DNG = 1914
    N_ISW = 26092
    N_CRW = 11427
    N_CHW = 5426
    N_GLJ = 19446
    N_ancDOG1 = 793
    N_ancDOG = 1999
    N_ancWLF1 = 1393
    N_ancWLF = 12627
    N_ancDW = 44993
    N_root = 18169

    # Estimated (calibrated) divergence times, from Table S12
    T_ancDOG1 = 12102
    T_ancDOG = 12795
    T_ancWLF1 = 13389
    T_ancWLF = 13455
    T_ancDW = 14874
    T_ancroot = 398262

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
        initial_size=N_BOX,
        name="BOX",
        description="Boxer",
    )
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
    model.add_population(
        initial_size=N_ancDOG1,
        name="ancDOG1",
        description="Ancestral Boxer and Basenji",
    )
    model.add_population(
        initial_size=N_ancDOG,
        name="ancDOG",
        description="Ancestral Boxer-Basenji and Dingo",
    )
    model.add_population(
        initial_size=N_ancWLF1,
        name="ancWLF1",
        description="Ancestral Israeli and Croatian wolf",
    )
    model.add_population(
        initial_size=N_ancWLF,
        name="ancWLF",
        description="Ancestral Israeli-Croatian and Chinese wolf",
    )
    model.add_population(
        initial_size=N_ancDW,
        name="ancDW",
        description="Ancestral dog and wolf",
    )
    model.add_population(
        initial_size=N_root,
        name="root",
        description="Ancestral root of dog, wolf, and jackal",
    )

    model.set_migration_rate(rate=m_ISW_BSJ, source="BSJ", dest="ISW")
    model.set_migration_rate(rate=m_BSJ_ISW, source="ISW", dest="BSJ")
    model.set_migration_rate(rate=m_CHW_DNG, source="DNG", dest="CHW")
    model.set_migration_rate(rate=m_DNG_CHW, source="CHW", dest="DNG")
    model.set_migration_rate(rate=m_GLJ_ISW, source="ISW", dest="GLJ")
    model.set_migration_rate(rate=m_ISW_BSJ, source="BSJ", dest="ISW")
    model.set_migration_rate(rate=m_ISW_GLJ, source="GLJ", dest="ISW")
    model.set_migration_rate(rate=m_GLJ_ancDW, source="ancDW", dest="GLJ")
    model.set_migration_rate(rate=m_ancDW_GLJ, source="GLJ", dest="ancDW")

    model.add_population_split(
        time=T_ancDOG1, derived=["BOX", "BSJ"], ancestral="ancDOG1"
    )
    model.add_population_split(
        time=T_ancDOG, derived=["ancDOG1", "DNG"], ancestral="ancDOG"
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
    model.add_population_split(
        time=T_ancroot, derived=["ancDW", "GLJ"], ancestral="root"
    )

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        citations=citations,
        generation_time=generation_time,
        model=model,
    )


_species.add_demographic_model(_dog_wolf_admixture())
