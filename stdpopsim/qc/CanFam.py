import msprime

import stdpopsim

_species = stdpopsim.get_species("CanFam")


def FreedmanEarlyWolfAdmixture():
    """
    Demographic model from Freedman et al (2014)
    https://doi.org/10.1371/journal.pgen.1004016
    The model is illustrated in Figure 5 and
    exact parameters given in Table S12.
    Further information is given in Text S8 and S9.
    """
    id = "QC-EarlyWolfAdmixture_6F14"

    generation_interval = 3  # page 3, Figure 5, and Table S12
    # Page 4 in Text S8 mentiones the use of a combination of
    # 2 years (for first 10,000 years) and 3 years (thereafter)
    # for Dingo and Basenji.
    # This was to account for different generation intervals between
    # domesticated dogs and wild/ancestral populations
    # when translating PSMC generations to calendar years.
    # 3 years was used for all other populations.
    # The rest of the text always mentions 3 years,
    # especially Table S12 which gives calibrated parameters
    # from the G-PhoCS analysis.

    mutation_rate = 1e-8  # page 3, Figure 5, Table S12, and page 3 in Text S9
    # Following Lindblad-Toh et al. (2005): _LindbladTohEtAl

    # Recombination rate not specified becase the G-PhoCS analysis
    # was done on "short" and "distant" loci. Meaning we will fallback
    # to species default recombination rate.
    # recombination_rate = ?
    #
    # Page 3 in Text S9 describes the filtering of data for G-PhoCS analysis.
    # They used 1 kb loci to reduce intra-locus recombination
    # in the time scale of dog and wolf evolution.
    # The loci were at least 30 kb apart to get sufficient inter-locus recombination and
    # hence reduce the correlation between the local genealogies at the loci.
    #
    # Additional notes on recombination rate:
    #
    # Page 7 in Text S8 describes MaCS simulations
    # with a model based on G-PhoCS analysis - to examine the
    # effect of gene-flow on estimates of Ne - to explain some
    # observed patterns.
    # This simulation used a recombination rate of 0.92 cM/Mb.
    # This value is the mean recombination rate estimated from
    # dog linkage map based on microsatellites (Wong et al., 2010).
    # Wong et al. (2010) A Comprehensive Linkage Map of the Dog Genome.
    # Genetics 184: 595-U436.
    #
    # Page 8 in Text S9 describes ms simulations to
    # examine effect of intra-locus recombination on estimates.
    # They simulated data under:
    # r = 0.0 cM/Mb,
    # 0.25 cM/Mb (their PSMC estimate), and
    # 0.92 cM/Mb (from Wong et al., 2010).
    # They found that divergence times from G-PhoCS were robust,
    # appart from Ne and divergence time at the root
    # (due to rec events since divergence from the goden jackal).

    # Estimated effective population sizes, from Table S12 (calibrated values, N)
    N_BSJ = 2_639
    N_DNG = 1_914
    N_ISW = 26_092
    N_CRW = 11_427
    N_CHW = 5_426
    N_GLJ = 19_446
    N_ancDOG1 = 793  # (BOX, BSJ)ancDOG1
    N_ancDOG = 1_999  # (DOG, DNG)ancDOG
    N_ancWLF1 = 1_393  # (ISW, CRW)ancWLF1
    N_ancWLF = 12_627  # (WLF1, CHW)WLF
    N_ancDW = 44_993  # (DOG, WLF)ancDW
    N_root = 18_169  # (DW, GLJ)root

    # Estimated times, from Table S12 (calibrated values, years)
    T_ancDOG1 = 12_102  # (BOX, BSJ)ancDOG1
    T_ancDOG = 12_795  # (DOG, DNG)ancDOG
    T_ancWLF1 = 13_389  # (ISW, CRW)ancWLF1
    T_ancWLF = 13_455  # (WLF1, CHW)WLF
    T_ancDW = 14_874  # (DOG, WLF)ancDW
    T_root = 398_262  # (DW, GLJ)root

    # Migration rates (m_FROM_TO), from Table S12
    m_ISW_BSJ = 0.18  # ISW -> BSJ (see Figure 5 to clarify direction)
    m_BSJ_ISW = 0.07  # BSJ -> ISW (see Figure 5 to clarify direction)
    m_CHW_DNG = 0.03
    m_DNG_CHW = 0.04
    m_GLJ_ISW = 0.00
    m_ISW_GLJ = 0.05
    m_GLJ_ancDW = 0.02
    m_ancDW_GLJ = 0.99

    # Demographic model based on the above parameters - populations
    model = msprime.Demography()
    model.add_population(initial_size=N_BSJ, name="BSJ", description="Basenji")
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
    #     description="Ancestral (BOX, BSJ)",
    # )
    # commenting ancDOG1 pop due to the removal of Boxer
    model.add_population(
        initial_size=N_ancDOG,
        name="ancDOG",
        description="Ancestral ((BOX, BSJ), DNG)",
    )
    model.add_population(
        initial_size=N_ancWLF1,
        name="ancWLF1",
        description="Ancestral (ISW, CRW)",
        # Not clear from Figure 5 or Table S12 that
        # ancWLF1 is the ancestor of ISW and CRW, but these two
        # populations are geogprahically closer (Israel and Croatia)
        # than with CHW (China).
    )
    model.add_population(
        initial_size=N_ancWLF,
        name="ancWLF",
        description="Ancestral ((ISW, CRW), CHW)",
    )
    model.add_population(
        initial_size=N_ancDW,
        name="ancDW",
        description="Ancestral (DOG, WLF)",
    )
    model.add_population(
        initial_size=N_root,
        name="root",
        description="Ancestral ((DOG, WLF), GLJ)root",
    )

    # Demographic model based on the above parameters - migrations
    # m_ISW_BSJ is ISW -> BSJ, but msprime does these backward in time,
    # so we flip the direction here
    model.set_migration_rate(source="BSJ", dest="ISW", rate=m_ISW_BSJ)
    model.set_migration_rate(source="ISW", dest="BSJ", rate=m_BSJ_ISW)
    model.set_migration_rate(source="DNG", dest="CHW", rate=m_CHW_DNG)
    model.set_migration_rate(source="CHW", dest="DNG", rate=m_DNG_CHW)
    model.set_migration_rate(source="ISW", dest="GLJ", rate=m_GLJ_ISW)
    model.set_migration_rate(source="GLJ", dest="ISW", rate=m_ISW_GLJ)

    # Demographic model based on the above parameters - splits/coalescences
    # model.add_population_split(
    #     time=T_ancDOG1, derived=["BSJ"], ancestral="ancDOG1"
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
        source="ancDW", dest="GLJ", rate=m_GLJ_ancDW, time=T_ancDW
    )
    model.add_migration_rate_change(
        source="GLJ", dest="ancDW", rate=m_ancDW_GLJ, time=T_ancDW
    )
    model.add_population_split(time=T_root, derived=["ancDW", "GLJ"], ancestral="root")

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=generation_interval,
        mutation_rate=mutation_rate,
        model=model,
    )


_species.get_demographic_model("EarlyWolfAdmixture_6F14").register_qc(
    FreedmanEarlyWolfAdmixture()
)
