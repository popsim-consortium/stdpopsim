import msprime
import stdpopsim

demography = msprime.Demography()

_species = stdpopsim.get_species("ApiMel")


_WallbergEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1038/ng.3077",
    year=2014,
    author="Wallberg et al.",
    reasons={stdpopsim.CiteReason.DEM_MODEL},
)


def model():
    id = "10subspecies"
    description = (
        "Ten subspecies within four lineages with a common ancestor in Northwest Africa"
    )
    long_description = """
    Demographic model modelled by Wallberg et al., 2014.
    """
    citations = [_WallbergEtAl.because(stdpopsim.CiteReason.DEM_MODEL)]

    demography = msprime.Demography()

    # Parameters from the paper
    # Ancestral size, before the split of the lineages
    Na = 3.48e05  # TODO: THIS IS STILL A MEAN
    ##################################
    # --- First split: lineages
    ##################################
    # Time of split of the lineages
    Ty_lineageSplit = 3e05
    # Generation time
    generation_time = 2  # TODO: Discuss this value!!!
    # Time of split of the lineages in generations
    Tg_lineageSplit = Ty_lineageSplit / generation_time

    # Size of the populations after the split
    N_A = 500184
    N_M = 210043
    # N_CO =

    # Further split into C and O
    # Ty_COSplit = 1.65e05
    # Tg_COSplit = Ty_COSplit / generation_time TODO: model the CO split
    N_C = 189552
    N_O = 298263

    ##################################
    # --- Second split: subspecies
    ##################################
    # Time of the split of the subspecies in years
    Ty_A1Split = 32709
    Ty_A2Split = 23187
    Ty_M1Split = 37845
    Ty_M2Split = 13069
    Ty_CSplit = 24633
    Ty_OSplit = 80646

    # Time of the split of the subspecies in generations
    Tg_A1Split = Ty_A1Split / generation_time
    Tg_A2Split = Ty_A2Split / generation_time
    Tg_M1Split = Ty_M1Split / generation_time
    Tg_M2Split = Ty_M2Split / generation_time
    Tg_CSplit = Ty_CSplit / generation_time
    Tg_OSplit = Ty_OSplit / generation_time

    # Sizes of the subspecie populations after the split
    N_A_scutellata = 400005
    N_A_capensis = 418821
    N_A_adansonii = 457253
    N_M_melliferaN = 157598
    N_M_melliferaS = 177484
    N_M_iberiensis = 217881
    N_C_carnica = 168783
    N_C_ligustica = 174353
    N_O_syriaca = 313262
    N_O_anatoliaca = 191419

    # Add the ancestral population
    demography.add_population(
        name="Ancestral", initial_size=Na, extra_metadata={id: "Ancestral"}
    )
    # Add the four lineages
    demography.add_population(
        name="A_lineage", initial_size=N_A, extra_metadata={id: "A_lineage"}
    )
    demography.add_population(
        name="M_lineage", initial_size=N_M, extra_metadata={id: "M_lineage"}
    )
    demography.add_population(
        name="C_lineage", initial_size=N_C, extra_metadata={id: "C_lineage"}
    )
    demography.add_population(
        name="O_lineage", initial_size=N_O, extra_metadata={id: "O_lineage"}
    )
    # Add the two populations for the internal splits in A and M lineage
    demography.add_population(
        name="A_capensis_scutellata",
        initial_size=N_A,
        extra_metadata={id: "A_capensis_scutellata"},
    )
    demography.add_population(
        name="M_melliferaN_melliferaS",
        initial_size=N_M,
        extra_metadata={id: "M_melliferaN_melliferaS"},
    )
    # Add the ten subspecies
    demography.add_population(
        name="A_adansonii",
        initial_size=N_A_adansonii,
        extra_metadata={id: "A_adansonii"},
    )
    demography.add_population(
        name="A_capensis", initial_size=N_A_capensis, extra_metadata={id: "A_capensis"}
    )
    demography.add_population(
        name="A_scutellata",
        initial_size=N_A_scutellata,
        extra_metadata={id: "A_scutellata"},
    )
    demography.add_population(
        name="M_iberiensis",
        initial_size=N_M_iberiensis,
        extra_metadata={id: "M_iberiensis"},
    )
    demography.add_population(
        name="M_melliferaN",
        initial_size=N_M_melliferaN,
        extra_metadata={id: "M_melliferaN"},
    )
    demography.add_population(
        name="M_melliferaS",
        initial_size=N_M_melliferaS,
        extra_metadata={id: "M_melliferaS"},
    )
    demography.add_population(
        name="C_carnica", initial_size=N_C_carnica, extra_metadata={id: "C_carnica"}
    )
    demography.add_population(
        name="C_ligustica",
        initial_size=N_C_ligustica,
        extra_metadata={id: "C_ligustica"},
    )
    demography.add_population(
        name="O_syriaca", initial_size=N_O_syriaca, extra_metadata={id: "O_syriaca"}
    )
    demography.add_population(
        name="O_anatoliaca",
        initial_size=N_O_anatoliaca,
        extra_metadata={id: "O_anatoliaca"},
    )

    # Now work from the bottom up
    # First the split of the sublineages
    # A lineage
    # IT DOES NOT WORK WITHOUT THE INTERMEDIATE STEP
    demography.add_population_split(
        time=Tg_A2Split,
        derived=["A_capensis", "A_scutellata"],
        ancestral="A_capensis_scutellata",
    )
    demography.add_population_split(
        time=Tg_A1Split,
        derived=["A_capensis_scutellata", "A_adansonii"],
        ancestral="A_lineage",
    )

    # M lineage
    demography.add_population_split(
        time=Tg_M2Split,
        derived=["M_melliferaN", "M_melliferaS"],
        ancestral="M_melliferaN_melliferaS",
    )
    demography.add_population_split(
        time=Tg_M1Split,
        derived=["M_melliferaN_melliferaS", "M_iberiensis"],
        ancestral="M_lineage",
    )

    # C lineage
    demography.add_population_split(
        time=Tg_CSplit, derived=["C_carnica", "C_ligustica"], ancestral="C_lineage"
    )

    # O lineage
    demography.add_population_split(
        time=Tg_OSplit, derived=["O_syriaca", "O_anatoliaca"], ancestral="O_lineage"
    )

    # Merge the four lineages into the ancestral population
    demography.add_population_split(
        time=Tg_lineageSplit,
        derived=["A_lineage", "M_lineage", "C_lineage", "O_lineage"],
        ancestral="Ancestral",
    )

    # Sort the demographic events
    demography.sort_events()

    # Create a demographic model
    demographicModel = stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        citations=citations,
        generation_time=generation_time,
        model=demography,
    )

    return demographicModel


_species.add_demographic_model(model())
