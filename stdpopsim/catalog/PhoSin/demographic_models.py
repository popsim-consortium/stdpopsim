import msprime
import stdpopsim

_species = stdpopsim.get_species("PhoSin")


def _2epoch():
    N_anc = 4485
    N_curr = 2807
    T = 2162
    populations = [
        stdpopsim.Population(
            id="Vaquita",
            description="Vaquita (Phocoena sinus)",
        )
    ]

    return stdpopsim.DemographicModel(
        id="Vaquita2Epoch_1R22",
        description="Vaquita two epoch model",
        long_description="""
            A two-epoch demographic model estimated using dadi from the site
            frequency spectrum at putatively neutrally evolving regions of
            the genome identified as those located >10 kb from coding
            sequences which did not overlap with CpG islands.
            Population genomic data obtained from 20 individuals sequenced
            at mean coverage 60x.
            Robinson et al. (2022) reports several inferred models in Supp
            Table S2. This is the 2-epoch model inferred by dadi, which is
            also depicted in Main Figure 1E.
            Size changes from N_anc to N_curr in time T.
        """,
        populations=populations,
        citations=[
            stdpopsim.Citation(
                author="Robinson et al.",
                year=2022,
                doi="https://doi.org/10.1126/science.abm1742",
                reasons={stdpopsim.CiteReason.DEM_MODEL},
            )
        ],
        generation_time=11.9,
        mutation_rate=5.83e-9,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_curr, metadata=populations[0].asdict()
            )
        ],
        demographic_events=[
            msprime.PopulationParametersChange(
                time=T, initial_size=N_anc, population_id=0
            )
        ],
    )


_species.add_demographic_model(_2epoch())
