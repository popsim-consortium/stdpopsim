import msprime
import numpy as np

import stdpopsim


_species = stdpopsim.get_species("PonAbe")


def LockePongo():
    id = "QC-TwoSpecies_2L11"
    populations = [
        stdpopsim.Population("Bornean", ""),
        stdpopsim.Population("Sumatran", ""),
    ]

    # This is a split-migration style model, with exponential growth or
    # decay allowed in each population after the split. They assumed a
    # generation time of 20 years and a mutation rate of 2e-8 per bp per gen
    generation_time = 20

    # Parameters given in Table S21-2
    Ne = 17934
    s = 0.592
    NB0 = s * Ne
    NS0 = (1 - s) * Ne
    NBF = 8805
    NSF = 37661
    mSB = 0.395 / 2 / Ne
    mBS = 0.239 / 2 / Ne
    T = 403149 / generation_time

    rB = np.log(NBF / NB0) / T
    rS = np.log(NSF / NS0) / T

    return stdpopsim.DemographicModel(
        id=id,
        description=id,
        long_description=id,
        generation_time=generation_time,
        populations=populations,
        # pop 0 is Bornean, pop 1 is Sumatran
        population_configurations=[
            msprime.PopulationConfiguration(initial_size=NBF, growth_rate=rB),
            msprime.PopulationConfiguration(initial_size=NSF, growth_rate=rS),
        ],
        migration_matrix=[[0, mBS], [mSB, 0]],
        demographic_events=[
            # merge, turn off migration, change size and growth rate
            msprime.MassMigration(source=1, destination=0, time=T, proportion=1),
            msprime.MigrationRateChange(time=T, rate=0),
            msprime.PopulationParametersChange(
                time=T, initial_size=Ne, growth_rate=0, population_id=0
            ),
        ],
        population_id_map=[
            {"Bornean": 0, "Sumatran": 1},
            {"Bornean": 0, "Sumatran": 1},
        ],
    )


_species.get_demographic_model("TwoSpecies_2L11").register_qc(LockePongo())
