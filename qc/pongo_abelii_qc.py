# This script is a QC implementation of the Pongo model
import msprime
import numpy as np
import stdpopsim.models as models

# Some generic populations to use for qc
population_sample_0 = models.Population("sampling_0",
                                        "Population that samples at time 0",
                                        0)


class LockePongo(models.DemographicModel):
    populations = [population_sample_0] * 2

    def __init__(self):
        # This is a split-migration style model, with exponential growth or
        # decay allowed in each population after the split. They assumed a
        # generation time of 20 years and a mutation rate of 2e-8 per bp per gen
        generation_time = 20

        # Parameters given in Table S21-2
        Ne = 17934
        s = 0.592
        NB0 = s*Ne
        NS0 = (1-s)*Ne
        NBF = 8805
        NSF = 37661
        mSB = 0.395 / 2 / Ne
        mBS = 0.239 / 2 / Ne
        T = 403149 / generation_time

        rB = np.log(NBF/NB0) / T
        rS = np.log(NSF/NS0) / T

        # pop 0 is Bornean, pop 1 is Sumatran
        self.population_configurations = [
            msprime.PopulationConfiguration(
                initial_size=NBF, growth_rate=rB),
            msprime.PopulationConfiguration(
                initial_size=NSF, growth_rate=rS)
        ]

        self.migration_matrix = [[0, mBS], [mSB, 0]]

        self.demographic_events = [
            # merge, turn off migration, change size and growth rate
            msprime.MassMigration(
                source=1, destination=0, time=T, proportion=1),
            msprime.MigrationRateChange(
                time=T, rate=0),
            msprime.PopulationParametersChange(
                time=T, initial_size=Ne, growth_rate=0, population_id=0)
        ]
