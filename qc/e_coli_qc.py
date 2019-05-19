import msprime
import numpy as np
import math
import stdpopsim.models as models

class LapierreConstant(models.Model):
    def __init__(self):
        super().__init__()
        # From paper https://doi.org/10.1093/molbev/msw048
        # Ne taken from Table 2
        Ne = 1.8e8

        # Single pop with initial size Ne
        self.population_configurations = [
            msprime.PopulationConfiguration(initial_size=Ne),
        ]
        # No migration
        self.migration_matrix = [[0]]
        # No demographic events
        self.demographic_events = []