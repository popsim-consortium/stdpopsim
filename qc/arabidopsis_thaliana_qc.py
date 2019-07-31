import msprime
import numpy as np
import math
import stdpopsim.models as models

class Durvasula2017MSMC(models.Model):
    def __init__(self):
        super().__init__()

        # Both of the following are directly
        # converted from MSMC output scaled by A.Thaliana 
        # mutation rate 7e-9 and 1 generation
        # per year.

        self.times = np.array([
            6.990000e+02, 2.796000e+03, 6.068000e+03, 9.894000e+03,
            1.437000e+04, 1.960600e+04, 2.573000e+04, 3.289400e+04,
            4.127500e+04, 5.107700e+04, 6.254400e+04, 7.595800e+04,
            9.164800e+04, 1.100010e+05, 1.314710e+05, 1.565840e+05,
            1.859600e+05, 2.203240e+05, 2.605200e+05, 3.075400e+05,
            3.625410e+05, 4.268790e+05, 5.021390e+05, 5.901730e+05,
            6.931510e+05, 8.136100e+05, 9.545170e+05, 1.119341e+06,
            1.312147e+06, 1.537686e+06, 1.801500e+06, 2.110100e+06])

        self.sizes = np.array([
            4.2252426e+07, 4.2252426e+07, 6.0323000e+04, 7.2174000e+04,
            4.0591000e+04, 2.1158000e+04, 2.1442000e+04, 3.9942000e+04,
            7.8908000e+04, 1.1113200e+05, 1.1074500e+05, 9.6283000e+04,
            8.7661000e+04, 8.3932000e+04, 8.3829000e+04, 9.1813000e+04,
            1.1164400e+05, 1.4345600e+05, 1.8157100e+05, 2.1733100e+05,
            2.4140000e+05, 2.4698400e+05, 2.3859300e+05, 2.2822200e+05,
            2.1775200e+05, 1.9801900e+05, 1.6521000e+05, 1.2179600e+05,
            1.2179600e+05, 7.3989000e+04, 7.3989000e+04, 7.3989000e+04])

        # The first 8 epochs are "masked" to 
        # the last Ne at 40kya due to 
        # the limitations of MSMC to infer 
        # population size in this range.
        #
        # Similarly, the last 2 entries 
        # are set to equal the third last. 
        #
        # Durvasula et al 2017 shows that
        # MSMC has power in A.Thaliana
        # between 40kya and 1.6Mya.

        self.sizes[:8] = self.sizes[8]
        self.sizes[30:32] = self.sizes[30]

        self.generation_times = 1.0
        self.demographic_events = []

        self.population_configurations = [
            msprime.PopulationConfiguration(initial_size=self.sizes[0])
        ]

        for i, t in enumerate(self.times):
            self.demographic_events.append(
                msprime.PopulationParametersChange(
                    time=t, initial_size=self.sizes[i], population_id=0))

        self.migration_matrix = [[0]]



            
