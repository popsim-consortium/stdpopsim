"""
Common infrastructure for specifying demographic models.
"""

import msprime


class Model(object):
    """
    Class representing a simulation model that can be run in msprime.
    """
    def __init__(self):
        self.population_configurations = None
        self.migration_matrix = None
        self.demographic_events = None

    def debug(self):
        # Use the demography debugger to print out the demographic history
        # that we have just described.
        dd = msprime.DemographyDebugger(
            population_configurations=self.population_configurations,
            migration_matrix=self.migration_matrix,
            demographic_events=self.demographic_events)
        dd.print_history()

    def asdict(self):
        return {
            "population_configurations": self.population_configurations,
            "migration_matrix": self.migration_matrix,
            "demographic_events": self.demographic_events}
