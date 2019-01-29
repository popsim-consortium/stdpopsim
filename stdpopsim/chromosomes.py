"""
Infrastructure for defining chromosome information for different species.
"""



class Chromosome(object):

    def __init__(self, length, mean_recombination_rate, mean_mutation_rate):
        self.length = length
        self.mean_recombination_rate = mean_recombination_rate
        self.mean_mutation_rate = mean_mutation_rate

    # TODO add methods to return recombination maps


