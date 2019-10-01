"""
Infrastructure for defining information about species' genomes.
"""
import logging

import attr


logger = logging.getLogger(__name__)


class Genome(object):
    """
    Class representing the genome for a species.

    .. todo:: Define the facilities that this object provides.
    """
    def __init__(self, chromosomes):
        self.chromosomes = chromosomes
        self.length = 0
        for chromosome in chromosomes:
            self.length += chromosome.length

    def __str__(self):
        s = "Chromosomes:\n"
        length_sorted = sorted(self.chromosomes, key=lambda x: -x.length)
        for chrom in length_sorted:
            s += "\t{}\n".format(chrom)
        return s

    def get_chromosome(self, name):
        """
        Returns the chromosome with the specified name.
        """
        for chrom in self.chromosomes:
            if chrom.name == name:
                return chrom
        raise ValueError("Chromosome not found")

    @property
    def mean_recombination_rate(self):
        """
        The length-weighted mean recombination rate across all chromosomes.
        """
        mean_recombination_rate = 0
        for chrom in self.chromosomes:
            normalized_weight = chrom.length / self.length
            cont = chrom.recombination_rate * normalized_weight
            mean_recombination_rate += cont
        return mean_recombination_rate

    @property
    def mean_mutation_rate(self):
        """
        The length-weighted mean mutation rate across all chromosomes.
        """
        mean_mutation_rate = 0
        for chrom in self.chromosomes:
            normalized_weight = chrom.length / self.length
            cont = chrom.mutation_rate * normalized_weight
            mean_mutation_rate += cont
        return mean_mutation_rate


@attr.s(frozen=True)
class Chromosome(object):
    """
    Class representing a single chromosome for a species.

    .. todo:: Define the facilities that this object provides.
    """
    name = attr.ib(type=str, kw_only=True)
    length = attr.ib(kw_only=True)
    recombination_rate = attr.ib(type=float, kw_only=True)
    mutation_rate = attr.ib(type=float, kw_only=True)


@attr.s(frozen=True)
class Contig(object):
    """
    Class representing a contiguous region of genome that is to be
    simulated. This contains all information required to simulate
    the region, including mutation, recombination rates, etc.


    .. todo:: Document properly.
    """
    recombination_map = attr.ib(default=None, kw_only=True)
    mutation_rate = attr.ib(default=None, type=float, kw_only=True)

    def __str__(self):
        s = (
            "Contig(length={:.2G}, recombination_rate={:.2G}, "
            "mutation_rate={:.2G})").format(
                self.recombination_map.get_length(),
                self.recombination_map.mean_recombination_rate,
                self.mutation_rate)
        return s
