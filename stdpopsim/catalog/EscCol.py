"""
Genome and demographic model definitions for Escherichia coli.
"""
import stdpopsim

###########################################################
#
# Genome definition
#
###########################################################

_lapierre_et_al = stdpopsim.Citation(
    author="Lapierre et al.",
    year="2016",
    doi="https://doi.org/10.1093/molbev/msw048")

_sezonov_et_al = stdpopsim.Citation(
    author="Sezonov et al.",
    year="2007",
    doi="https://doi.org/10.1128/JB.01368-07")

_perfeito_et_al = stdpopsim.Citation(
    author="Perfeito et al.",
    year="2007",
    doi="https://doi.org/10.1126/science.1142284")

_kibota_and_lynch = stdpopsim.Citation(
    author="Kibota and Lynch",
    year="1996",
    doi="https://doi.org/10.1038/381694a0")

_blattner_et_al = stdpopsim.Citation(
    author="Blattner et al.",
    year="1997",
    doi="10.1126/science.277.5331.1453")

_chromosomes = []
_chromosomes.append(stdpopsim.Chromosome(
        id=None,
        length=4641652,
        # Lapierre et al. (2016) refer to:
        #  Genomic adaptive mutation rate: 1e-5, Perfeito et al. (2007), and
        #  Genomic deleterious mutation rate: 2e−4, Kibota and Lynch (1996).
        mutation_rate=1e-5+2e-4,
        recombination_rate=0.0))
# mean_conversion_rate=8.9e-11 # not implemented yet!
# mean_conversion_length=542 # not implemented yet!

#: :class:`stdpopsim.Genome` definition for E. Coli.
# Chromosome length data is based on strain K-12.

_genome = stdpopsim.Genome(
        chromosomes=_chromosomes,
        mutation_rate_citations=[
            _perfeito_et_al.because(stdpopsim.CiteReason.MUT_RATE),
            _kibota_and_lynch.because(stdpopsim.CiteReason.MUT_RATE),
            ],
        assembly_citations=[
            _blattner_et_al.because(stdpopsim.CiteReason.ASSEMBLY)])

_species = stdpopsim.Species(
    id="EscCol",
    name="Escherichia coli",
    common_name="E. coli",
    genome=_genome,
    generation_time=0.00003805175,  # 1.0 / (525600 min/year / 20 min/gen)
    generation_time_citations=[
        _sezonov_et_al.because(stdpopsim.CiteReason.GEN_TIME)],
    population_size=1.8e8,
    population_size_citations=[
        _lapierre_et_al.because(stdpopsim.CiteReason.POP_SIZE)])

stdpopsim.register_species(_species)
