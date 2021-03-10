import collections

import stdpopsim
from . import genome_data


_LiAndStephan = stdpopsim.Citation(
    author="Li et al.",
    year=2006,
    doi="https://doi.org/10.1371/journal.pgen.0020166",
    reasons={stdpopsim.CiteReason.GEN_TIME, stdpopsim.CiteReason.POP_SIZE},
)

_SchriderEtAl = stdpopsim.Citation(
    author="Schrider et al.",
    year=2013,
    doi="https://doi.org/10.1534/genetics.113.151670",
)

_DosSantosEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1093/nar/gku1099",
    year=2015,
    author="dos Santos et al.",
    reasons={stdpopsim.CiteReason.ASSEMBLY},
)

_genome_wide_estimate = 8.4e-9  # WRONG, underestimate used in S&S!

_recombination_rate_data = collections.defaultdict(lambda: _genome_wide_estimate)
# Set some exceptions for non-recombining chrs.
_recombination_rate_data["Y"] = 0
_recombination_rate_data["mitochondrion_genome"] = 0

_chromosomes = []
for name, data in genome_data.data["chromosomes"].items():
    _chromosomes.append(
        stdpopsim.Chromosome(
            id=name,
            length=data["length"],
            synonyms=data["synonyms"],
            mutation_rate=5.49e-9,  # citation: _SchriderEtAl
            recombination_rate=_recombination_rate_data[name],
        )
    )


_genome = stdpopsim.Genome(
    chromosomes=_chromosomes,
    assembly_name=genome_data.data["assembly_name"],
    assembly_accession=genome_data.data["assembly_accession"],
    citations=[_SchriderEtAl.because(stdpopsim.CiteReason.MUT_RATE), _DosSantosEtAl],
)

_species = stdpopsim.Species(
    id="DroMel",
    ensembl_id="drosophila_melanogaster",
    name="Drosophila melanogaster",
    common_name="D. melanogaster",
    genome=_genome,
    generation_time=0.1,
    population_size=1720600,
    citations=[_LiAndStephan],
)

stdpopsim.register_species(_species)
