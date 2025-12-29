import stdpopsim

from . import genome_data

# Recombination rates
# From linkage map in Mas-Gomez et al. (2025) Table 1
# (divided the genetic length per chromosome in cM with physical length)
_recombination_rate = {
    "1": 1.97e-8,
    "2": 1.51e-8,
    "3": 1.89e-8,
    "4": 2.22e-8,
    "5": 2.54e-8,
    "6": 1.98e-8,
    "7": 2.61e-8,
    "8": 2.12e-8,
    "MT": 0,
    "Pt": 0,
}

# Generic and chromosome-specific ploidy
_species_ploidy = 2
_ploidy = {
    "1": _species_ploidy,
    "2": _species_ploidy,
    "3": _species_ploidy,
    "4": _species_ploidy,
    "5": _species_ploidy,
    "6": _species_ploidy,
    "7": _species_ploidy,
    "8": _species_ploidy,
    "MT": 1,
    "Pt": 1,
}


# Evolutionary Genomics of Peach and Almond Domestication
_VelascoEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1534/g3.116.032672",
    year=2016,
    author="Velasco et al.",
    reasons={
        stdpopsim.CiteReason.GEN_TIME,  # page 3987
        stdpopsim.CiteReason.POP_SIZE,  # Figure S8, page 3987-8
        stdpopsim.CiteReason.MUT_RATE,  # page 3987
    },
)

# A phased genome of the highly heterozygous ‘Texas’ almond
# uncovers patterns of allele-specific expression linked to
# heterozygous structural variants
_CastaneraEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1093/hr/uhae106",
    year=2024,
    author="Castanera et al.",
    reasons={stdpopsim.CiteReason.ASSEMBLY},
)

# Integration of linkage mapping, QTL analysis, RNA-Seq data, and
# Genome-Wide Association Studies (GWAS) to explore relative
# flowering traits in almond
_MasGomezEtAl = stdpopsim.Citation(
    doi="https://doi.org/10.1016/j.hpj.2025.04.013",
    year=2025,
    author="Mas-Gomez et al.",
    reasons={stdpopsim.CiteReason.REC_RATE},  # Table 1
)

# Whole-genome sequence and methylome profiling of the almond
# [Prunus dulcis (Mill.) D.A. Webb] cultivar ‘Nonpareil’
_DAmicoWillmanEtal = stdpopsim.Citation(
    doi="https://doi.org/10.1093/g3journal/jkac065",
    year=2022,
    author="D'Amico-Willman et al.",
    reasons={stdpopsim.CiteReason.ASSEMBLY},  # mitochondrial and chloroplast
)

_overall_rate = 1e-8  # per generation, Velasco et al. (2016), page 3987
_mutation_rate = {
    "1": _overall_rate,
    "2": _overall_rate,
    "3": _overall_rate,
    "4": _overall_rate,
    "5": _overall_rate,
    "6": _overall_rate,
    "7": _overall_rate,
    "8": _overall_rate,
    "MT": _overall_rate,
    "Pt": _overall_rate,
}


_genome = stdpopsim.Genome.from_data(
    genome_data.data,
    recombination_rate=_recombination_rate,
    mutation_rate=_mutation_rate,
    ploidy=_ploidy,
    citations=[
        _VelascoEtAl,  # MUT_RATE
        _MasGomezEtAl,  # REC_RATE
        _CastaneraEtAl,  # ASSEMBLY
        _DAmicoWillmanEtal,  # ASSEMBLY
    ],
)
stdpopsim.utils.append_common_synonyms(_genome)

_species = stdpopsim.Species(
    id="PruDul",
    ensembl_id="prunus_dulcis",
    name="Prunus dulcis",
    common_name="Prunus dulcis",
    genome=_genome,
    generation_time=10,  # Velasco et al. (2016), page 3987
    ploidy=_species_ploidy,
    population_size=216627,  # Velasco et al. (2016), Figure S8, page 3987-8
    citations=[
        _VelascoEtAl,  # GEN_TIME, POP_SIZE
    ],
)

stdpopsim.register_species(_species)
