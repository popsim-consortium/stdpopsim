import stdpopsim
from . import genome_data


_FrezalAndFelix2015 = stdpopsim.Citation(
    author="Frézal & Félix",
    year=2015,
    doi="https://doi.org/10.7554/eLife.05849",
    reasons={stdpopsim.CiteReason.GEN_TIME},
)


_genome1998 = stdpopsim.Citation(
    doi="https://doi.org/10.1126/science.282.5396.2012",
    year=1998,
    author="The C. elegans Sequencing Consortium*",
    reasons={stdpopsim.CiteReason.ASSEMBLY},
)


_Rockman2009 = stdpopsim.Citation(
    doi="https://doi.org/10.1371/journal.pgen.1000419",
    year=2009,
    author="Rockman & Kruglyak",
    reasons={stdpopsim.CiteReason.REC_RATE},
)

_KonradEtAl2019 = stdpopsim.Citation(
    doi="https://doi.org/10.1534/genetics.119.302054",
    year=2019,
    author="Konrad et al.",
    reasons={stdpopsim.CiteReason.MUT_RATE},
)

_KonradEtAl2017 = stdpopsim.Citation(
    author="Konrad et al.",
    doi="https://doi.org/10.1093/molbev/msx051",
    year=2017,
    reasons={stdpopsim.CiteReason.MUT_RATE},
)


_BarriereAndFelix2005 = stdpopsim.Citation(
    doi="https://doi.org/10.1016/j.cub.2005.06.022",
    year=2005,
    author="Barrière & Félix",
    reasons={stdpopsim.CiteReason.POP_SIZE},
)


_Cutter2006 = stdpopsim.Citation(
    doi="https://doi.org/10.1534/genetics.105.048207",
    year=2006,
    author="Cutter",
    reasons={stdpopsim.CiteReason.POP_SIZE},
)

_SivasundarAndHey2003 = stdpopsim.Citation(
    doi="https://doi.org/10.1093/genetics/163.1.147",
    year=2003,
    author="Sivasundar & Hey",
    reasons={stdpopsim.CiteReason.POP_SIZE},
)

_SivasundarAndHey2005 = stdpopsim.Citation(
    doi="https://doi.org/10.1016/j.cub.2005.08.034",
    year=2005,
    author="Sivasundar & Hey",
    reasons={stdpopsim.CiteReason.POP_SIZE},
)

# Mean chromosomal rates, calculated from the
# Rockman and Kruglyak 2009 map. !!! The recombination
# rate is multiplied by 0.001 because the outcrosssing
# rate is 0.1%
_recombination_rate_data = {
    "I": 3.1216265402124167e-11,
    "II": 3.529290802315087e-11,
    "III": 3.906598767640363e-11,
    "IV": 2.712098077556377e-11,
    "V": 2.4705737572511805e-11,
    "X": 2.9472374817864404e-11,
    "MtDNA": 0,
}

genomic_dna_mu = 1.84e-9
_mutation_rate_data = {
    "I": genomic_dna_mu,
    "II": genomic_dna_mu,
    "III": genomic_dna_mu,
    "IV": genomic_dna_mu,
    "V": genomic_dna_mu,
    "X": genomic_dna_mu,
    "MtDNA": 1.05e-7,
}

# Generic and chromosome-specific ploidy
_species_ploidy = 2
_ploidy = {
    "I": _species_ploidy,
    "II": _species_ploidy,
    "III": _species_ploidy,
    "IV": _species_ploidy,
    "V": _species_ploidy,
    "X": _species_ploidy,
    "MtDNA": 1,
}

_genome = stdpopsim.Genome.from_data(
    genome_data=genome_data.data,
    recombination_rate=_recombination_rate_data,
    mutation_rate=_mutation_rate_data,
    ploidy=_ploidy,
    citations=[
        _genome1998,
        _KonradEtAl2019.because(stdpopsim.CiteReason.MUT_RATE),
        _KonradEtAl2017.because(stdpopsim.CiteReason.MUT_RATE),
        _Rockman2009.because(stdpopsim.CiteReason.REC_RATE),
    ],
)

_species = stdpopsim.Species(
    id="CaeEle",
    ensembl_id="caenorhabditis_elegans",
    name="Caenorhabditis elegans",
    common_name="C. elegans",
    genome=_genome,
    generation_time=0.01,  # One of the estimates reported
    # by the paper is ~150 generations per year (0.00666).
    # We settled on 0.01, because generations are expected
    # to be fewer in nature due to dauer larvae stage.
    population_size=10000,  # Lots of estimates available,
    # they are all over the place. Settling on 10,000.
    #    selfing_rate,=0.999,
    ploidy=_species_ploidy,
    citations=[
        _FrezalAndFelix2015.because(stdpopsim.CiteReason.GEN_TIME),
        _BarriereAndFelix2005.because(stdpopsim.CiteReason.POP_SIZE),
        _Cutter2006.because(stdpopsim.CiteReason.POP_SIZE),
    ],
)

stdpopsim.register_species(_species)
