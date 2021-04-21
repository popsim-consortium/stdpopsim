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

# Mean chromosomal rates, calculated from the Rockman and Kruglyak 2009 map.
# !!! The recombination rate is multiplied by 0.001 because the outcrosssing rate is 0.1%
_recombination_rate_data = {
    "I": 3.330038e-11,
    "II": 3.999342e-11,
    "III": 4.484974e-11,
    "IV": 2.417689e-11,
    "V": 2.722476e-11,
    "X": 3.447911e-11,
    "MtDNA": 0,
}


_chromosomes = []
for name, data in genome_data.data["chromosomes"].items():
    _chromosomes.append(
        stdpopsim.Chromosome(
            id=name,
            length=data["length"],
            synonyms=data["synonyms"],
            mutation_rate=1.84e-9,  # _Konrad et al. de-nove mutation rate,
            # it's not uniform and it's much better to use a mutation map.
            # mutation_rate=_mutation_rate_data[name],
            recombination_rate=_recombination_rate_data[name],
        )
    )


_genome = stdpopsim.Genome(
    chromosomes=_chromosomes,
    assembly_name=genome_data.data["assembly_name"],
    assembly_accession=genome_data.data["assembly_accession"],
    citations=[
        _genome1998,
        _KonradEtAl2019.because(stdpopsim.CiteReason.MUT_RATE),
        _KonradEtAl2017.because(stdpopsim.CiteReason.MUT_RATE),
        _Rockman2009.because(stdpopsim.CiteReason.REC_RATE),
    ],
)

_species = stdpopsim.Species(
    id="CaeEle",
    ensembl_id="",
    name="Caenorhabditis elegans",
    common_name="C. elegans",
    genome=_genome,
    generation_time=0.01,  # the generation time in the lab ~150
    # generation per year (0.00666), it should be less in the wild
    population_size=10000,
    #    selfing_rate,=0.999,
    citations=[
        _FrezalAndFelix2015.because(stdpopsim.CiteReason.GEN_TIME),
        _BarriereAndFelix2005.because(stdpopsim.CiteReason.POP_SIZE),
        _Cutter2006.because(stdpopsim.CiteReason.POP_SIZE),
    ],
)

stdpopsim.register_species(_species)
