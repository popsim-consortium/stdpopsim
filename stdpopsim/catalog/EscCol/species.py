import stdpopsim
from . import genome_data

_hartl_et_al = stdpopsim.Citation(
    author="Hartl, Moriyama, and Sawyer",
    year=1994,
    # doesn't have a doi
    doi="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1206133/",
)

_sezonov_et_al = stdpopsim.Citation(
    author="Sezonov et al.", year=2007, doi="https://doi.org/10.1128/JB.01368-07"
)

_wielgoss_et_al = stdpopsim.Citation(
    author="Wielgoss et al.", year=2011, doi="https://doi.org/10.1534/g3.111.000406"
)

_blattner_et_al = stdpopsim.Citation(
    author="Blattner et al.",
    year=1997,
    doi="https://doi.org/10.1126/science.277.5331.1453",
)

_didelot_et_al = stdpopsim.Citation(
    author="Didelot et al.",
    year=2012,
    doi="https://doi.org/10.1186/1471-2164-13-256",
)

_species_ploidy = 1
_chromosomes = []
for name, data in genome_data.data["chromosomes"].items():
    _chromosomes.append(
        stdpopsim.Chromosome(
            id=name,
            length=data["length"],
            synonyms=data["synonyms"],
            # Wielgoss et al. (2011) calculated for strain REL606,
            # from synonymous substitutions over 40,000 generations.
            mutation_rate=8.9e-11,
            ploidy=_species_ploidy,
            recombination_rate=8.9e-11,
            gene_conversion_length=542,
        )
    )


_ploidy_data = {str(i.id): _species_ploidy for i in _chromosomes}

_mutation_rate = 8.9e-11
_mutation_rate_data = {str(i.id): _mutation_rate for i in _chromosomes}

_recombination_rate = 8.9e-11
_recombination_rate_data = {str(i.id): _recombination_rate for i in _chromosomes}

_gene_conversion_length = 542
_gene_conversion_data = {str(i.id): _gene_conversion_length for i in _chromosomes}
_genome = stdpopsim.Genome.from_data(
    genome_data=genome_data.data,
    mutation_rate=_mutation_rate_data,
    recombination_rate=_recombination_rate_data,
    gene_conversion_length=_gene_conversion_data,
    bacterial_recombination=True,
    ploidy=_ploidy_data,
    citations=[
        _wielgoss_et_al.because(
            {stdpopsim.CiteReason.MUT_RATE, stdpopsim.CiteReason.GENE_CONVERSION}
        ),
        _blattner_et_al.because(stdpopsim.CiteReason.ASSEMBLY),
        _didelot_et_al.because(stdpopsim.CiteReason.GENE_CONVERSION),
    ],
)
stdpopsim.utils.append_common_synonyms(_genome)

_species = stdpopsim.Species(
    id="EscCol",
    name="Escherichia coli",
    common_name="E. coli",
    # We use the K-12 strain, because the parameters we're using more
    # closely match this strain than the ensembl default (HUSEC2011).
    ensembl_id="escherichia_coli_str_k_12_substr_mg1655_gca_000005845",
    genome=_genome,
    # E. coli K-12 strain MG1655 "doubling time during steady-state growth in
    # Luria-Bertani broth was 20 min".
    generation_time=0.00003805175,  # 1.0 / (525600 min/year / 20 min/gen)
    # Hartl et al. calculated Ne for "natural isolates of E. coli",
    # assuming mu=5e-10 (from Drake 1991).
    population_size=1.8e8,
    ploidy=_species_ploidy,
    citations=[
        _sezonov_et_al.because(stdpopsim.CiteReason.GEN_TIME),
        _hartl_et_al.because(stdpopsim.CiteReason.POP_SIZE),
    ],
)

stdpopsim.register_species(_species)
