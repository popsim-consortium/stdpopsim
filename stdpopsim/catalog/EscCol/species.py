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
            recombination_rate=8.9e-11,
            gene_conversion_length=542,
        )
    )

_genome = stdpopsim.Genome(
    chromosomes=_chromosomes,
    bacterial_recombination=True,
    assembly_name=genome_data.data["assembly_name"],
    assembly_accession=genome_data.data["assembly_accession"],
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
    citations=[
        _sezonov_et_al.because(stdpopsim.CiteReason.GEN_TIME),
        _hartl_et_al.because(stdpopsim.CiteReason.POP_SIZE),
    ],
)

stdpopsim.register_species(_species)
