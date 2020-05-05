import stdpopsim

# Lengths are based on CanFam3.1, and recombination rates are the means
# calculated from the Campbell2016 map. ChrX isn't included in that genetic
# map, so the weighted genome-wide average is provided here instead.
_chromosome_data = """\
chr1	122678785	7.636001498077e-09
chr2	85426708	8.798516284201015e-09
chr3	91889043	8.000868549297668e-09
chr4	88276631	8.052302321538563e-09
chr5	88915250	9.344333839828404e-09
chr6	77573801	8.192191165856811e-09
chr7	80974532	7.293465322857858e-09
chr8	74330416	8.291312822214086e-09
chr9	61074082	9.287722798450121e-09
chr10	69331447	9.107151930130527e-09
chr11	74389097	7.639449804041734e-09
chr12	72498081	7.761061128067527e-09
chr13	63241923	8.413015225299971e-09
chr14	60966679	9.028123091776737e-09
chr15	64190966	7.856754982030383e-09
chr16	59632846	8.614063697877811e-09
chr17	64289059	9.718834385033715e-09
chr18	55844845	1.029925446520415e-08
chr19	53741614	1.0425051705950442e-08
chr20	58134056	9.990971109010329e-09
chr21	50858623	1.0339033506095636e-08
chr22	61439934	8.61504863805126e-09
chr23	52294480	9.1266350068389e-09
chr24	47698779	1.1146043069685935e-08
chr25	51628933	1.1543746646856006e-08
chr26	38964690	1.2084571786110922e-08
chr27	45876710	1.1260328390864541e-08
chr28	41182112	1.2463587044656355e-08
chr29	41845238	1.1013629677730708e-08
chr30	40214260	1.1687642441683382e-08
chr31	39895921	1.1397713284329192e-08
chr32	38810281	1.1555927931648279e-08
chr33	31377067	1.3339402745926785e-08
chr34	42124431	1.0483812411227089e-08
chr35	26524999	1.4299102611645524e-08
chr36	30810995	1.187517782077471e-08
chr37	30902991	1.3834580623461596e-08
chr38	23914537	1.4363726512881696e-08
chrX	123869142	9.506483722244087e-09
"""

_LindbladTohEtAl = stdpopsim.Citation(
        # Genome sequence, comparative analysis and haplotype structure of the
        # domestic dog.
        author="Lindblad-Toh et al.",
        year=2005,
        doi="https://doi.org/10.1038/nature04338")

_SkoglundEtAl = stdpopsim.Citation(
        # Ancient wolf genome reveals an early divergence of domestic dog
        # ancestors and admixture into high-latitude breeds.
        author="Skoglund et al.",
        year=2015,
        doi="https://doi.org/10.1016/j.cub.2015.04.019")

_FranzEtAl = stdpopsim.Citation(
        # Genomic and archaeological evidence suggest a dual origin of
        # domestic dogs.
        author="Franz et al.",
        year=2016,
        doi="https://doi.org/10.1126/science.aaf3161")

_CampbellEtAl = stdpopsim.Citation(
        # A Pedigree-Based Map of Recombination in the Domestic Dog Genome.
        author="Campbell et al.",
        year=2016,
        doi="https://doi.org/10.1534/g3.116.034678")

_chromosomes = []
for line in _chromosome_data.splitlines():
    name, length, mean_rr = line.split()
    _chromosomes.append(stdpopsim.Chromosome(
        id=name, length=int(length),
        mutation_rate=4e-9,  # based on non-CpG sites only
        recombination_rate=float(mean_rr)))

_genome = stdpopsim.Genome(
    chromosomes=_chromosomes,
    mutation_rate_citations=[
        _SkoglundEtAl.because(stdpopsim.CiteReason.MUT_RATE),
        _FranzEtAl.because(stdpopsim.CiteReason.MUT_RATE),
        ],
    recombination_rate_citations=[
        _CampbellEtAl.because(stdpopsim.CiteReason.REC_RATE)
        ],
    assembly_citations=[
        _LindbladTohEtAl.because(stdpopsim.CiteReason.ASSEMBLY)
        ],
    )

_species = stdpopsim.Species(
    id="CanFam",
    name="Canis familiaris",
    common_name="Dog",
    genome=_genome,
    generation_time=3,
    generation_time_citations=[
        # Everyone uses 3 years because everyone else uses it.
        # It's likely higher, at least in wolves:
        # https://pubs.er.usgs.gov/publication/70187564
        ],
    population_size=13000,  # ancestral dog size
    population_size_citations=[
        _LindbladTohEtAl.because(stdpopsim.CiteReason.POP_SIZE)
        ],
    )

stdpopsim.register_species(_species)

_gm = stdpopsim.GeneticMap(
    species=_species,
    id="Campbell2016_CanFam3_1",
    description="Pedigree-based crossover map from 237 individuals",
    long_description="""
        Sex-averaged crossover frequency map based on 163,400 autosomal SNPs
        genotyped in a pedigree of 237 Labrador Retriever x Greyhound crosses.
        Genotypes were phased without respect to the pedigree, using SHAPEIT2,
        recombinations were called using duoHMM, and genetic distances were
        obtained using Haldane's map function.
        """,
    url="https://stdpopsim.s3-us-west-2.amazonaws.com/genetic_maps/"
        "CanFam/dog_genetic_maps.tar.gz",
    file_pattern="{id}_average_canFam3.1.txt",
    citations=[
            _CampbellEtAl.because(stdpopsim.CiteReason.GEN_MAP)
        ],
    )
_species.add_genetic_map(_gm)
