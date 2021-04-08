import stdpopsim

_species = stdpopsim.get_species("CanFam")

_CampbellEtAl = stdpopsim.Citation(
    # A Pedigree-Based Map of Recombination in the Domestic Dog Genome.
    author="Campbell et al.",
    year=2016,
    doi="https://doi.org/10.1534/g3.116.034678",
)

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
    sha256="585afb424615e2fb0825d807db0b10fe1c797a6dbb804ecbb3fef5e8387d194f",
    file_pattern="chr{id}_average_canFam3.1.txt",
    citations=[_CampbellEtAl.because(stdpopsim.CiteReason.GEN_MAP)],
)
_species.add_genetic_map(_gm)
