import stdpopsim


_hapmap2007 = stdpopsim.Citation(
    doi="https://doi.org/10.1038/nature06258",
    year=2007,
    author="The International HapMap Consortium",
)

_species = stdpopsim.get_species("HomSap")

_gm = stdpopsim.GeneticMap(
    species=_species,
    id="HapMapII_GRCh37",
    description="HapMap Phase II lifted over to GRCh37",
    long_description="""
        This genetic map is from the Phase II Hapmap project
        and based on 3.1 million genotyped SNPs
        from 270 individuals across four populations (YRI, CEU, CHB and JPT).
        Genome wide recombination rates were estimated using LDHat.
        This version of the HapMap genetic map was lifted over to GRCh37
        (and adjusted in regions where the genome assembly had rearranged)
        for use in the 1000 Genomes project. Please see the README file on
        the 1000 Genomes download site for details of these adjustments.
        ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20110106_recombination_hotspots
        """,
    url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/genetic_maps/"
        "HomSap/HapmapII_GRCh37_RecombinationHotspots.tar.gz"
    ),
    sha256="80f22d9e6cb0e497074ed1bc277e765fa9d8e22f21b2f66c3b10286520f6b68f",
    file_pattern="genetic_map_GRCh37_chr{id}.txt",
    citations=[_hapmap2007.because(stdpopsim.CiteReason.GEN_MAP)],
)
_species.add_genetic_map(_gm)

_gm = stdpopsim.GeneticMap(
    species=_species,
    id="HapMapII_GRCh38",
    description="HapMap Phase II lifted over to GRCh37 then lifted over to GRCh38",
    long_description="""
        This genetic map is from the Phase II Hapmap project
        and based on 3.1 million genotyped SNPs
        from 270 individuals across four populations (YRI, CEU, CHB and JPT).
        Genome wide recombination rates were estimated using LDHat.
        This version is lifted over to GRCh38 using liftover from the HapMap Phase II
        map previously lifted over to GRCh37.
        Liftover was performed using the liftOver_catalog.py script from
        stdpopsim/maintainance. Exact command used is as follows:
        `python <path_to_stdpopsim>/stdpopsim/maintenance/liftOver_catalog.py
        --species HomSap --map HapMapII_GRCh37
        --chainFile <path_to_chainfiles>/chainfiles/hg19ToHg38.over.chain.gz
        --validationChain <path_to_chainfiles>/chainfiles/hg38ToHg19.over.chain.gz
        --winLen 1000 --useAdjacentAvg --retainIntermediates --gapThresh 1000000`
        """,
    url=(
        "https://stdpopsim.s3.us-west-2.amazonaws.com/genetic_maps/"
        "HomSap/HapMapII_GRCh38.tar.gz"
    ),
    sha256="497512ed1c0f8a40e9aa13696049a9f8c3cb062e898921cfd7d85ce9d14c4baa",
    file_pattern="genetic_map_Hg38_chr{id}.txt",
    citations=[_hapmap2007.because(stdpopsim.CiteReason.GEN_MAP)],
)
_species.add_genetic_map(_gm)


_gm = stdpopsim.GeneticMap(
    species=_species,
    id="DeCodeSexAveraged_GRCh36",
    description="Sex averaged map from deCode family study",
    long_description="""
        This genetic map is from the deCode study of recombination
        events in 15,257 parent-offspring pairs from Iceland.
        289,658 phased autosomal SNPs were used to call recombinations
        within these families, and recombination rates computed from the
        density of these events. This is the combined male and female
        (sex averaged) map. See
        https://www.decode.com/addendum/ for more details.""",
    url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/genetic_maps/"
        "HomSap/decode_2010_sex-averaged_map.tar.gz"
    ),
    sha256="a17f491756f863c1cc61240c0ce4f383d8085c68b45edd3fed677bb9692af23d",
    file_pattern="genetic_map_decode_2010_sex-averaged_chr{id}.txt",
    citations=[
        stdpopsim.Citation(
            year=2010,
            author="Kong et al",
            doi="https://doi.org/10.1038/nature09525",
            reasons={stdpopsim.CiteReason.GEN_MAP},
        )
    ],
)
_species.add_genetic_map(_gm)


_gm = stdpopsim.GeneticMap(
    species=_species,
    id="DeCodeSexAveraged_GRCh38",
    description="Sex averaged map from deCode family study",
    long_description="""
        This genetic map is from the deCode study of recombination
        events in 15,257 parent-offspring pairs from Iceland.
        289,658 phased autosomal SNPs were used to call recombinations
        within these families, and recombination rates computed from the
        density of these events. This is the combined male and female
        (sex averaged) map. See
        https://www.decode.com/addendum/ for more details.
        This map is further lifted over from the original GRCh36 to
        GRCh38 using liftover.
        Liftover was performed using the liftOver_catalog.py script from
        stdpopsim/maintainance. Exact command used is as follows:
        `python <path_to_stdpopsim>/stdpopsim/maintenance/liftOver_catalog.py
        --species HomSap --map DeCodeSexAveraged_GRCh36
        --chainFile <path_to_chainfiles>/chainfiles/hg18ToHg38.over.chain.gz
        --winLen 1000 --useAdjacentAvg --retainIntermediates --gapThresh 1000000`
        Validation chain file does not exist for liftover between these two
        assemblies hence validation was performed manually separately.
        """,
    url=(
        "https://stdpopsim.s3.us-west-2.amazonaws.com/genetic_maps/"
        "HomSap/DeCodeSexAveraged_GRCh38.tar.gz"
    ),
    sha256="431feece641a720021b1834b97867cbc9842617a11aff6cdd4138af4ac910459",
    file_pattern="genetic_map_Hg38_chr{id}.txt",
    citations=[
        stdpopsim.Citation(
            year=2010,
            author="Kong et al",
            doi="https://doi.org/10.1038/nature09525",
            reasons={stdpopsim.CiteReason.GEN_MAP},
        )
    ],
)
_species.add_genetic_map(_gm)


# 26 populations in 1000 Genomes
for pop, sha256 in zip(
    [
        "ACB",
        "ASW",
        "BEB",
        "CDX",
        "CEU",
        "CHB",
        "CHS",
        "CLM",
        "ESN",
        "FIN",
        "GBR",
        "GIH",
        "GWD",
        "IBS",
        "ITU",
        "JPT",
        "KHV",
        "LWK",
        "MSL",
        "MXL",
        "PEL",
        "PJL",
        "PUR",
        "STU",
        "TSI",
        "YRI",
    ],
    """
588be112b6aac9f82bfb64ea828809051cdf021a627450b2d8c37d7909486903
e7017a82e92e87f443cd21c7f1c3ce24168fa5efaa48f8448c8e095d13aebcaf
c6f69e6a847af6f88fbf1c82024af6fd5eeb35c167e3a721d9de6d6ff146a3ba
37cd4e8ac16775b9c8126dc09ab90840925f47138675d6a5942dd6264d089f6c
fcfacc1858bff5a64480982bab0821d47c4c56262310dd85290f907b8ea6b081
3836a7853b2e128b6bdbeb0e5eb1ba582975304a9efa10b742be856efdc520a9
e0ba174edf944e96f327d799010a5be8ecc110968be6b27fb8f6901d5a413cbb
969f04541d90b70e68f75c0edcfffeb00bee2cd97cf7ad598439f18a8d3a0885
6310f3860e4ccfcedddc173df2a03400d654ddd83c78eb88ec92e104d00bf40b
bb7ac5df072b8b2a4567825c5ee4382b67cae5290f5879a1d7c3f743f244b3f6
798daacb41ae538017303f103c39cdb88add7d88eb964922fe6c81406d8aaa64
28e4552e485b55b0bc392485cd389a155fa83ee3556f221b3b3377d2c996cbdf
f2d57c835ac61763112c15acd9cc8240af7c7fb1d7cfa7543cb65a565da2f62b
bbc49af3646b47dd5b58b8ae28716dcce695eb4a9cb223e261f45345b1a051b9
050b965efb0b4c9197a289521406586f4207dd4833be9fd429d118f7d061788a
e7284218b4bf3e1487b9194c9f243c5dff88f5edbd40247118225a9acc08c76e
74cd0bc6bc12c2b59f928c6a3d60773b6118625375266ffba0c89da61c376f55
dd3ef07641ad7c8d3bfa74a98a5820f1b70ea74c3e82103663bbe61931bf8579
ad5064ee575a10f890c81b9bdbd55128197582db1c9b153b575f5e7fbe542b0d
f37b50b940816d9c12429a4b71709b18569e77c1ea4c78d5fed19141751822c3
024ac242da742d006a210cbf77c5b8457735a6690dc78cc01686deafe7cb7beb
871b518e3b9ca92838e206125e4a1187690018ec5bddeaf39eb2e6f6be5705f1
00f87ca1ad6430b272251d3ce5fa93b44a0a50d6472237bcb5dbad2dbcd567b8
fa9e61492de4c160cf987ac6a8cc3a77efdea751f0f5aa9303a294a37f81cbd7
e3cf06041c6ffceb29fbc617c86966948203e166c3399d036c95915e3d5ebdcd
2989aecaa086af37cdef9ee2245f12000b3f29c08b5d1a71744578a684fe377f
    """.strip().splitlines(),
):
    _gm = stdpopsim.GeneticMap(
        species=_species,
        id="Pyrho{}_GRCh38".format(pop),
        description="Pyrho population-specific map for {}".format(pop),
        long_description="""
        This genetic map was inferred using individuals from the {}
        population from Phase 3 of the 1000 Genomes Project. Rates
        were estimated using pyrho (https://github.com/popgenmethods/pyrho)
        while using population-specific population size history estimates
        obtained from smc++ (https://github.com/popgenmethods/smcpp).
        Genetic maps are only available for the 22 autosomes.
        See https://doi.org/10.1126/sciadv.aaw9206 for more
        details.""".format(
            pop
        ),
        url=(
            "https://stdpopsim.s3-us-west-2.amazonaws.com/genetic_maps/"
            "HomSap/Pyrho{}_GRCh38.tar.gz".format(pop)
        ),
        sha256=sha256,
        file_pattern=("Pyrho{}_GRCh38_chr{{id}}.txt".format(pop)),
        citations=[
            stdpopsim.Citation(
                year=2019,
                author="Spence and Song",
                doi="https://doi.org/10.1126/sciadv.aaw9206",
                reasons={stdpopsim.CiteReason.GEN_MAP},
            )
        ],
    )
    _species.add_genetic_map(_gm)
