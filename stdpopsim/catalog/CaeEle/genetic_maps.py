import stdpopsim

_species = stdpopsim.get_species("CaeEle")

_gm = stdpopsim.GeneticMap(
    species=_species,
    id="RockmanRIAIL_ce11",
    description="Genetic map from recombinant inbred advanced intercross lines",
    long_description="""
        The authors genotyped 1454 nuclear SNP markers in 236 recombinant
        inbred advanced intercross lines (RIAILs).
        The genetic distances were estimated in r/qtl using the Haldane map
        function, treating observed recombination fractions as though they had
        been observed in a backcross. The marker density is sufficiently high
        that the exact form of map function employed has little effect on
        estimated genetic distances. The tip domains of each chromosome were
        defined by including all markers between the chromosome ends and the
        first recombination breakpoint observed in the RIAILs.
        The genetic map corresponds to the assembly ce11 (GCA_000002985.3).
        """,
    url="http://sesame.uoregon.edu/~ateterina/rockman2009_maps.tgz",
    sha256="ef0efa0aec3aa8fcd9800830d83d80bb2525788cf7b28105b7292072fedad8fb",
    file_pattern="genetic_map/C.elegans.Rockman.Kruglyak.2009.{id}.hapmap.txt",
    citations=[
        stdpopsim.Citation(
            author="Rockman & Kruglyak",
            doi="https://doi.org/10.1371/journal.pgen.1000419",
            year=2009,
            reasons={stdpopsim.CiteReason.GEN_MAP},
        )
    ],
)
_species.add_genetic_map(_gm)
