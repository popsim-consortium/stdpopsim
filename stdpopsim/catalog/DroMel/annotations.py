import stdpopsim

_species = stdpopsim.get_species("DroMel")
# exons
_an = stdpopsim.Annotation(
    species=_species,
    id="FlyBase_BDGP6.32.51_exons",
    description="FlyBase exon annotations on BDGP6",
    url=(
        "http://ftp.ebi.ac.uk/ensemblgenomes/pub/current/metazoa/"
        "gff3/drosophila_melanogaster/"
        "Drosophila_melanogaster.BDGP6.32.51.gff3.gz"
    ),
    gff_sha256="d882d9a2af1c090ad69b4c81e54b809506f7a8d5fdd90597c6ed05c79ad502bc",
    intervals_url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/"
        "annotations/DroMel/FlyBase_BDGP6.32.51_exons.tar.gz"
    ),
    intervals_sha256="680ee9b0a565f85c561cd3672927cb2bc8649405a750532f91e1717b3ed8b993",
    citations=[
        stdpopsim.Citation(
            year=2014,
            author="Hoskins et al",
            doi="https://doi.org/10.1101/gr.185579.114",
            reasons={stdpopsim.CiteReason.ANNOTATION},
        )
    ],
    file_pattern="flybase_exons_{id}.txt",
    annotation_source="FlyBase",
    annotation_type="exon",
)
_species.add_annotations(_an)

_an2 = stdpopsim.Annotation(
    species=_species,
    id="FlyBase_BDGP6.32.51_CDS",
    description="FlyBase CDS annotations on BDGP6",
    url=(
        "http://ftp.ebi.ac.uk/ensemblgenomes/pub/current/metazoa/"
        "gff3/drosophila_melanogaster/"
        "Drosophila_melanogaster.BDGP6.32.51.gff3.gz"
    ),
    gff_sha256="d882d9a2af1c090ad69b4c81e54b809506f7a8d5fdd90597c6ed05c79ad502bc",
    intervals_url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/"
        "annotations/DroMel/FlyBase_BDGP6.32.51_CDS.tar.gz"
    ),
    intervals_sha256="5f202f454e7d0051f863a5146e775c85e4f3e39bab434cc96701d46026eb7364",
    citations=[
        stdpopsim.Citation(
            year=2014,
            author="Hoskins et al",
            doi="https://doi.org/10.1101/gr.185579.114",
            reasons={stdpopsim.CiteReason.ANNOTATION},
        )
    ],
    file_pattern="flybase_CDS_{id}.txt",
    annotation_source="FlyBase",
    annotation_type="CDS",
)
_species.add_annotations(_an2)
