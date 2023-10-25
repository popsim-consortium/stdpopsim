import stdpopsim

_species = stdpopsim.get_species("DroMel")
# exons
_an = stdpopsim.Annotation(
    species=_species,
    id="FlyBase_BDGP6.32.51_exons",
    description="FlyBase exon annotations on BDGP6",
    url=(
        "http://ftp.ensembl.org/pub/release-104/"
        "gff3/drosophila_melanogaster/"
        "Drosophila_melanogaster.BDGP6.32.104.gff3.gz"
    ),
    gff_sha256="d882d9a2af1c090ad69b4c81e54b809506f7a8d5fdd90597c6ed05c79ad502bc",
    intervals_url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/"
        "annotations/DroMel/FlyBase_BDGP6.32.51_exons.tar.gz"
    ),
    intervals_sha256="665df7639725ca18bde85001bbdc08fbef4450c52a8788e02ba76206e2e03d50",
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
        "http://ftp.ensembl.org/pub/release-104/"
        "gff3/drosophila_melanogaster/"
        "Drosophila_melanogaster.BDGP6.32.104.gff3.gz"
    ),
    gff_sha256="d882d9a2af1c090ad69b4c81e54b809506f7a8d5fdd90597c6ed05c79ad502bc",
    intervals_url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/"
        "annotations/DroMel/FlyBase_BDGP6.32.51_CDS.tar.gz"
    ),
    intervals_sha256="3f84b69f7ca558f16cc9b787ee74ba5d16685cb80b2ba8a6ddca4b7b0c4486e0",
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
