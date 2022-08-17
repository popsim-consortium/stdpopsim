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
    intervals_sha256="f09b4684505c7c8c86c7739d632c79927b8d329fd60fd55ea3f610b944bb5856",
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
    intervals_sha256="b993c8fc997e1c7ecdad626b7eeceae724cc0e0e477d8ab2f186866a6a0def15",
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
