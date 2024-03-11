import stdpopsim

_species = stdpopsim.get_species("AraTha")

_an = stdpopsim.Annotation(
    species=_species,
    id="araport_11_exons",
    description="Araport11 exon annotations on TAIR10",
    url=(
        "http://ftp.ensemblgenomes.org/pub/plants/release-54/"
        "gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.54.gff3.gz"
    ),
    gff_sha256="72dbbe75631b0c42499c65712911dc24c0da4a265275e415c4732ffe21c00425",
    intervals_url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/"
        "annotations/AraTha/araport_11_v1_exons.tar.gz"
    ),
    intervals_sha256="5702bb2aecd3e7345ce40fece6236a2483f17d500e0ceabab237d6210749c87c",
    citations=[
        stdpopsim.Citation(
            year=2017,
            author="Cheng et al",
            doi="http://dx.doi.org/10.1111/tpj.13415",
            reasons={stdpopsim.CiteReason.ANNOTATION},
        )
    ],
    file_pattern="araport_exons_{id}.txt",
    annotation_source="araport11",
    annotation_type="exon",
)
_species.add_annotations(_an)

# add CDS
_an2 = stdpopsim.Annotation(
    species=_species,
    id="araport_11_CDS",
    description="Araport11 exon annotations on TAIR10",
    url=(
        "http://ftp.ensemblgenomes.org/pub/plants/release-54/"
        "gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.54.gff3.gz"
    ),
    gff_sha256="72dbbe75631b0c42499c65712911dc24c0da4a265275e415c4732ffe21c00425",
    intervals_url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/"
        "annotations/AraTha/araport_11_v1_CDS.tar.gz"
    ),
    intervals_sha256="2d5b7c81502aea0791fde1eb2731747713822e8957e56a67ea9f5865fceae070",
    citations=[
        stdpopsim.Citation(
            year=2017,
            author="Cheng et al",
            doi="http://dx.doi.org/10.1111/tpj.13415",
            reasons={stdpopsim.CiteReason.ANNOTATION},
        )
    ],
    file_pattern="araport_CDS_{id}.txt",
    annotation_source="araport11",
    annotation_type="CDS",
)
_species.add_annotations(_an2)
