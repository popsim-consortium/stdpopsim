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
        "annotations/AraTha/araport_11_exons.tar.gz"
    ),
    intervals_sha256="235105850a365f7f171c459d4c5ee50483e0c17e3b2b0232412d22addca9915f",
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
        "annotations/AraTha/araport_11_CDS.tar.gz"
    ),
    intervals_sha256="c1e1c17a0bf3591e91a4cf85a0ad964d1e9205cc2788c40f855870c589cacca7",
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
