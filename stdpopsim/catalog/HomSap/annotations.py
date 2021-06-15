import stdpopsim

_species = stdpopsim.get_species("HomSap")

_an = stdpopsim.Annotation(
    species=_species,
    id="Ensembl_GRCh38_104_gff3",
    description="Ensembl GFF3 annotations on GRCh38",
    url=(
        "ftp://ftp.ensembl.org/pub/release-104/"
        "gff3/homo_sapiens/Homo_sapiens.GRCh38.104.gff3.gz"
    ),
    gff_sha256="313ad46bd4af78b45b9f5d8407bbcbd3f87f4be0747060e84b3b5eb931530ec1",
    intervals_url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/"
        "annotations/HomSap/HomSap.GRCh38.104.tar.gz"
    ),
    intervals_sha256="1a2de33beaf2dada376654769db8370e649ff6f67dd0ec74287544eb52f850b2",
    citations=[
        stdpopsim.Citation(
            year=2018,
            author="Hunt et al",
            doi="https://doi.org/10.1093/database/bay119",
            reasons={stdpopsim.CiteReason.ANNOTATION},
        )
    ],
    file_pattern="ensembl_havana_exons_{id}.txt",
)
_species.add_annotations(_an)
