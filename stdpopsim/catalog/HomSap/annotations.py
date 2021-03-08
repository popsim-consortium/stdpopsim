import stdpopsim

_species = stdpopsim.get_species("HomSap")

_an = stdpopsim.Annotation(
    species=_species,
    id="Ensembl_GRCh38_gff3",
    description="Ensembl GFF3 annotations on GRCh38",
    url=(
        "ftp://ftp.ensembl.org/pub/release-101/"
        "gff3/homo_sapiens/Homo_sapiens.GRCh38.101.gff3.gz"
    ),
    zarr_url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/annotations/HomSap.GRCh38.zip"
    ),
    zarr_sha256="929c52dc5e5172d4a96a776d066f014f8f207f1477a21f7b88626f81a7142d0b",
    citations=[
        stdpopsim.Citation(
            year=2018,
            author="Hunt et al",
            doi="https://doi.org/10.1093/database/bay119",
            reasons={stdpopsim.CiteReason.ANNOTATION},
        )
    ],
)
_species.add_annotations(_an)
