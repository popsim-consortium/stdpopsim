import stdpopsim

_species = stdpopsim.get_species("HomSap")

_an = stdpopsim.Annotation(
    species=_species,
    id="ensembl_havana_104_exons",
    description="Ensembl Havana exon annotations on GRCh38",
    url=(
        "ftp://ftp.ensembl.org/pub/release-104/"
        "gff3/homo_sapiens/Homo_sapiens.GRCh38.104.gff3.gz"
    ),
    gff_sha256="313ad46bd4af78b45b9f5d8407bbcbd3f87f4be0747060e84b3b5eb931530ec1",
    intervals_url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/"
        "annotations/HomSap/ensembl_havana_104_exons.tar.gz"
    ),
    intervals_sha256="b9df30b3a37cdd26ec625fd80968ff2e412810c045a11dddc458dc606c702c96",
    citations=[
        stdpopsim.Citation(
            year=2018,
            author="Hunt et al",
            doi="https://doi.org/10.1093/database/bay119",
            reasons={stdpopsim.CiteReason.ANNOTATION},
        )
    ],
    file_pattern="ensembl_havana_exons_{id}.txt",
    annotation_source="ensembl_havana",
    annotation_type="exon",
)
_species.add_annotations(_an)

# add CDS
_an2 = stdpopsim.Annotation(
    species=_species,
    id="ensembl_havana_104_CDS",
    description="Ensembl Havana CDS annotations on GRCh38",
    url=(
        "ftp://ftp.ensembl.org/pub/release-104/"
        "gff3/homo_sapiens/Homo_sapiens.GRCh38.104.gff3.gz"
    ),
    gff_sha256="313ad46bd4af78b45b9f5d8407bbcbd3f87f4be0747060e84b3b5eb931530ec1",
    intervals_url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/"
        "annotations/HomSap/ensembl_havana_104_CDS.tar.gz"
    ),
    intervals_sha256="237801b39642b91e733e155dc62960c0b097d3c144eea408327e9dc5f1ac84ae",
    citations=[
        stdpopsim.Citation(
            year=2018,
            author="Hunt et al",
            doi="https://doi.org/10.1093/database/bay119",
            reasons={stdpopsim.CiteReason.ANNOTATION},
        )
    ],
    file_pattern="ensembl_havana_CDS_{id}.txt",
    annotation_source="ensembl_havana",
    annotation_type="CDS",
)
_species.add_annotations(_an2)
