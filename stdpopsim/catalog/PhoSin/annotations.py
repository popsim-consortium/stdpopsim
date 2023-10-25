import stdpopsim

_species = stdpopsim.get_species("PhoSin")

_an = stdpopsim.Annotation(
    species=_species,
    id="Phocoena_sinus.mPhoSin1.pri.110_exons",
    description="Vaquita exon annotations from Morin et al 2020 using NCBI's pipeline",
    url=(
        "https://ftp.ensembl.org/pub/release-110/"
        "gff3/phocoena_sinus/Phocoena_sinus.mPhoSin1.pri.110.chr.gff3.gz"
    ),
    gff_sha256="dbc8640f8f2e046d063c45f160f39edcf729fc39549726603e389cb7b6ce542e",
    intervals_url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/"
        "annotations/PhoSin/Phocoena_sinus.mPhoSin1.pri.110_exons.tar.gz"
    ),
    intervals_sha256="ec95fd1418b8a9a53cd6486fd6a37e699ecd1fce2635df85d8ea4360d53bfbf4",
    citations=[
        stdpopsim.Citation(
            year=2020,
            author="Morin et al",
            doi="https://doi.org/10.1111/1755-0998.13284",
            reasons={stdpopsim.CiteReason.ANNOTATION},
        )
    ],
    file_pattern="Phocoena_sinus.mPhoSin1.pri.110_exon_{id}.txt",
    annotation_source="ensembl",
    annotation_type="exon",
)
_species.add_annotations(_an)

# add CDS
_an2 = stdpopsim.Annotation(
    species=_species,
    id="Phocoena_sinus.mPhoSin1.pri.110_CDS",
    description="Vaquita CDS annotations from Morin et al 2020 using NCBI's pipeline",
    url=(
        "https://ftp.ensembl.org/pub/release-110/"
        "gff3/phocoena_sinus/Phocoena_sinus.mPhoSin1.pri.110.chr.gff3.gz"
    ),
    gff_sha256="dbc8640f8f2e046d063c45f160f39edcf729fc39549726603e389cb7b6ce542e",
    intervals_url=(
        "https://stdpopsim.s3-us-west-2.amazonaws.com/"
        "annotations/PhoSin/Phocoena_sinus.mPhoSin1.pri.110_CDS.tar.gz"
    ),
    intervals_sha256="3e30657d2db3ac39136684c56d91b6fd04a03112b139b5dd1a2a1b276f456a44",
    citations=[
        stdpopsim.Citation(
            year=2020,
            author="Morin et al",
            doi="https://doi.org/10.1111/1755-0998.13284",
            reasons={stdpopsim.CiteReason.ANNOTATION},
        )
    ],
    file_pattern="Phocoena_sinus.mPhoSin1.pri.110_CDS_{id}.txt",
    annotation_source="ensembl",
    annotation_type="CDS",
)
_species.add_annotations(_an2)
