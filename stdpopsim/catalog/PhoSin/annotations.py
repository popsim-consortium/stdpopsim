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
    intervals_sha256="12e767581efc46c5494f800e17011d3e80eb000950f51a98308249afdeb68f87",
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
    intervals_sha256="f44859c4ce99e93f1bacce1da1e620ed41eb2fe1bd75aafe8ac99d54cb6e7a7a",
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
