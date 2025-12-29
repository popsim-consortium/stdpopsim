import stdpopsim

from . import genome_data

# Recombination rates for autosomes are obtained from
# Wall et. al. 2022 GBE
# The autosomes recombination rates were computed as simple averages
# across all windows.
# Recombination for X chromosome has not been estimated
# and is assigned the mean value of other chromosomes
_recombination_rate = {
    "1": 9.926379e-09,
    "2": 9.605435e-09,
    "3": 9.022377e-09,
    "4": 9.825128e-09,
    "5": 9.579804e-09,
    "6": 1.049788e-08,
    "7": 1.118884e-08,
    "8": 1.108988e-08,
    "9": 1.132883e-08,
    "10": 1.175322e-08,
    "11": 1.184026e-08,
    "12": 1.082400e-08,
    "13": 1.246772e-08,
    "14": 1.274188e-08,
    "15": 1.260836e-08,
    "16": 1.476158e-08,
    "17": 1.524101e-08,
    "18": 1.368410e-08,
    "19": 1.303735e-08,
    "20": 1.677201e-08,
    "X": 1.18898e-08,
    "Y": 0.0,
}

# Generic and chromosome-specific ploidy
_species_ploidy = 2
_ploidy = {
    "1": _species_ploidy,
    "2": _species_ploidy,
    "3": _species_ploidy,
    "4": _species_ploidy,
    "5": _species_ploidy,
    "6": _species_ploidy,
    "7": _species_ploidy,
    "8": _species_ploidy,
    "9": _species_ploidy,
    "10": _species_ploidy,
    "11": _species_ploidy,
    "12": _species_ploidy,
    "13": _species_ploidy,
    "14": _species_ploidy,
    "15": _species_ploidy,
    "16": _species_ploidy,
    "17": _species_ploidy,
    "18": _species_ploidy,
    "19": _species_ploidy,
    "20": _species_ploidy,
    "X": _species_ploidy,
    "Y": 1,
}

_batra2020 = stdpopsim.Citation(
    author="Batra et. al.",
    year=2020,
    doi="https://doi.org/10.1093/gigascience/giaa134",
    reasons={stdpopsim.CiteReason.ASSEMBLY},
)

_wall2022 = stdpopsim.Citation(
    author="Wall et. al.",
    year=2022,
    doi="https://doi.org/10.1093/gbe/evac040",
    reasons={stdpopsim.CiteReason.REC_RATE},
)

_wu2020 = stdpopsim.Citation(
    author="Wu et. al.",
    year=2020,
    doi="https://doi.org/10.1371/journal.pbio.3000838",
    reasons={stdpopsim.CiteReason.MUT_RATE},
)

# mutation rate from Wu et. al. 2020 PLoS Biology
# recombination rates from Wall et. al. 2022 GBE
_chromosomes = []
for name, data in genome_data.data["chromosomes"].items():
    _chromosomes.append(
        stdpopsim.Chromosome(
            id=name,
            length=data["length"],
            synonyms=data["synonyms"],
            mutation_rate=5.7e-9,
            recombination_rate=_recombination_rate[name],
            ploidy=_ploidy[name],
        )
    )

_genome = stdpopsim.Genome(
    chromosomes=_chromosomes,
    assembly_name=genome_data.data["assembly_name"],
    assembly_accession=genome_data.data["assembly_accession"],
    citations=[
        _batra2020.because(stdpopsim.CiteReason.ASSEMBLY),
        _wall2022.because(stdpopsim.CiteReason.REC_RATE),
        _wu2020.because(stdpopsim.CiteReason.MUT_RATE),
    ],
)

_species = stdpopsim.Species(
    id="PapAnu",
    ensembl_id="papio_anubis",
    name="Papio anubis",
    common_name="Olive baboon",
    separate_sexes=True,
    genome=_genome,
    generation_time=11,  # Generation time from Wu et al section
    # "Inferring split times of humans and baboons"
    population_size=335505,  # Most recent from Wall et al demographic model
    # included in demographic_models.py in this directory
    ploidy=_species_ploidy,
    citations=[
        _wall2022.because(stdpopsim.CiteReason.POP_SIZE),
        _wu2020.because(stdpopsim.CiteReason.GEN_TIME),
    ],
)

stdpopsim.register_species(_species)
