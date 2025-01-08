import pytest

import stdpopsim
from tests import test_species


class TestSpeciesData(test_species.SpeciesTestBase):

    species = stdpopsim.get_species("CanFam")

    def test_ensembl_id(self):
        assert self.species.ensembl_id == "canis_lupus_familiaris"

    def test_name(self):
        assert self.species.name == "Canis familiaris"

    def test_common_name(self):
        assert self.species.common_name == "Dog"

    def test_qc_population_size(self):
        assert self.species.population_size == 13000

    def test_qc_generation_time(self):
        assert self.species.generation_time == 3


class TestGenomeData(test_species.GenomeTestBase):

    genome = stdpopsim.get_species("CanFam").genome

    # Recombination rates are chromosome-wide averages; so to get this I did:
    # wget \
    # https://github.com/cflerin/dog_recombination/raw/master/dog_genetic_maps.tar.gz
    # tar -xvzf dog_genetic_maps.tar.gz
    # for x in *average*.txt;
    # do y=${x#chr};
    # echo "$y $(tail -n 1 $x) $(head -n 2 $x | tail -n 1)" \
    # | cut -f 2- -d ' ' | sort -k 1 -g;
    # done

    raw_recomb_data = {
        # chr: start, end, cM
        "1": (4283592, 122309715, 90.124765204022),
        "2": (3621442, 85062551, 71.6560923739889),
        "3": (5604604, 91556345, 68.7688581324283),
        "4": (5840941, 87934673, 66.1043548767368),
        "5": (1243143, 88673195, 81.6975593521555),
        "6": (455434, 77489595, 63.1078573213385),
        "7": (180153, 80809723, 58.8068972791936),
        "8": (2763496, 72510424, 57.829359843644),
        "9": (876259, 60812630, 55.6672399393065),
        "10": (2125046, 69293175, 61.1710355665598),
        "11": (4087888, 74253347, 53.6025502008049),
        "12": (82400, 72115946, 55.9056753777465),
        "13": (4067434, 62932928, 49.52362972668),
        "14": (7309849, 60600364, 48.1113329044173),
        "15": (4913124, 64007939, 46.4293482163414),
        "16": (6692748, 58967916, 45.0301626969261),
        "17": (5285642, 63501532, 56.5790593487334),
        "18": (3203856, 55355125, 53.7119190114314),
        "19": (3189264, 53349320, 52.2921177373366),
        "20": (4356904, 58000062, 53.5947241774075),
        "21": (4450666, 50719350, 47.8373474158951),
        "22": (2513263, 61217407, 50.5739055815169),
        "23": (1203392, 52291577, 46.6263217656862),
        "24": (1029209, 47233919, 51.4999687682348),
        "25": (6038820, 51469123, 52.4435907921901),
        "26": (2407850, 38657286, 43.8058911548029),
        "27": (444525, 42191669, 47.0086550820712),
        "28": (4526049, 40963512, 45.4141491786946),
        "29": (913671, 41543120, 44.7477705296246),
        "30": (5433983, 39826282, 40.1964893459467),
        "31": (580882, 39466279, 44.3204605953315),
        "32": (114049, 37999255, 43.7798710211649),
        "33": (447033, 30994965, 40.7491168003185),
        "34": (624977, 41979553, 43.3553617129835),
        "35": (1164664, 26257078, 35.8799002559891),
        "36": (251708, 30523428, 35.9482057940702),
        "37": (1364851, 30583437, 40.4226883720546),
        "38": (246006, 23695770, 33.6825996887619),
    }
    mean_rates = []
    for k in raw_recomb_data:
        mean_rates.append(
            [
                k,
                0.01
                * raw_recomb_data[k][2]
                / (raw_recomb_data[k][1] - raw_recomb_data[k][0]),
            ]
        )
    mean_rates.append(["MT", 0.0])
    # weighted genome-wide average used for X
    total_cM = sum([x[2] for x in raw_recomb_data.values()])
    total_bp = sum([(x[1] - x[0]) for x in raw_recomb_data.values()])
    mean_rates.append(["X", 0.01 * total_cM / total_bp])

    @pytest.mark.parametrize(["name", "rate"], mean_rates)
    def test_recombination_rate(self, name, rate):
        assert rate == pytest.approx(
            self.genome.get_chromosome(name).recombination_rate,
            rel=0.01,
        )

    _genome_mutation_rate = 4e-9

    @pytest.mark.parametrize(
        ["name", "rate"],
        {
            "1": _genome_mutation_rate,
            "2": _genome_mutation_rate,
            "3": _genome_mutation_rate,
            "4": _genome_mutation_rate,
            "5": _genome_mutation_rate,
            "6": _genome_mutation_rate,
            "7": _genome_mutation_rate,
            "8": _genome_mutation_rate,
            "9": _genome_mutation_rate,
            "10": _genome_mutation_rate,
            "11": _genome_mutation_rate,
            "12": _genome_mutation_rate,
            "13": _genome_mutation_rate,
            "14": _genome_mutation_rate,
            "15": _genome_mutation_rate,
            "16": _genome_mutation_rate,
            "17": _genome_mutation_rate,
            "18": _genome_mutation_rate,
            "19": _genome_mutation_rate,
            "20": _genome_mutation_rate,
            "21": _genome_mutation_rate,
            "22": _genome_mutation_rate,
            "23": _genome_mutation_rate,
            "24": _genome_mutation_rate,
            "25": _genome_mutation_rate,
            "26": _genome_mutation_rate,
            "27": _genome_mutation_rate,
            "28": _genome_mutation_rate,
            "29": _genome_mutation_rate,
            "30": _genome_mutation_rate,
            "31": _genome_mutation_rate,
            "32": _genome_mutation_rate,
            "33": _genome_mutation_rate,
            "34": _genome_mutation_rate,
            "35": _genome_mutation_rate,
            "36": _genome_mutation_rate,
            "37": _genome_mutation_rate,
            "38": _genome_mutation_rate,
            "X": _genome_mutation_rate,
            "MT": _genome_mutation_rate,
        }.items(),
    )
    def test_mutation_rate(self, name, rate):
        assert rate == pytest.approx(self.genome.get_chromosome(name).mutation_rate)

    @pytest.mark.parametrize("chrom", [chrom for chrom in genome.chromosomes])
    def test_chromosome_ploidy(self, chrom):
        if chrom.id in ["MT"]:
            assert chrom.ploidy == 1
        else:
            assert chrom.ploidy == 2
