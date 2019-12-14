"""
Tests for the utils module.
"""
import unittest


from stdpopsim import utils


class TestValidDemographicModelId(unittest.TestCase):
    """
    Tests for the is_valid_demographic_model_id function.
    """

    def test_empty_string(self):
        self.assertFalse(utils.is_valid_demographic_model_id(""))

    def test_contains_spaces(self):
        bad_ids = [
            "CamelCase 4X19", " CamelCase4X19", "CamelCase4X19"]
        for bad_id in bad_ids:
            self.assertFalse(utils.is_valid_demographic_model_id(bad_id))

    def test_too_many_sections(self):
        bad_ids = [
            "_CamelCase_4X19", "CamelCase_4X19_", "CamelCase__4X19",
            "CamelCase_4X19_1234_5678"]
        for bad_id in bad_ids:
            self.assertFalse(utils.is_valid_demographic_model_id(bad_id))

    def test_bad_name(self):
        bad_ids = [
            "camelCase_4X19", "1CamelCase_4X19", "Camel-Case_4X19"]
        for bad_id in bad_ids:
            self.assertFalse(utils.is_valid_demographic_model_id(bad_id))

    def test_bad_populations(self):
        bad_ids = [
            "CamelCase_0X19", "CamelCase_0.1X19", "CamelCase_OneX19"]
        for bad_id in bad_ids:
            self.assertFalse(utils.is_valid_demographic_model_id(bad_id))

    def test_bad_author(self):
        bad_ids = [
            "CamelCase_1019", "CamelCase_1-19", "CamelCase_1;19"]
        for bad_id in bad_ids:
            self.assertFalse(utils.is_valid_demographic_model_id(bad_id))

    def test_bad_year(self):
        bad_ids = [
            "CamelCase_10X1", "CamelCase_1X", "CamelCase_1X2001",
            "CamelCase_1Xfive"]
        for bad_id in bad_ids:
            self.assertFalse(utils.is_valid_demographic_model_id(bad_id))

    def test_good_ids(self):
        good_ids = [
            "CamelCase_10X10", "CamelCase_1Y00", "CamelCase_100A05",
            "Camel_1X10", "C_1Y00", "C01234_1Y00"]
        for good_id in good_ids:
            self.assertTrue(utils.is_valid_demographic_model_id(good_id))


class TestValidSpeciesId(unittest.TestCase):
    """
    Tests for the is_valid_demographic_model_id function.
    """

    def test_empty_string(self):
        self.assertFalse(utils.is_valid_species_id(""))

    def test_contains_spaces(self):
        bad_ids = ["Cam Cas", " CamCas", "CamCas "]
        for bad_id in bad_ids:
            self.assertFalse(utils.is_valid_species_id(bad_id))

    def test_bad_ids(self):
        bad_ids = [
            "camCas", "camcas", "CAMCAS", "C1mCas", "Camcas", "CamelCase"]
        for bad_id in bad_ids:
            self.assertFalse(utils.is_valid_species_id(bad_id))

    def test_good_ids(self):
        bad_ids = ["CamCas", "HomSap"]
        for bad_id in bad_ids:
            self.assertTrue(utils.is_valid_species_id(bad_id))


class TestValidSpeciesName(unittest.TestCase):
    """
    Tests for the is_valid_demographic_model_id function.
    """

    def test_empty_string(self):
        self.assertFalse(utils.is_valid_species_name(""))

    def test_extra_spaces(self):
        bad_names = ["Pan  pan", " Pan pan", "Pan pan "]
        for bad_name in bad_names:
            self.assertFalse(utils.is_valid_species_name(bad_name))

    def test_bad_names(self):
        bad_names = [
            "PAN PAN", "Pan Pan", "Panpan", "PanPan" "pan pan", "Pan p0n"
            "Pan, pan"]
        for bad_name in bad_names:
            self.assertFalse(utils.is_valid_species_name(bad_name))

    def test_good_names(self):
        good_names = ["Pan pan", "Homo sapiens"]
        for good_name in good_names:
            self.assertTrue(utils.is_valid_species_name(good_name))

    @unittest.skip("Implement more flexible species name")
    def test_three_or_more(self):
        bad_names = ["Pan pan pan"]
        for bad_name in bad_names:
            self.assertTrue(utils.is_valid_species_name(bad_name))


class TestValidSpeciesCommonName(unittest.TestCase):
    """
    Tests for the is_valid_demographic_model_id function.
    """

    def test_empty_string(self):
        self.assertFalse(utils.is_valid_species_common_name(""))

    def test_bad_common_names(self):
        # TODO add more examples when we define restrictions (#330)
        bad_common_names = ["0", "human"]
        for bad_common_name in bad_common_names:
            self.assertFalse(utils.is_valid_species_common_name(bad_common_name))

    def test_good_common_names(self):
        # TODO add more examples when we define restrictions (#330)
        good_names = ["Human", "Stuff and things"]
        for good_name in good_names:
            self.assertTrue(utils.is_valid_species_common_name(good_name))
