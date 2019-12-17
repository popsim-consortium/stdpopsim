"""
Miscellaneous utilities.
"""
import re


def is_valid_demographic_model_id(model_id):
    """
    Returns True if the specified string is a valid demographic model ID. This must
    be a string with the following pattern:

    {CamelCaseName}_{num populations}{First letter author name}{2 digit year}.
    """
    regex = re.compile(r"[A-Z][A-Za-z0-9]*_[1-9]\d*[A-Z]\d\d")
    return regex.fullmatch(model_id) is not None


def is_valid_species_id(species_id):
    """
    Returns True if the specified string is a valid species ID. This must
    be a 6 letter CamelCase identifier.
    """
    regex = re.compile(r"[A-Z][a-z]{,2}[A-Z][a-z]{,2}")
    return regex.fullmatch(species_id) is not None


def is_valid_species_name(name):
    """
    Returns True if the specified string is a valid species name. This
    must be two or more words with first letter capitalised.
    """
    # FIXME this only supports two words for now. See #329
    regex = re.compile(r"[A-Z][a-z]+ [a-z]+")
    return regex.fullmatch(name) is not None


def is_valid_species_common_name(common_name):
    """
    Returns True if the specified string is a valid species common name. This
    must start with a capital letter.
    """
    # FIXME any sensible restrictions we can make on common names? See #330.
    regex = re.compile(r"[A-Z].*")
    return regex.fullmatch(common_name) is not None
