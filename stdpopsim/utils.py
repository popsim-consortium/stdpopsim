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
