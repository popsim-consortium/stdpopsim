"""
Test python package versioning
"""

from packaging.version import Version

from stdpopsim import _version


class TestPythonVersion:
    """
    Test that the version is PEP440 compliant
    """

    def test_version(self):
        assert str(Version(_version.version)) == _version.version
