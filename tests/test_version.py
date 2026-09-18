"""
Tests for the package version.
"""

from packaging.version import Version

import stdpopsim


class TestVersion:
    def test_version_is_pep440(self):
        # A version that parses and round-trips unchanged is PEP 440 compliant.
        assert str(Version(stdpopsim.__version__)) == stdpopsim.__version__
