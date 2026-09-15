"""
Tests for the package version attribute.
"""

import importlib
import sys

import stdpopsim


class TestVersion:
    def test_version_is_string(self):
        assert isinstance(stdpopsim.__version__, str)
        assert len(stdpopsim.__version__) > 0

    def test_fallback_without_version_file(self):
        # setuptools_scm writes stdpopsim/_version.py when the package is
        # installed. A bare source checkout has no such file, and the package
        # must still import. A None entry in sys.modules makes the import
        # raise ImportError, which mimics the missing file. The attribute on
        # the package must go too, because `from . import _version` returns
        # an existing attribute without consulting sys.modules.
        saved = sys.modules.get("stdpopsim._version")
        sys.modules["stdpopsim._version"] = None
        if hasattr(stdpopsim, "_version"):
            delattr(stdpopsim, "_version")
        try:
            importlib.reload(stdpopsim)
            assert stdpopsim.__version__ == "undefined"
        finally:
            if saved is None:
                del sys.modules["stdpopsim._version"]
            else:
                sys.modules["stdpopsim._version"] = saved
            importlib.reload(stdpopsim)
        assert stdpopsim.__version__ != "undefined"
