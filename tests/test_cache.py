"""
Tests for the cache management code.
"""
import pathlib
import os

import appdirs

import stdpopsim
import tests


class TestSetCacheDir(tests.CacheWritingTest):
    """
    Tests the set_cache_dir function.
    """
    paths = [
        "/somefile", "/some/other/file/", "relative/path", "relative/path/"]

    def test_paths(self):
        for test in self.paths:
            stdpopsim.set_cache_dir(test)
            self.assertEqual(stdpopsim.get_cache_dir(), pathlib.Path(test))
            stdpopsim.set_cache_dir(pathlib.Path(test))
            self.assertEqual(stdpopsim.get_cache_dir(), pathlib.Path(test))

    def test_none(self):
        stdpopsim.set_cache_dir(None)
        cache_dir = pathlib.Path(appdirs.user_cache_dir("stdpopsim", "popgensims"))
        self.assertEqual(stdpopsim.get_cache_dir(), cache_dir)

    def test_environment_var(self):
        try:
            for test in self.paths:
                os.environ["STDPOPSIM_CACHE"] = test
                stdpopsim.set_cache_dir()
                self.assertEqual(stdpopsim.get_cache_dir(), pathlib.Path(test))
        finally:
            os.environ.pop("STDPOPSIM_CACHE")
