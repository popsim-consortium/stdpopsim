"""
Tests for the cache management code.
"""

import os
import pathlib
import tempfile
import tarfile

import appdirs
import pytest

import stdpopsim
from stdpopsim import utils
import tests


class TestSetCacheDir(tests.CacheWritingTest):
    """
    Tests the set_cache_dir function.
    """

    paths = ["/somefile", "/some/other/file/", "relative/path", "relative/path/"]

    def test_paths(self):
        for test in self.paths:
            stdpopsim.set_cache_dir(test)
            assert stdpopsim.get_cache_dir() == pathlib.Path(test)
            stdpopsim.set_cache_dir(pathlib.Path(test))
            assert stdpopsim.get_cache_dir() == pathlib.Path(test)

    def test_none(self):
        stdpopsim.set_cache_dir(None)
        cache_dir = pathlib.Path(appdirs.user_cache_dir("stdpopsim", "popgensims"))
        assert stdpopsim.get_cache_dir() == cache_dir

    def test_environment_var(self):
        try:
            for test in self.paths:
                os.environ["STDPOPSIM_CACHE"] = test
                stdpopsim.set_cache_dir()
                assert stdpopsim.get_cache_dir() == pathlib.Path(test)
        finally:
            os.environ.pop("STDPOPSIM_CACHE")


class TestCachedData(tests.CacheWritingTest):
    def test_caching(self):
        for extract in (True, False):
            with tempfile.TemporaryDirectory() as tmpdir:
                tmpdir = pathlib.Path(tmpdir)
                with utils.cd(tmpdir):
                    filename = "test.foo"
                    with open(filename, "w") as f:
                        print("foo", file=f)
                    tar = tmpdir / "test.tgz"
                    with tarfile.open(tar, "w:gz") as tf:
                        tf.add(filename)

                sha256 = utils.sha256(tar)
                cache = stdpopsim.CachedData(
                    namespace="test",
                    url=tar.resolve().as_uri(),
                    sha256=sha256,
                    extract=extract,
                )
                assert not (cache.is_cached())
                assert not (cache.is_valid())
                cache.download()
                assert cache.is_cached()
                assert cache.is_valid()

                # try to download with incorrect checksum
                cache.sha256 = "1234"
                assert cache.is_cached()
                assert not (cache.is_valid())
                with pytest.raises(ValueError):
                    # checksum mismatch
                    cache.download()
                assert not (cache.is_cached())
                assert not (cache.is_valid())

                # fix the checksum and download again
                cache.sha256 = sha256
                cache.download()
                assert cache.is_cached()
                assert cache.is_valid()

    def test_multiple_threads_downloading(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = pathlib.Path(tmpdir)
            with utils.cd(tmpdir):
                filename = "test.foo"
                with open(filename, "w") as f:
                    print("foo", file=f)
                tar = tmpdir / "test.tgz"
                with tarfile.open(tar, "w:gz") as tf:
                    tf.add(filename)

            cache = stdpopsim.CachedData(
                namespace="test",
                url=tar.resolve().as_uri(),
                sha256=utils.sha256(tar),
                extract=True,
            )
            cache.download()
            # Trick the download code into thinking there's several happening
            # concurrently
            cache.is_cached = lambda: False
            with pytest.warns(UserWarning, match="multiple processes downloading"):
                cache.download()
