"""
Package definition for tests. Defined to allow cross-importing.
"""
import tempfile

import stdpopsim


class CacheReadingTest:
    """
    This should be used as the superclass of all tests that use downloaded data.
    Rather than using the standard userdir for downloaded data, we use this
    local directory that can be easily removed and modified without fear of
    interfering with production code.
    """

    cache_dir = "test_cache"
    saved_cache_dir = None
    saved_urls = {}

    @classmethod
    def setup_class(cls):
        cls.saved_cache_dir = stdpopsim.get_cache_dir()
        stdpopsim.set_cache_dir(cls.cache_dir)

    @classmethod
    def teardown_class(cls):
        stdpopsim.set_cache_dir(cls.saved_cache_dir)


class CacheWritingTest:
    """
    This should be used as the superclass of all tests that alter the
    downloaded data cache in any non-standard way.
    """

    saved_cache_dir = None
    tmp_cache_dir = None

    @classmethod
    def setup_class(cls):
        cls.saved_cache_dir = stdpopsim.get_cache_dir()
        cls.tmp_cache_dir = tempfile.TemporaryDirectory()
        stdpopsim.set_cache_dir(cls.tmp_cache_dir.name)

    @classmethod
    def teardown_class(cls):
        stdpopsim.set_cache_dir(cls.saved_cache_dir)
        del cls.tmp_cache_dir
