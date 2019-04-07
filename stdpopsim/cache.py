"""
Cache handling for downloaded data.
"""
import pathlib
import logging
import os

import appdirs

logger = logging.getLogger(__name__)

_cache_dir = None


def set_cache_dir(cache_dir=None):
    """
    The cache_dir is the directory in which stdpopsim stores and checks for
    downloaded data. If the specified cache_dir is not None, this value is
    converted to a pathlib.Path instance, which is used as the cache directory.
    If cache_dir is None (the default), the cache directory is set either from
    the environment variable `STDPOPSIM_CACHE` if it exists, or set to the
    default location using the :mod:`appdirs` module.

    No checks for existance, writability, etc. are performed by this function.
    """
    if cache_dir is None:
        cache_dir = os.environ.get("STDPOPSIM_CACHE", None)
    if cache_dir is None:
        cache_dir = appdirs.user_cache_dir("stdpopsim", "popgensims")
    global _cache_dir
    _cache_dir = pathlib.Path(cache_dir)
    logging.info(f"Set cache_dir to {_cache_dir}")


def get_cache_dir():
    """
    Returns the directory used to cache material downloaded by stdpopsim as a
    pathlib.Path instance. Defaults to a directory 'stdpopsim' in a user cache directory
    (e.g., ~/.cache/stdopsim on Unix flavours). See the :func:`.set_cache_dir` function
    for how this value can be set.
    """
    return _cache_dir


set_cache_dir()
