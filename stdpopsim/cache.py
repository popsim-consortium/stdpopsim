"""
Cache handling for downloaded data.
"""

import pathlib
import logging
import os
import urllib.parse
import tempfile
import warnings

import appdirs
import attr

from . import utils

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
    logger.info(f"Set cache_dir to {_cache_dir}")


def get_cache_dir():
    """
    Returns the directory used to cache material downloaded by stdpopsim as a
    pathlib.Path instance. Defaults to a directory 'stdpopsim' in a user cache directory
    (e.g., ~/.cache/stdopsim on Unix flavours). See the :func:`.set_cache_dir` function
    for how this value can be set.
    """
    return _cache_dir


set_cache_dir()


@attr.s(kw_only=True)
class CachedData:
    """
    Downloadable data that will be cached locally.

    The downloadable should be a single file. The local cache may
    correspond to this same file, or to its extracted contents. In the latter
    case, the downloaded file will be removed after archive extraction.
    The downloaded file is compared against the expected SHA256 checksum,
    and if correct, the checksum is then also stored locally.

    :ivar str namespace: The namespace under which the cache will be stored.
        This will be converted into a folder, by constructing folders in the
        cache corresponding to each component of the namespace.
        E.g. if we're on a unix system with cache under ``/path/to/cache``, and
        ``namespace="foo/bar"``, the cached data will live under
        ``/path/to/cache/foo/bar``.
    :ivar str url: The URL of the data to be cached.
    :ivar str sha256: The SHA256 checksum of the downloaded file.
    :ivar bool extract: True if the downloaded file is a tarball that should be
        extracted into the cached namespace, False otherwise.
    """

    namespace = attr.ib(type=str)
    url = attr.ib(type=str)
    sha256 = attr.ib(type=str)
    extract = attr.ib(type=bool)

    def __attrs_post_init__(self):
        u = urllib.parse.urlparse(self.url)
        self._basename = pathlib.PurePath(u.path).name

    @property
    def sha256_file(self):
        return get_cache_dir() / self.namespace / f"{self._basename}.sha256"

    @property
    def cache_path(self):
        # the cache path could be a folder or a file, depending on self.extract
        path = get_cache_dir() / self.namespace
        if not self.extract:
            path = path / self._basename
        return path

    def is_cached(self):
        """
        Returns True if the data is cached locally.
        """
        return self.cache_path.exists()

    def is_valid(self):
        """
        Returns True if the cached data matches the checksum.
        """
        is_valid = False
        if self.is_cached() and self.sha256_file.exists():
            with open(self.sha256_file, "r") as f:
                cached_sha256 = f.read().strip()
            is_valid = self.sha256 == cached_sha256
        return is_valid

    def download(self):
        """
        Downloads the file from the source URL and stores it in the cache.
        If the local cache already exists, it is first removed.
        """
        if self.is_cached():
            logger.info(f"Clearing cache {self.cache_path}")
            with tempfile.TemporaryDirectory(dir=get_cache_dir()) as tempdir:
                # Atomically move to a temporary directory, which will be automatically
                # deleted on exit.
                dest = pathlib.Path(tempdir) / "will_be_deleted"
                os.rename(self.cache_path, dest)

        self.cache_path.parent.mkdir(parents=True, exist_ok=True)

        logger.info(f"Downloading {self.url}")
        # os.rename will not work on some Unixes if the source and dest are on
        # different file systems. Keep the tempdir in the same directory as
        # the destination to ensure it's on the same file system.
        with tempfile.TemporaryDirectory(dir=get_cache_dir()) as tempdir:
            tempdir = pathlib.Path(tempdir)
            local_path = tempdir / "downloaded"
            utils.download(self.url, local_path)

            logger.debug("Checking SHA256")
            download_sha256 = utils.sha256(local_path)
            if download_sha256 != self.sha256:
                # TODO: use a more appropriate exception here.
                raise ValueError(
                    f"Expected SHA256={self.sha256}, but downloaded file has"
                    f"{download_sha256}."
                )

            if self.extract:
                extract_dir = tempdir / "extracted"
                extract_dir.mkdir()
                logger.debug(f"Extracting {local_path}")
                utils.untar(local_path, extract_dir)
                local_path = extract_dir

            # If this has all gone OK up to here we can now move the
            # data into the cache location. This should minimise the
            # chances of having malformed data in the cache.
            logger.info(f"Saving to {self.cache_path}")
            # os.rename is atomic, and will raise an OSError if the destination
            # is a directory and already exists. Therefore, if we see the map
            # exists we assume that some other process has already downloaded
            # it, and raise a warning.
            # If the source and destination are regular files (such as when
            # self.extract==False), the destination will be silently replaced
            # on unix systems, but FileExistsError will be raised on windows.
            try:
                os.rename(local_path, self.cache_path)
            except (OSError, FileExistsError):
                warnings.warn(
                    "Error occured renaming map directory. Are multiple processes "
                    "downloading this map at the same time?"
                )
                return

            # Write out the checksum.
            with open(self.sha256_file, "w") as f:
                print(self.sha256, file=f)
