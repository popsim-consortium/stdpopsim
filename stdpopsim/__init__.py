# Main entry point for stdpopsim

__version__ = "undefined"
try:
    from . import _version
    __version__ = _version.version
except ImportError:
    pass

from . genetic_maps import *  # NOQA
from . models import *  # NOQA
from . genomes import *  # NOQA
from . import h_sapiens  # NOQA
from . import pongo  # NOQA
