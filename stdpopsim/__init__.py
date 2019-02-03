# Main entry point for stdpopsim

__version__ = "undefined"
try:
    from . import _version
    __version__ = _version.version
except ImportError:
    pass

# Internal modules. Import here to flatten the namespace.
from . genetic_maps import *  # NOQA
from . models import *  # NOQA
from . genomes import *  # NOQA

# Species definitions.
from . import homo_sapiens  # NOQA
from . import pongo  # NOQA
