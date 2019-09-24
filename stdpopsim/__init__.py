# Main entry point for stdpopsim

__version__ = "undefined"
try:
    from . import _version
    __version__ = _version.version
except ImportError:
    pass

# Cache handling for downloaded data.
from . cache import *  # NOQA

# Internal modules. Import here to flatten the namespace.
from . genetic_maps import *  # NOQA
from . models import *  # NOQA
from . genomes import *  # NOQA

# TODO what should we really be exporting here as the public API?

# Species definitions.
from . import homo_sapiens  # NOQA
from . import pongo  # NOQA
from . import arabidopsis_thaliana  # NOQA
from . import e_coli  # NOQA
from . import drosophila_melanogaster  # NOQA
