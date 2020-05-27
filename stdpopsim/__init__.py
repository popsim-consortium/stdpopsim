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
from . species import *  # NOQA
from . genomes import *  # NOQA
from . cache import *  # NOQA
from . citations import *  # NOQA
from . engines import *  # NOQA

# Add imports for all defined species here.
# We import these here to build the catalog, but the internal functions
# defined are not part of the external API.
from .catalog import AraTha  # NOQA
from .catalog import CanFam  # NOQA
from .catalog import DroMel  # NOQA
from .catalog import EscCol  # NOQA
from .catalog import HomSap  # NOQA
from .catalog import PonAbe  # NOQA

from . import qc  # NOQA

from . slim_engine import *  # NOQA
