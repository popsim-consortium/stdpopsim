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
from . engines import * # NOQA

# We import these here to build the catalog, but the internal classes
# defined are not part of the external API.
from .catalog import homo_sapiens  # NOQA
from .catalog import pongo_pygmaeus  # NOQA
from .catalog import arabidopsis_thaliana  # NOQA
from .catalog import e_coli  # NOQA
from .catalog import drosophila_melanogaster  # NOQA
from .catalog import canis_familiaris  # NOQA

from . slim_engine import * # NOQA
