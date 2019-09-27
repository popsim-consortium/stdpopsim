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
from . cache import *  # NOQA

# We import these here to build the catalog, but the internal classes
# defined are not part of the external API.
from . import homo_sapiens  # NOQA
from . import pongo  # NOQA
from . import arabidopsis_thaliana  # NOQA
from . import e_coli  # NOQA
from . import drosophila_melanogaster  # NOQA
