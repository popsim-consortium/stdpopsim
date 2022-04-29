# Main entry point for stdpopsim

__version__ = "undefined"
try:
    from . import _version

    __version__ = _version.version
except ImportError:
    pass

import pathlib as _pathlib

# Internal modules. Import here to flatten the namespace.
from .genetic_maps import *  # NOQA
from .models import *  # NOQA
from .species import *  # NOQA
from .genomes import *  # NOQA
from .annotations import *  # NOQA
from .dfe import *  # NOQA
from .cache import *  # NOQA
from .citations import *  # NOQA
from .engines import *  # NOQA
from .warning_categories import *  # NOQA

# Extensions.
from . import ext  # NOQA

# We import catalog here, but the internal functions
# defined are not part of the external API.
from .catalog import *  # NOQA

for sp in all_species():
    path = _pathlib.Path(__path__[0]) / "catalog" / "demographic-models" / sp.id
    for yaml_file in path.glob("*.yaml"):
        dm = DemographicModel.from_yaml(yaml_file)
        sp.add_demographic_model(dm)

from . import qc  # NOQA

from .slim_engine import *  # NOQA
