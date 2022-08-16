import importlib
import pathlib

from . import ensembl_info  # noqa: F401
import stdpopsim

# Import all species definitions in the catalog.
for path in pathlib.Path(__path__[0]).glob("*"):
    module_name = path.parts[-1]
    if module_name[0].isupper():
        importlib.import_module("stdpopsim.catalog." + module_name)

for species in stdpopsim.all_species():
    path = pathlib.Path(__path__[0]) / "demographic-models" / species.id
    for yaml_file in path.glob("*.yaml"):
        dm = stdpopsim.DemographicModel.from_yaml(yaml_file)
        species.add_demographic_model(dm)
