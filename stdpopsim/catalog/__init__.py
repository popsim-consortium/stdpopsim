import pathlib


# Import all species definitions in the catalog.
__all__ = []
for path in pathlib.Path(__path__[0]).glob("*"):
    module_name = path.parts[-1]
    if module_name[0].isupper():
        __all__.append(module_name)
