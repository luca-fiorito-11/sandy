import logging
import sys
from importlib import import_module
from types import ModuleType
from typing import Dict, Any

testdir = "tests"


class ShutdownHandler(logging.Handler):
    """
    Trigger exit on errors.
    """

    def emit(self, record):
        logging.shutdown()
        sys.exit(1)


class DuplicateFilter(object):
    """
    Define a filter which keeps track of what was logged, and attach it to
    your logger for the duration of a loop.
    """

    def __init__(self):
        self.msgs = set()

    def filter(self, record):
        rv = record.msg not in self.msgs
        self.msgs.add(record.msg)
        return rv


class Error(Exception):
    pass


class ConditionalFormatter(logging.Formatter):
    """Change format dynamically based on log level."""
    
    FORMATS = {
        logging.INFO: logging.Formatter('%(message)s'),  # No level prefix for INFO
        'default': logging.Formatter('%(levelname)s:  %(message)s'),  # Default format
    }
    
    def format(self, record):
        formatter = self.FORMATS.get(record.levelno, self.FORMATS['default'])
        return formatter.format(record)

# Setup logging
logger = logging.getLogger()
handler = logging.StreamHandler()
handler.setFormatter(ConditionalFormatter())  # Apply the custom formatter
logger.addHandler(handler)
logger.setLevel(logging.INFO)


__version__ = '1.1.0'



# ================================================================
# LAZY PUBLIC API
# ================================================================

# List of available sandy submodules (lazy-loaded)
_SUBMODULES  = [
    # modules
    "constants",
    "cov",
    "decay",
    "edistr",
    "endf6",
    "energy_grids",
    "errorr",
    "fy",
    "gendf",
    "gls",
    "libraries",
    "lpc",
    "njoy",
    "pert",
    "records",
    "samples",
    "sampling",
    "settings",
    "spectra",
    "tsl",
    "utils",
    "xs",
    "zam",

    # folders
    "aleph2",
    "mcnp",
    "sections", 
]



"""
_PUBLIC: Top-level lazy public API mapping for the `sandy` package.

Purpose
-------
The `_PUBLIC` dictionary defines the official, stable, user-facing API of the
`sandy` package.  Keys are the attribute names that users can access via:

    import sandy
    sandy.Name

Values are *dotted module paths* (relative to the `sandy` package) that tell
the top-level `__getattr__` resolver where to import each name from.

For example:

    "Xs": "xs"
    "read_mf1": "sections.mf1"

means that:

    sandy.Xs         -> loads sandy/xs.py     and returns Xs
    sandy.read_mf1   -> loads sandy/sections/mf1.py and returns read_mf1

Why this exists
---------------
Python packages do *not* automatically re-export objects from submodules, and
importing all sandy submodules eagerly would make `import sandy` extremely slow
because some submodules import heavy dependencies (numpy, pandas, scipy).

Instead, we provide a **flat, friendly public API** while keeping imports lazy.
This means:

- Users get a clean high-level interface (`sandy.Xs`, `sandy.Endf6`, ...).
- Heavy modules are imported *only when the user accesses them*.
- Subpackages (`sandy.sections`, `sandy.mcnp`, etc.) remain lightweight.
- No circular imports occur at startup.

How `_PUBLIC` interacts with `__getattr__`
------------------------------------------
`sandy.__getattr__` checks:

1. If the attribute is a submodule → lazy-load it.
2. If the attribute is in `_PUBLIC` → import the module named in `_PUBLIC[name]`
   and return the requested attribute.

This means that **calling `sandy.Name` does not import unrelated modules**, only
the exact module declared in `_PUBLIC`.

Naming guidelines
-----------------
- `_PUBLIC` should include only *classes* or *important functions* intended for
  user consumption.
- Do not expose internal helpers, filenames, or modules (e.g. no `mf1`, no
  `records`, etc.).
- Subpackages may re-export their own symbols (e.g. `sandy.mcnp.Mctal`) via
  their own `__init__.py`, but this table governs *top-level* access only.

Testing the public API
----------------------

An example test is in `__getattr__`

Summary
-------
`_PUBLIC` is the single authoritative definition of what the user-facing API
of `sandy` looks like. It provides:

- a stable public surface,
- lazy loading of heavy modules,
- predictable import behavior,
- reduced startup time,
- clean separation between internal structure and external API.

Maintain this table carefully: the top-level behavior of `sandy` depends on it.
"""
_PUBLIC: Dict[str, str] = {
    # ENDF-6 format sections
    "read_mf1": "sections.mf1",
    "write_mf1": "sections.mf1",

    "read_mf2": "sections.mf2",
    "write_mf2": "sections.mf2",

    "read_mf3": "sections.mf3",
    "write_mf3": "sections.mf3",

    "read_mf4": "sections.mf4",
    "write_mf4": "sections.mf4",

    "read_mf5": "sections.mf5",
    "write_mf5": "sections.mf5",

    "read_mf6": "sections.mf6",
    "write_mf6": "sections.mf6",

    "read_mf7": "sections.mf7",
    "write_mf7": "sections.mf7",

    "read_mf8": "sections.mf8",
    "write_mf8": "sections.mf8",

    "read_mf9": "sections.mf9",
    "write_mf9": "sections.mf9",

    "read_mf10": "sections.mf10",
    "write_mf10": "sections.mf10",

    "read_mf31": "sections.mf31",
    "read_mf32": "sections.mf32",
    "read_mf33": "sections.mf33",
    "read_mf34": "sections.mf34",
    "read_mf35": "sections.mf35",
    "read_mf40": "sections.mf40",

    # Covariances
    "CategoryCov": "cov",
    "corr2cov": "cov",

    # Radioactive Decay Data
    "DecayData": "decay",
    "decay_modes": "decay",
    "BranchingRatio": "decay",
    "HalfLife": "decay",
    "DecayEnergy": "decay",

    # Energy distributions
    "Edistr": "edistr",

    # ENDF6
    "Endf6": "endf6",
    "get_endf6_file": "endf6",
    
    # ERRORR
    "Errorr": "errorr",
    
    # Fission yields
    "Fy": "fy",
    "fy_cea_pu239th": "fy",
    "fy_cea_pu239th_corr": "fy",
    "fy_cea_u235th": "fy",
    "fy_cea_u235th_corr": "fy",
    "get_cea_fy": "fy",

    # GENDF
    "Gendf": "gendf",
    
    # Legendre polynomial coefficients
    "Lpc": "lpc",
    
    # NJOY
    "get_njoy": "njoy",

    # Perturbations
    "Pert": "pert",

    # Sampling
    "Samples": "samples",

    # Thermal scattering laws
    "Tsl": "tsl",

    # Cross-sections
    "Xs": "xs",

    # ALEPH2
    "read_fy_output": "aleph2.fy_output",
    "read_matrix_output": "aleph2.matrix_output",
    "OutputFile": "aleph2.output_file",
    "read_output": "aleph2.output_file",
    "read_xs_output": "aleph2.xs_output",
    "AlephFile": "aleph2.xsfile",

    # MCNP
    "MctalTally": "mcnp.mctal",
    "MshtTally": "mcnp.meshtal",

}


# internal cache for loaded modules to avoid repeated imports
_loaded_modules: Dict[str, ModuleType] = {}


def __getattr__(name: str) -> Any:
    """
    Lazy attribute resolution.

    1. If `name` is a submodule, import and return it.
    2. If `name` is a re-exported symbol, import its module and return it.
    
    Examples
    --------

    Below is a minimal smoke test that attempts to import every public symbol from
    the top-level namespace. It is intended to help contributors verify that:
    
    - All `_PUBLIC` entries point to valid modules.
    - The attribute exists in the indicated module.
    - The lazy resolver functions correctly.

    >>> import sandy
    >>> def test_public_api_smoke():
    ...     for name in sandy._PUBLIC:
    ...         obj = getattr(sandy, name)
    ...         assert obj is not None, f"Failed to import sandy.{name}"
    >>> test_public_api_smoke()

    This test does **not** check semantics; it only ensures that every public
    symbol resolves cleanly through the lazy loader.

    Important: This test imports *all* public symbols, so it will pull in heavy
    dependencies during the test phase—this is expected.

    """
    # Case 1: sandy.xs, sandy.cov, sandy.sections, ...
    if name in _SUBMODULES:
        full = f"sandy.{name}"               # FULL dotted path
        mod = _loaded_modules.get(full)
        if mod is None:
            mod = import_module(full)
            _loaded_modules[full] = mod
        return mod

    # Case 2: re-export public symbols
    if name in _PUBLIC:
        full = f"sandy.{_PUBLIC[name]}"      # FULL dotted path
        mod = _loaded_modules.get(full)
        if mod is None:
            mod = import_module(full)
            _loaded_modules[full] = mod
        return getattr(mod, name)

    raise AttributeError(f"module 'sandy' has no attribute '{name}'")



def __dir__():
    """Make autocomplete work nicely."""
    base = set(globals().keys())
    base.update(_SUBMODULES)
    base.update(_PUBLIC.keys())
    return sorted(base)
