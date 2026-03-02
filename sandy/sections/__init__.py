"""
Lazy-loading initializer for `sandy.sections`.

Why this file exists
--------------------
The `sections` package contains many ENDF-6 MF submodules (mf1, mf2, mf3, ...).
In normal Python packages, importing `sandy.sections` does *not* automatically
import submodules such as `sandy.sections.mf1`. Therefore, attribute access like:

    sandy.sections.mf1

would fail unless the submodule was explicitly imported elsewhere.

Why we implement __getattr__
----------------------------
To preserve fast imports and avoid pulling in all MF modules eagerly,
we keep `__init__.py` minimal and use `__getattr__` to load a submodule
*only when it is actually requested*. For example:

    sandy.sections.mf1       -> import sandy.sections.mf1 on demand
    sandy.sections.mf1.read_mf1

This keeps the `sections` namespace lightweight and avoids unnecessary
imports during `import sandy`, while still allowing fully-qualified
access to each MF module.

Why __all__ is defined
----------------------
The `__all__` list declares which MF submodules belong to this package.
It allows tools such as IDEs, autocomplete engines, and static analyzers
to discover valid attributes. It also ensures that __dir__() reports a
complete list of available MF modules.

In summary:
-----------
- No MF modules are imported at package-import time.
- Accessing sandy.sections.mfX triggers lazy loading of that module.
- This avoids circular imports and preserves sandy's fast startup.
- Fully-qualified imports (sandy.sections.mf1) work as expected.
- Top-level API logic in `sandy.__init__` remains unaffected.
"""

from importlib import import_module

__all__ = [
    "mf1", "mf2", "mf3", "mf4", "mf5", "mf6", "mf7", "mf8", "mf9",
    "mf10", "mf31", "mf32", "mf33", "mf34", "mf35", "mf40",
]

def __getattr__(name):
    """
    Lazy attribute resolution. 
    
    To keep imports fast and avoid loading all MF modules eagerly, this module
    implements ``__getattr__`` (PEP 562). When an attribute such as
    ``sandy.sections.mf1`` is accessed, the corresponding module is imported
    *on demand*.

    Examples
    --------

    Below is a minimal smoke test that ensures all declared MF modules can be
    imported through the lazy loader.
    
    >>> import sandy.sections as sec
    >>> for name in sec.__all__:
    ...     obj = getattr(sec, name)
    ...     assert obj is not None, f"Failed to import sections.{name}"
    
    This test does *not* validate semantics. It only that lazy loading works and
    ``__all__`` is complete.

    """
    if name in __all__:
        return import_module(f"sandy.sections.{name}")
    raise AttributeError(f"module 'sandy.sections' has no attribute '{name}'")

def __dir__():
    return sorted(list(globals().keys()) + __all__)
