"""
Lazy-loading initializer for `sandy.mcnp`.

See `sandy.sections.__init__`
"""

from importlib import import_module

__all__ = [
    "mctal", "meshtal", "output_file",
]

def __getattr__(name):
    """
    Lazy attribute resolution. 
    
    To keep imports fast and avoid loading all MF modules eagerly, this module
    implements ``__getattr__`` (PEP 562). When an attribute such as
    ``sandy.mcnp.mctal`` is accessed, the corresponding module is imported
    *on demand*.

    Examples
    --------

    Below is a minimal smoke test that ensures all declared modules in `__all__` can be
    imported through the lazy loader.
    
    >>> import sandy.mcnp as mcnp
    >>> for name in mcnp.__all__:
    ...     obj = getattr(mcnp, name)
    ...     assert obj is not None, f"Failed to import sections.{name}"
    
    This test does *not* validate semantics. It only that lazy loading works and
    ``__all__`` is complete.

    """
    if name in __all__:
        return import_module(f"sandy.mcnp.{name}")
    raise AttributeError(f"module 'sandy.mcnp' has no attribute '{name}'")

def __dir__():
    return sorted(list(globals().keys()) + __all__)