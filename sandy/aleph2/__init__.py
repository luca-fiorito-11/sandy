"""
Lazy-loading initializer for `sandy.aleph2`.

See `sandy.sections.__init__`
"""

from importlib import import_module

__all__ = [
    "fy_output", "matrix_output", "output_file", "xs_output", "xsfile"
]

def __getattr__(name):
    """
    Lazy attribute resolution. 
    
    To keep imports fast and avoid loading all MF modules eagerly, this module
    implements ``__getattr__`` (PEP 562). When an attribute such as
    ``sandy.al2.output_file`` is accessed, the corresponding module is imported
    *on demand*.

    Examples
    --------

    Below is a minimal smoke test that ensures all declared modules in `__all__` can be
    imported through the lazy loader.
    
    >>> import sandy.aleph2 as al2
    >>> for name in mcnp.__all__:
    ...     obj = getattr(al2, name)
    ...     assert obj is not None, f"Failed to import sections.{name}"
    
    This test does *not* validate semantics. It only that lazy loading works and
    ``__all__`` is complete.

    """
    if name in __all__:
        return import_module(f"sandy.aleph2.{name}")
    raise AttributeError(f"module 'sandy.aleph2' has no attribute '{name}'")

def __dir__():
    return sorted(list(globals().keys()) + __all__)