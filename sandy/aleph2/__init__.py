"""
Lazy-loading initializer for `sandy.aleph2`.

See `sandy.sections.__init__`
"""

from importlib import import_module

__all__ = [
    "fy_output", "matrix_output", "output_file", "xs_output", "xsfile"
]

def __getattr__(name):
    if name in __all__:
        return import_module(f"sandy.aleph2.{name}")
    raise AttributeError(f"module 'sandy.aleph2' has no attribute '{name}'")

def __dir__():
    return sorted(list(globals().keys()) + __all__)