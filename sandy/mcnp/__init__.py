"""
Lazy-loading initializer for `sandy.mcnp`.

See `sandy.sections.__init__`
"""

from importlib import import_module

__all__ = [
    "mctal", "mesthal", "output_file",
]

def __getattr__(name):
    if name in __all__:
        return import_module(f"sandy.sections.{name}")
    raise AttributeError(f"module 'sandy.sections' has no attribute '{name}'")

def __dir__():
    return sorted(list(globals().keys()) + __all__)