"""
Collection of utilities, functions and classes that are requetsed in all code
components.
"""
from itertools import zip_longest

from sandy.shared import pad_from_beginning, \
                         pad_from_beginning_fast, \
                         uniform_loggrid

__author__ = "Luca Fiorito"
__all__ = [
    "grouper",
    "pad_from_beginning",
    "pad_from_beginning_fast",
    "uniform_loggrid",
    ]


def grouper(iterable, n, fillvalue=None):
    """
    Collect data into fixed-length chunks or blocks
    """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)
