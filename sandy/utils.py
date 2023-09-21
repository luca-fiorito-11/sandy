"""
Collection of utilities, functions and classes that are requetsed in all code
components.
"""
from itertools import zip_longest
import numpy as np

from sandy.shared import pad_from_beginning, \
                         pad_from_beginning_fast, \
                         uniform_loggrid

__author__ = "Luca Fiorito"
__all__ = [
    "grouper",
    "pad_from_beginning",
    "pad_from_beginning_fast",
    "uniform_loggrid",
    "get_seed",
    ]


def grouper(iterable, n, fillvalue=None):
    """
    Collect data into fixed-length chunks or blocks
    """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

def get_seed():
    upper_limit = 2**32  # limit from numpy.random.seed
    return int(np.random.rand() * (upper_limit - 1))
    