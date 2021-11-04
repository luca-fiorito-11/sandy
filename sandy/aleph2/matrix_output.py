"""
This module contains all classes and functions necessary to work with
the ALEPH-2 reaction rate matrix files produced as partial outputs in a burnup
calculation.
"""

__author__ = "Luca Fiorito"
__all__ = [
    "read_matrix_output",
    ]

import pandas as pd


def read_matrix_output(file):
    """
    Read ALEPH-2 reaction rate matrix output file.
    Entries are saved without an index and with the following colums:
        - J: column index in the transition matrix
        -"NUCLIDE_J: nuclide identifier of the daughter nuclide
        - DAUGHTER: ZAM of the daughter nuclide
        - I: row index in the transition matrix
        - NUCLIDE_I: nuclide identifier of the parent nuclide
        - PARENT: ZAM of the parent nuclide
        - CONC: concentration of the parent nuclide in $at/cm^3$
        - RR: reaction rate in $1/s$
        - JAC: jacobian in $at/cm^3/s$

    Parameters
    ----------
    file : `str`
        ALEPH-2 reaction rate matrix output file.

    Returns
    -------
    df : `pandas.DataFrame`
        dataframe with reaction rates.
    """
    widths = [10, 12, 10, 10, 14, 11, 15, 16, 19]
    names = [
        "J",
        "NUCLIDE_J",
        "DAUGHTER",
        "I",
        "NUCLIDE_I",
        "PARENT",
        "CONC",
        "RR",
        "JAC",
        ]
    df = pd.read_fwf(
        file,
        widths=widths,
        skiprows=2,
        names=names,
        )
    return df
