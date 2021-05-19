"""
This module contains all classes and functions necessary to work with
the ALEPH-2 cross section files produced as partial outputs in a burnup
calculation.
"""

__author__ = "Luca Fiorito"
__all__ = [
    "read_xs_output",
    ]

import pandas as pd


def read_xs_output(file):
    """
    Read ALEPH-2 cross section output file.

    Parameters
    ----------
    file : `str`
        ALEPH-2 cross section output file.

    Returns
    -------
    df : `pandas.DataFrame`
        dataframe with cross sections.
    """
    df = pd.read_csv(file,
                     sep='\s{2,}',
                     engine="python",
                     skiprows=1,
                     index_col=0,
                     )
    df.index.name = "NUCLIDE"
    df.index = df.index.astype(int)
    df.columns.name = "REACTION"
    return df
