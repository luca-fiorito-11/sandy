"""
This module contains all classes and functions necessary to work with
the ALEPH-2 fission yield files produced as partial outputs in a burnup
calculation.
"""

__author__ = "Luca Fiorito"
__all__ = [
    "read_fy_output",
    ]

import pandas as pd


def read_fy_output(file):
    """
    Read ALEPH-2 fission yield output file.

    Parameters
    ----------
    file : `str`
        ALEPH-2 fission yield output file.

    Returns
    -------
    df : `pandas.DataFrame`
        dataframe with fission yields.
    """
    df = pd.read_csv(file, sep='\s+', skiprows=1, index_col=0).iloc[:-1]
    df.index.name = "DAUGHTER"
    df.index = df.index.astype(int)
    df.columns.name = "PARENT"
    df.columns = df.columns.astype(int)
    return df
