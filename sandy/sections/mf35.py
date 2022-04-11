# -*- coding: utf-8 -*-
"""
This module contains only one public functions:

    * `read_mf34`

Function `read` reads a MF34/MT section from a string and produces a content
object with a dictionary-like structure.
The content object can be accessed using most of the keywords specified in
the ENDF6 manual for this specific MF section.
"""
import sandy

__author__ = "Aitor Bengoechea"


def read_mf35(tape, mat, mt):
    """
    Parse MAT/MF=35/MT section from `sandy.Endf6` object and return
    structured content in nested dcitionaries.

    Parameters
    ----------
    tape : `sandy.Endf6`
        endf6 object containing requested section
    mat : `int`
        MAT number
    mt : `int`
        MT number

    Returns
    -------
    `dict`
        Content of the ENDF-6 tape structured as nested `dict`.

    Examples
    --------
    >>> tape = sandy.get_endf6_file("jeff_33",'xs',922380)
    >>> read_mf35(tape, mat=9237, mt=18)["SUB"][1]["FKK"][0:15]
    [9.63626e-37,
     6.16874e-36,
     8.95968e-36,
     3.04725e-35,
     1.95073e-34,
     2.8333e-34,
     9.63626e-34,
     6.16874e-33,
     8.95968e-33,
     3.0472500000000003e-32,
     1.95073e-31,
     2.8333e-31,
     9.63625e-31,
     6.16871e-30,
     8.959570000000001e-30]
    """
    mf = 35
    df = tape._get_section_df(mat, mf, mt)
    out = {"MAT": mat,
           "MF": mf,
           "MT": mt}
    i = 0
    C, i = sandy.read_cont(df, i)
    out.update({"ZA": C.C1,
                "AWR": C.C2,
                "NK": C.N1,
                "SUB": {}})
    for k in range(out["NK"]):
        L, i = sandy.read_list(df, i)
        D = {"ELO": L.C1,
             "EHI": L.C2,
             "LS": L.L1,
             "LB": L.L2,
             "NE": L.N2,
             "EK": L.B[:L.N2],
             "FKK": L.B[L.N2:]}
        out["SUB"].update({k+1: D})
    return out
