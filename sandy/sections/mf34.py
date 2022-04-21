# -*- coding: utf-8 -*-
"""
This module contains a single public function:

    * `read_mf34`

Function `read_mf34` reads a MF34/MT section from a string and produces a 
content object with a dictionary-like structure.
The content object can be accessed using most of the keywords specified in
the ENDF6 manual for this specific MF section.

"""
import sandy

__author__ = "Aitor Bengoechea"

def read_mf34(tape, mat, mt):
    """
    Parse MAT/MF=34/MT section from `sandy.Endf6` object and return
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

    Notes
    -----
    .. note:: The covariance pattern as a function of incident energy values
              (0,1,2,5 and 6 are allowed) are defined as for File 33

    Examples
    --------
    >>> tape = sandy.get_endf6_file("jeff_33", 'xs', 922380)
    >>> read_mf34(tape, mat=9237, mt=2)["REAC"][(0,2)]["P"][(1, 1)]["NI"][0]["FKK"][0:15]
    [1.51558,
     1.51558,
     1.51558,
     0.316347,
     0.376019,
     0.362778,
     0.161342,
     0.0891301,
     0.0634872,
     0.00469711,
     -0.0321085,
     -0.0320439,
     -0.0216885,
     -0.0128297,
     -0.00264976]
    """
    mf = 34
    df = tape._get_section_df(mat, mf, mt)
    out = {"MAT": mat,
           "MF": mf,
           "MT": mt}
    i = 0
    C, i = sandy.read_cont(df, i)
    out.update({"ZA": C.C1,
                "AWR": C.C2,
                "LTT": C.L2,  # Legendre coefficient covariances starting coefficient
                "NMT1": C.N2,  # Number of subsections
                "REAC": {}})
    for j in range(out["NMT1"]):
            C, i = sandy.read_cont(df, i)
            mat1 = C.L1
            mt1 = C.L2
            sub = {"NL": C.N1,  # Number of Legendre coefficients for the MT reaction
                   "NL1": C.N2,  # Number of Legendre coefficients for the MT1 reaction
                   "P": {}}
            nss = C.N1 * (C.N1 + 1) // 2 if C.L2 == out["MT"] else C.N1 * C.N2
            for k in range(nss):
                C, i = sandy.read_cont(df, i)
                l = C.L1  # Index of the Legendre coefficient for reaction MT
                l1 = C.L2  # Index of the Legendre coefficient for reaction MT1
                ni = C.N2
                ssub = {"LCT": C.N1,  # Flag to specify coordinate system
                        "NI": {}}
                for m in range(ni):
                    L, i = sandy.read_list(df, i)
                    sssub = {"LS": L.L1,  # Indicate if the matrix is symmetric
                             "LB": L.L2,  # Flag to indicate the covariance pattern as a function of incident energy.
                             "NT": L.NPL,  # Total number of items
                             "NE": L.N2}  # Energy points
                    if sssub["LB"] in range(5):  # Flag to indicate the covariance pattern as a function of incident energy.
                        if sssub["LS"] == 0:
                            sssub.update({"EK": L.B[::2],
                                          "FK": L.B[1::2]})
                        else:
                            Nk = sssub["NE"] - sssub["LS"]
                            ARRk = L.B[:Nk]
                            ARRl = L.B[Nk:]
                            sssub.update({"EK": ARRk[:Nk/2],
                                          "FK": ARRk[Nk/2:],
                                          "EL": ARRl[:sssub["LS"]],
                                          "FL": ARRl[sssub["LS"]:]})
                    elif sssub["LB"] == 5:
                        sssub.update({"EK": L.B[:L.N2],
                                      "FKK": L.B[L.N2:]})
                    elif sssub["LB"] == 6:
                        sssub.update({"NT": L.NPL,
                                      "NER": L.N2,
                                      "NEC": (L.NPL-1)//L.N2})
                        sssub.update({"EK": L.B[:sssub["NER"]]})
                        sssub.update({"EL": L.B[sssub["NER"]:sssub["NER"]+sssub["NEC"]]})
                        sssub.update({"FKL": L.B[sssub["NER"]+sssub["NEC"]:]})
                    ssub["NI"].update({m: sssub})
                sub["P"].update({(l, l1): ssub})
            out["REAC"].update({(mat1, mt1): sub})
    return out
