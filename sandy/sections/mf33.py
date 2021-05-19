# -*- coding: utf-8 -*-
"""
"""

import sandy

__author__ = "Luca Fiorito"
__all__ = [
        "read_mf33",
        ]


def read_mf33(tape, mat, mt, mf=33):
    """
    Parse MAT/MF=8/MT section for MF31/33 covariance matrix from
    `sandy.Endf6` object and return structured content in nested dcitionaries.

    Parameters
    ----------
    tape : `sandy.Endf6`
        endf6 object containing requested section
    mat : `int`
        MAT number
    mtt : `int`
        MT number
    mf : `int`, optional, default is `33`
        MF number, it also allows parsing nubar covariances with MF=31

    Returns
    -------
    `dict`
        Content of the ENDF-6 tape structured as nested `dict`.

    Notes
    -----
    .. note:: this function can parse MF=31 sections if `mf=31` is passed
              as argument
    """
    df = tape._get_section_df(mat, mf, mt)
    out = {
            "MAT": mat,
            "MF": mf,
            "MT": mt,
            }
    i = 0
    C, i = sandy.read_cont(df, i)
    out["ZA"] = C.C1
    out["AWR"] = C.C2
    out["MTL"] = C.L2
    nsub = C.N2
    subs = {}
    for j in range(nsub):
        C, i = sandy.read_cont(df, i)
        xmf1 = C.C1
        xlfs1 = C.C2
        mat1 = C.L1
        mt1 = C.L2
        sub = {"XMF1": xmf1, "XLFS1": xlfs1, "MAT1": mat1, "MT1": mt1}
        nc = C.N1  # number of NC-type sections
        ni = C.N2  # number of NI-type sections
        ncdict = {}
        for k in range(nc):
            C, i = sandy.read_cont(df, i)
            lty = C.L2
            subsub = {"LTY": lty}
            L, i = sandy.read_list(df, i)
            if lty == 0:
                subsub["E1"] = L.C1,
                subsub["E2"] = L.C2,
                subsub["CI"] = L.B[:L.N2],
                subsub["XMTI"] = L.B[L.N2:]
            elif lty in (1, 2, 3):
                subsub["E1"] = L.C1
                subsub["E2"] = L.C2
                subsub["MATS"] = L.L1
                subsub["MTS"] = L.L2
                subsub["XMFS"] = L.B[0]
                subsub["XLFSS"] = L.B[1]
                subsub["EI"] = L.B[2:2 + L.N2]
                subsub["WEI"] = L.B[2 + L.N2:]
            else:
                raise AssertionError
            ncdict[k] = subsub
        sub["NC"] = ncdict
        nidict = {}
        for k in range(ni):
            L, i = sandy.read_list(df, i)
            lb = L.L2
            subsub = {"LB": lb}
            if lb in [0, 1, 2, 3, 4]:
                subsub["LT"] = L.L1
                subsub["NT"] = L.NPL
                subsub["NP"] = L.N2
                if subsub["LT"] == 0:
                    subsub["EK"] = L.B[::2]
                    subsub["FK"] = L.B[1::2]
                else:  # found in 26-Fe-54g.jeff33, MF33/MT103, LB=4 LT=1
                    nk = subsub["NP"] - subsub["LT"]
                    arrk = L.B[:nk*2]
                    subsub["EK"] = arrk[::2]
                    subsub["FK"] = arrk[1::2]
                    arrl = L.B[nk*2:]
                    subsub["EL"] = arrl[::2]
                    subsub["FL"] = arrl[1::2]
            elif subsub["LB"] == 5:
                subsub["LS"] = L.L1
                subsub["NT"] = L.NPL
                subsub["NE"] = ne = L.N2
                subsub["EK"] = L.B[:ne]
                subsub["FKK"] = L.B[ne:]
            elif subsub["LB"] == 6:
                subsub["NT"] = L.NPL
                subsub["NER"] = ner = L.N2
                subsub["NEC"] = nec = (L.NPL - 1) // L.N2
                subsub["EK"] = L.B[:ner]
                subsub["EL"] = L.B[ner:ner + nec]
                subsub["FKL"] = L.B[ner + nec:]
            elif subsub["LB"] in (8, 9):
                subsub["LT"] = L.L1
                subsub["NT"] = L.NPL
                subsub["NP"] = np = L.N2
                subsub["EK"] = L.B[:np]
                subsub["FK"] = L.B[np:]
            else:
                raise AssertionError("LB not found")
            nidict[k] = subsub
        sub["NI"] = nidict
        subs[mat1*1000 + mt1] = sub
    if subs:
        out["SUB"] = subs
    return out
