import pdb

import sandy

__author__ = "Luca Fiorito"
__all__ = [
        "read_mf3",
        "write_mf3"
        ]


def read_mf3(tape, mat, mt):
    """
    Parse MAT/MF=3/MT section from `sandy.Endf6` object and return structured
    content in nested dcitionaries.

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
    """
    mf = 3
    df = tape._get_section_df(mat, mf, mt)
    out = {
            "MAT": mat,
            "MF": mf,
            "MT": mt,
            }
    i = 0
    C, i = sandy.read_cont(df, i)
    add = {
            "ZA": C.C1,
            "AWR": C.C2,
            "PFLAG": C.L2,
            }
    out.update(add)
    T, i = sandy.read_tab1(df, i)
    add = {
            "QM": T.C1,
            "QI": T.C2,
            "LR": T.L2,
            "NBT": T.NBT,
            "INT": T.INT,
            "E": T.x,
            "XS": T.y,
            }
    out.update(add)
    return out


def write_mf3(sec):
    """
    Given the content of a MF3 section as nested dictionaries, write it
    to string.

    Returns
    -------
    `str`
        Multiline string reproducing the content of a ENDF-6 section.

    Notes
    -----
    .. note:: The end-of-line records MAT, MF, MT and line number are added at
              the end of each line.

    .. important:: The string does not endf with a newline symbol `\n`.
    """
    lines = sandy.write_cont(
            sec["ZA"],
            sec["AWR"],
            0,
            sec["PFLAG"],
            0,
            0,
            )
    lines += sandy.write_tab1(
            sec["QM"],
            sec["QI"],
            0,
            sec["LR"],
            sec["NBT"],
            sec["INT"],
            sec["E"],
            sec["XS"],
            )
    return "\n".join(sandy.write_eol(lines, sec["MAT"], 3, sec["MT"]))


def _read_errorr(text):
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    out = {"MAT" : MAT, "MF" : MF, "MT" : MT}
    i = 0
    L, i = read_list(str_list, i)
    out.update({"XS" : L.B})
    return out


def _read_groupr(text):
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    out = {"MAT" : MAT, "MF" : MF, "MT" : MT}
    i = 0
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2, "NL" : C.L1, "NZ" : C.L2, "LRFLAG" : C.N1, "NGN" : C.N2})
    groups = {}
    for ig in range(out["NGN"]):
        L, i = read_list(str_list, i)
        group = {"TEMPIN" : L.C1, "NG2" : L.L1, "IG2LO" : L.L2, "IG" : L.N2, "DATA" : L.B}
        groups[L.N2] = group
    out["GROUPS"] = groups
    return out