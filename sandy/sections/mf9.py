"""
This module contains only two public functions:

    * `read_mf9`
    * `write_mf9`

Function `read` reads a MF9/MT section from a string and produces a content
object with a dictionary-like structure.
The content object can be accessed using most of the keywords specified in
the ENDF6 manual for this specific MF section.

Function `write_mf9` writes a content object for a MF9/MT section into a
string.
MAT, MF, MT and line numbers are also added (each line ends with a `\n`).
"""
__author__ = "Aitor Bengoechea"
__all__ = [
        "read_mf9",
        "write_mf9",
        ]

import sandy


def read_mf8(tape, mat, mt):
    """
    Parse MAT/MF=9/MT section from `sandy.Endf6` object and return
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
    """
    mf = 9
    df = tape._get_section_df(mat, mf, mt)
    out = {"MAT": mat, "MF": mf, "MT": mt}
    i = 0
    C, i = sandy.read_cont(df, i)
    out.update({
        "ZA": C.C1,
        "AWR": C.C2,
        "LIS": C.L1,
        "NS": C.N1,
        })
    T, i = sandy.read_tab1(df, i)
    out.update({
            "QM": T.C1,
            "QI": T.C2,
            "IZAP": T.L1,
            "LFS": T.L2,
            "NBT": T.NBT,
            "INT": T.INT,
            "E": T.x,
            "Y": T.y,
            })
    return out


def write_mf9(sec):
    """
    Given the content of a MF9 section as nested dictionaries, write it
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
            sec["LIS"],
            0,
            sec["NS"],
            0,
            )
    lines += sandy.write_tab1(
            sec["QM"],
            sec["QI"],
            sec["IZAP"],
            sec["LFS"],
            sec["NBT"],
            sec["INT"],
            sec["E"],
            sec["Y"],
            )
    return "\n".join(sandy.write_eol(lines, sec["MAT"], 9, sec["MT"]))
