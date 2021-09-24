"""
This module contains only two public functions:
<<<<<<< HEAD
    * `read_mf9`
    * `write_mf9`
=======

    * `read_mf9`
    * `write_mf9`

>>>>>>> ce08eca97e010ab1996b71348064829541a5facf
Function `read` reads a MF9/MT section from a string and produces a content
object with a dictionary-like structure.
The content object can be accessed using most of the keywords specified in
the ENDF6 manual for this specific MF section.
<<<<<<< HEAD
=======

>>>>>>> ce08eca97e010ab1996b71348064829541a5facf
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


def read_mf9(tape, mat, mt):
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
    out: `dict`
        Content of the ENDF-6 tape structured as nested `dict`.

    Examples
    --------
    Endf-6 tape structured 'dict' of Radiactive capture of Am-241 from
    the ENDFB-VII.1 library to obtain the multiplicities for production
    of radioactive nuclides
    >>> tape = sandy.get_endf6_file("endfb_71", 'xs', 952410)
    >>> read_mf9(tape, 9543, 102)
    {'MAT': 9543,
     'MF': 9,
     'MT': 102,
     'ZA': 95241.0,
     'AWR': 238.986,
     'LIS': 0,
     'LFS': {0: {'QM': 5539101.0,
       'QI': 5539101.0,
       'IZAP': 95242,
       'NBT': [9],
       'INT': [3],
       'E': array([1.000000e-05, 3.690000e-01, 1.000000e+03, 1.000000e+05,
              6.000001e+05, 1.000000e+06, 2.000000e+06, 4.000001e+06,
              3.000000e+07]),
       'Y': array([0.9    , 0.9    , 0.8667 , 0.842  , 0.81533, 0.74382, 0.5703 ,
              0.52   , 0.52   ])},
      2: {'QM': 5539101.0,
       'QI': 5490471.0,
       'IZAP': 95242,
       'NBT': [9],
       'INT': [3],
       'E': array([1.000000e-05, 3.690000e-01, 1.000000e+03, 1.000000e+05,
              6.000001e+05, 1.000000e+06, 2.000000e+06, 4.000001e+06,
              3.000000e+07]),
       'Y': array([0.1    , 0.1    , 0.1333 , 0.158  , 0.18467, 0.25618, 0.4297 ,
              0.48   , 0.48   ])}}}
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
        })
    subsections = {}
    for hx in range(C.N1):
        T, i = sandy.read_tab1(df, i)
        LFS = T.L2
        add = {
                "QM": T.C1,
                "QI": T.C2,
                "IZAP": T.L1,
                "NBT": T.NBT,
                "INT": T.INT,
                "E": T.x,
                "Y": T.y,
              }
        subsections[LFS] = add
    out["LFS"] = subsections
    return out


def write_mf9(sec):
    """
    Given the content of a MF9 section as nested dictionaries, write it
    to string.

    Parameters
    ----------
    sec : 'dic'
        Content of the ENDF-6 tape structured as nested `dict`.

    Notes
    -----
    .. note:: The end-of-line records MAT, MF, MT and line number are added at
              the end of each line.

    .. important:: The string does not endf with a newline symbol `\n`.

    Examples
    --------
    String reproducing the content of a ENDF-6 section for Radiactive capture
    of Am-241 from the ENDFB-VII.1 library to obtain the multiplicities for
    production of radioactive nuclides
    >>> tape = sandy.get_endf6_file("endfb_71", 'xs', 952410)
    >>> sec = read_mf9(tape, 9543, 102)
    >>> text = write_mf9(sec)
    >>> print(text)
     95241.0000 238.986000          0          0          2          09543 9102    1
     5539101.00 5539101.00      95242          0          1          99543 9102    2
              9          3                                            9543 9102    3
     1.000000-5 9.000000-1 3.690000-1 9.000000-1 1000.00000 8.667000-19543 9102    4
     100000.000 8.420000-1 600000.100 8.153300-1 1000000.00 7.438200-19543 9102    5
     2000000.00 5.703000-1 4000001.00 5.200000-1 30000000.0 5.200000-19543 9102    6
     5539101.00 5490471.00      95242          2          1          99543 9102    7
              9          3                                            9543 9102    8
     1.000000-5 1.000000-1 3.690000-1 1.000000-1 1000.00000 1.333000-19543 9102    9
     100000.000 1.580000-1 600000.100 1.846700-1 1000000.00 2.561800-19543 9102   10
     2000000.00 4.297000-1 4000001.00 4.800000-1 30000000.0 4.800000-19543 9102   11
    """

    lines = sandy.write_cont(
            sec["ZA"],
            sec["AWR"],
            sec["LIS"],
            0,
            len(sec["LFS"]),
            0,
            )

    for hz in range(len(sec["LFS"])):
        LFS = sec["LFS"]
        key = list(LFS.keys())
        lines += sandy.write_tab1(
                LFS[key[hz]]["QM"],
                LFS[key[hz]]["QI"],
                LFS[key[hz]]["IZAP"],
                key[hz],
                LFS[key[hz]]["NBT"],
                LFS[key[hz]]["INT"],
                LFS[key[hz]]["E"],
                LFS[key[hz]]["Y"],
                )
    return "\n".join(sandy.write_eol(lines, sec["MAT"], 9, sec["MT"]))
