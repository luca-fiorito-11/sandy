"""
This module contains only two public functions:

    * `read_mf10`
    * `write_mf10`

Function `read` reads a MF10/MT section from a string and produces a content
object with a dictionary-like structure.
The content object can be accessed using most of the keywords specified in
the ENDF6 manual for this specific MF section.

Function `write_mf10` writes a content object for a MF10/MT section into a
string.
MAT, MF, MT and line numbers are also added (each line ends with a `\n`).
"""
__author__ = "Aitor Bengoechea"
__all__ = [
        "read_mf10",
        "write_mf10",
        ]

import sandy
import numpy as np

def read_mf10(tape, mat, mt):
    """
    Parse MAT/MF=10/MT section from `sandy.Endf6` object and return
    structured content in nested dictionary.

    Parameters
    ----------
    tape : `sandy.Endf6`
        endf6 object containing requested section.
    mat : `int`
        MAT number.
    mt : `int`
        MT number.

    Returns
    -------
    out : `dict`
        Content of the ENDF-6 tape structured as nested `dict`.

    Examples
    --------
    >>> tape = sandy.get_endf6_file("jeff_33", 'xs', 410930)
    >>> test = read_mf10(tape, 4125, 16)
    >>> read_mf10_test(test)
    """
    mf = 10
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
                    "XS": T.y,
                    }
        subsections[LFS] = add
    out["LFS"] = subsections
    return out


def write_mf10(sec):
    """
    Given the content of a MF9 section as nested dictionaries, write it
    to string.

    Parameters
    ----------
    sec : 'dict'
        Multiline string reproducing the content of a ENDF-6 section.

    Returns
    -------
    `str`
        Multiline string reproducing the content of a ENDF-6 section.

    Notes
    -----
    .. note:: The end-of-line records MAT, MF, MT and line number are added at
              the end of each line.

    .. important:: The string does not endf with a newline symbol `\n`.

    Examples
    --------
    String reproducing the content of a ENDF-6 section for (z,2n)
    of Nb-93 from the JEFF-33 library to obtain the cross sections for
    production of radiactive nuclides
    >>> tape = sandy.get_endf6_file("jeff_33", 'xs', 410930)
    >>> sec = read_mf10(tape, 4125, 16)
    >>> text = write_mf10(sec)
    >>> print(text)
     41093.0000 92.1082700          0          0          2          0412510 16    1
    -8830560.00-8830560.00      41092          0          1         31412510 16    2
             31          2                                            412510 16    3
     8926430.00 0.00000000 9000000.00 1.941860-2 9500000.00 1.621040-1412510 16    4
     10000000.0 3.089020-1 10500000.0 4.351160-1 11000000.0 5.174740-1412510 16    5
     11500000.0 5.739980-1 12000000.0 6.246930-1 12500000.0 6.593280-1412510 16    6
     13000000.0 6.818000-1 13500000.0 6.988180-1 14000000.0 7.179610-1412510 16    7
     14500000.0 7.367761-1 15000000.0 7.507430-1 16000000.0 7.747931-1412510 16    8
     17000000.0 7.894241-1 18000000.0 7.807111-1 19000000.0 7.326431-1412510 16    9
     20000000.0 6.574070-1 21000000.0 5.754551-1 22000000.0 4.940910-1412510 16   10
     23000000.0 4.301710-1 24000000.0 3.861880-1 25000000.0 3.484020-1412510 16   11
     26000000.0 3.224040-1 27000000.0 2.991570-1 28000000.0 2.820340-1412510 16   12
     29000000.0 2.640060-1 30000000.0 2.530630-1 30000000.0 0.00000000412510 16   13
      200000000 0.00000000                                            412510 16   14
    -8830560.00-8966060.00      41092          1          1         30412510 16   15
             30          2                                            412510 16   16
     9063400.00 0.00000000 9500000.00 3.129770-2 10000000.0 1.272110-1412510 16   17
     10500000.0 2.131390-1 11000000.0 2.852680-1 11500000.0 3.443660-1412510 16   18
     12000000.0 3.896360-1 12500000.0 4.214350-1 13000000.0 4.415390-1412510 16   19
     13500000.0 4.529980-1 14000000.0 4.591280-1 14500000.0 4.601300-1412510 16   20
     15000000.0 4.603000-1 16000000.0 4.520400-1 17000000.0 4.450740-1412510 16   21
     18000000.0 4.207860-1 19000000.0 3.720790-1 20000000.0 3.237210-1412510 16   22
     21000000.0 2.723630-1 22000000.0 2.328870-1 23000000.0 2.067260-1412510 16   23
     24000000.0 1.869740-1 25000000.0 1.722300-1 26000000.0 1.604680-1412510 16   24
     27000000.0 1.507150-1 28000000.0 1.421910-1 29000000.0 1.353680-1412510 16   25
     30000000.0 1.291160-1 30000000.0 0.00000000  200000000 0.00000000412510 16   26
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
                LFS[key[hz]]["XS"],
                )
    return "\n".join(sandy.write_eol(lines, sec["MAT"], 10, sec["MT"]))


def read_mf10_test(text):
    from numpy import array
    desired = {'MAT': 4125,
                     'MF': 10,
                     'MT': 16,
                     'ZA': 41093.0,
                     'AWR': 92.10827,
                     'LIS': 0,
                     'LFS': {0: {'QM': -8830560.0,
                                 'QI': -8830560.0,
                                 'IZAP': 41092,
                                 'NBT': [31],
                                 'INT': [2],
                                 'E': array([8.92643e+06, 9.00000e+06, 9.50000e+06, 1.00000e+07, 1.05000e+07,
                                             1.10000e+07, 1.15000e+07, 1.20000e+07, 1.25000e+07, 1.30000e+07,
                                             1.35000e+07, 1.40000e+07, 1.45000e+07, 1.50000e+07, 1.60000e+07,
                                             1.70000e+07, 1.80000e+07, 1.90000e+07, 2.00000e+07, 2.10000e+07,
                                             2.20000e+07, 2.30000e+07, 2.40000e+07, 2.50000e+07, 2.60000e+07,
                                             2.70000e+07, 2.80000e+07, 2.90000e+07, 3.00000e+07, 3.00000e+07,
                                             2.00000e+08]),
                                 'XS': array([0.       , 0.0194186, 0.162104 , 0.308902 , 0.435116 , 0.517474 ,
                                              0.573998 , 0.624693 , 0.659328 , 0.6818   , 0.698818 , 0.717961 ,
                                              0.7367761, 0.750743 , 0.7747931, 0.7894241, 0.7807111, 0.7326431,
                                              0.657407 , 0.5754551, 0.494091 , 0.430171 , 0.386188 , 0.348402 ,
                                              0.322404 , 0.299157 , 0.282034 , 0.264006 , 0.253063 , 0.       ,
                                              0.       ])},
                             1: {'QM': -8830560.0,
                                 'QI': -8966060.0,
                                 'IZAP': 41092,
                                 'NBT': [30],
                                 'INT': [2],
                                 'E': array([9.0634e+06, 9.5000e+06, 1.0000e+07, 1.0500e+07, 1.1000e+07,
                                             1.1500e+07, 1.2000e+07, 1.2500e+07, 1.3000e+07, 1.3500e+07,
                                             1.4000e+07, 1.4500e+07, 1.5000e+07, 1.6000e+07, 1.7000e+07,
                                             1.8000e+07, 1.9000e+07, 2.0000e+07, 2.1000e+07, 2.2000e+07,
                                             2.3000e+07, 2.4000e+07, 2.5000e+07, 2.6000e+07, 2.7000e+07,
                                             2.8000e+07, 2.9000e+07, 3.0000e+07, 3.0000e+07, 2.0000e+08]),
                                 'XS': array([0.       , 0.0312977, 0.127211 , 0.213139 , 0.285268 , 0.344366 ,
                                              0.389636 , 0.421435 , 0.441539 , 0.452998 , 0.459128 , 0.46013  ,
                                              0.4603   , 0.45204  , 0.445074 , 0.420786 , 0.372079 , 0.323721 ,
                                              0.272363 , 0.232887 , 0.206726 , 0.186974 , 0.17223  , 0.160468 ,
                                              0.150715 , 0.142191 , 0.135368 , 0.129116 , 0.       , 0.       ])}}}
    return np.testing.assert_equal(text,desired, err_msg='The program did not pass the test')
