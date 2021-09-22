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


def read_mf10(tape, mat, mt):
    """
    Parse MAT/MF=10/MT section from `sandy.Endf6` object and return
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
    >>>read_mf10(sandy.get_endf6_file("jeff_40t0", 'xs', 410930),4125,16)
    {'MAT': 4125,
     'MF': 10,
     'MT': 16,
     'ZA': 41093.0,
     'AWR': 92.10827,
     'LIS': 0,
     'NS': 2,
     'subsection1': {'QM': -8830870.0,
                     'QI': -8830870.0,
                     'IZAP': 41092,
                     'LFS': 0,
                     'NR': [36],
                     'NP': [2],
                     'E': array([8.92674e+06, 9.00000e+06, 9.50000e+06, 1.00000e+07, 1.05000e+07,
                                 1.10000e+07, 1.15000e+07, 1.20000e+07, 1.25000e+07, 1.30000e+07,
                                 1.35000e+07, 1.40000e+07, 1.45000e+07, 1.50000e+07, 1.55000e+07,
                                 1.60000e+07, 1.65000e+07, 1.70000e+07, 1.75000e+07, 1.80000e+07,
                                 1.85000e+07, 1.90000e+07, 1.95000e+07, 2.00000e+07, 2.10000e+07,
                                 2.20000e+07, 2.30000e+07, 2.40000e+07, 2.50000e+07, 2.60000e+07,
                                 2.70000e+07, 2.80000e+07, 2.90000e+07, 3.00000e+07, 3.00000e+07,
                                 2.00000e+08]),
                     'XS': array([0.        , 0.0194457 , 0.08620951, 0.273596  , 0.447849  ,
                                  0.546337  , 0.600883  , 0.6416661 , 0.663568  , 0.672792  ,
                                  0.678453  , 0.689907  , 0.70193   , 0.716936  , 0.722198  ,
                                  0.733738  , 0.7491971 , 0.7579651 , 0.766395  , 0.759905  ,
                                  0.746675  , 0.7225931 , 0.69735   , 0.656834  , 0.574932  ,
                                  0.493704  , 0.429843  , 0.385963  , 0.348244  , 0.322294  ,
                                  0.299083  , 0.281999  , 0.263989  , 0.253054  , 0.        ,
                                  0.        ])},
     'subsection2': {'QM': -8830870.0,
                     'QI': -8966370.0,
                     'IZAP': 41092,
                     'LFS': 1,
                     'NR': [35],
                     'NP': [2],
                     'E': array([9.06371e+06, 9.50000e+06, 1.00000e+07, 1.05000e+07, 1.10000e+07,
                                 1.15000e+07, 1.20000e+07, 1.25000e+07, 1.30000e+07, 1.35000e+07,
                                 1.40000e+07, 1.45000e+07, 1.50000e+07, 1.55000e+07, 1.60000e+07,
                                 1.65000e+07, 1.70000e+07, 1.75000e+07, 1.80000e+07, 1.85000e+07,
                                 1.90000e+07, 1.95000e+07, 2.00000e+07, 2.10000e+07, 2.20000e+07,
                                 2.30000e+07, 2.40000e+07, 2.50000e+07, 2.60000e+07, 2.70000e+07,
                                 2.80000e+07, 2.90000e+07, 3.00000e+07, 3.00000e+07, 2.00000e+08]),
                     'XS': array([0.       , 0.0312977, 0.127211 , 0.213139 , 0.285268 , 0.344366 ,
                                  0.389636 , 0.421435 , 0.441539 , 0.452998 , 0.459128 , 0.46013  ,
                                  0.460299 , 0.45819  , 0.45204  , 0.44854  , 0.444947 , 0.433968 ,
                                  0.420582 , 0.395455 , 0.371834 , 0.346106 , 0.323454 , 0.272154 ,
                                  0.232743 , 0.206621 , 0.186896 , 0.172196 , 0.160433 , 0.15069  ,
                                  0.142175 , 0.135358 , 0.129102 , 0.       , 0.       ])}}   
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
                "NS": C.N1,
    })
    for hx in range(C.N1):
        add = {}
        T, i = sandy.read_tab1(df, i)
        add = {
                    "QM": T.C1,
                    "QI": T.C2,
                    "IZAP": T.L1,
                    "LFS": T.L2,
                    "NR": T.NBT,
                    "NP": T.INT,
                    "E": T.x,
                    "XS": T.y,
                    }
        out["subsection" + str(hx+1)] = add
    return out


def write_mf10(sec):
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

    Examples
    --------
    String reproducing the content of a ENDF-6 section for (z,2n)
    of Nb-93 from the JEFF-4.0T0 library to obtain the cross sections for
    production of radiactive nuclides
    >>> write_mf10(read_mf10(sandy.get_endf6_file("jeff_40t0", 'xs', 410930),4125,16))
    ' 41093.0000 92.1082700          0          0          2          0412510 16    1\n-8830870.00-8830870.00      41092          0          1         36412510 16    2\n         36          2                                            412510 16    3\n 8926740.00 0.00000000 9000000.00 1.944570-2 9500000.00 8.620951-2412510 16    4\n 10000000.0 2.735960-1 10500000.0 4.478490-1 11000000.0 5.463370-1412510 16    5\n 11500000.0 6.008830-1 12000000.0 6.416661-1 12500000.0 6.635680-1412510 16    6\n 13000000.0 6.727920-1 13500000.0 6.784530-1 14000000.0 6.899070-1412510 16    7\n 14500000.0 7.019300-1 15000000.0 7.169360-1 15500000.0 7.221980-1412510 16    8\n 16000000.0 7.337380-1 16500000.0 7.491971-1 17000000.0 7.579651-1412510 16    9\n 17500000.0 7.663950-1 18000000.0 7.599050-1 18500000.0 7.466750-1412510 16   10\n 19000000.0 7.225931-1 19500000.0 6.973500-1 20000000.0 6.568340-1412510 16   11\n 21000000.0 5.749320-1 22000000.0 4.937040-1 23000000.0 4.298430-1412510 16   12\n 24000000.0 3.859630-1 25000000.0 3.482440-1 26000000.0 3.222940-1412510 16   13\n 27000000.0 2.990830-1 28000000.0 2.819990-1 29000000.0 2.639890-1412510 16   14\n 30000000.0 2.530540-1 30000000.0 0.00000000  200000000 0.00000000412510 16   15\n-8830870.00-8966370.00      41092          1          1         35412510 16   16\n         35          2                                            412510 16   17\n 9063710.00 0.00000000 9500000.00 3.129770-2 10000000.0 1.272110-1412510 16   18\n 10500000.0 2.131390-1 11000000.0 2.852680-1 11500000.0 3.443660-1412510 16   19\n 12000000.0 3.896360-1 12500000.0 4.214350-1 13000000.0 4.415390-1412510 16   20\n 13500000.0 4.529980-1 14000000.0 4.591280-1 14500000.0 4.601300-1412510 16   21\n 15000000.0 4.602990-1 15500000.0 4.581900-1 16000000.0 4.520400-1412510 16   22\n 16500000.0 4.485400-1 17000000.0 4.449470-1 17500000.0 4.339680-1412510 16   23\n 18000000.0 4.205820-1 18500000.0 3.954550-1 19000000.0 3.718340-1412510 16   24\n 19500000.0 3.461060-1 20000000.0 3.234540-1 21000000.0 2.721540-1412510 16   25\n 22000000.0 2.327430-1 23000000.0 2.066210-1 24000000.0 1.868960-1412510 16   26\n 25000000.0 1.721960-1 26000000.0 1.604330-1 27000000.0 1.506900-1412510 16   27\n 28000000.0 1.421750-1 29000000.0 1.353580-1 30000000.0 1.291020-1412510 16   28\n 30000000.0 0.00000000  200000000 0.00000000                      412510 16   29'
    """
    lines = sandy.write_cont(
            sec["ZA"],
            sec["AWR"],
            sec["LIS"],
            0,
            sec["NS"],
            0,
            )
    for hz in range(sec["NS"]):
        subsection = sec["subsection"+str(hz+1)]
        lines += sandy.write_tab1(
                subsection["QM"],
                subsection["QI"],
                subsection["IZAP"],
                subsection["LFS"],
                subsection["NR"],
                subsection["NP"],
                subsection["E"],
                subsection["XS"],
                )
    return "\n".join(sandy.write_eol(lines, sec["MAT"], 10, sec["MT"]))
