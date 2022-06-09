"""
"This module contains only two public functions:

    * `read_mf5`
    * `write_mf5`

Function `read` reads a MF5/MT section from a string and produces a content
object with a dictionary-like structure.
The content object can be accessed using most of the keywords specified in
the ENDF6 manual for this specific MF section.

Function `write_mf5` writes a content object for a MF5/MT section into a
string.
MAT, MF, MT and line numbers are also added (each line ends with a `\n`).
"""
import sandy

__author__ = "Luca Fiorito"
__all__ = [
        "read_mf5",
        ]

allowed_mt = (
        451,
        452,
        455,
        456,
        458,
        )
mf = 5


def read_mf5(tape, mat, mt):
    """
    Parse MAT/MF=5/MT section from `sandy.Endf6` object and return structured
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
    df = tape._get_section_df(mat, mf, mt)
    out = {
            "MAT": mat,
            "MF": mf,
            "MT": mt,
            }
    i = 0
    C, i = sandy.read_cont(df, i)
    # Number of partial energy distributions. There will be one subsection
    # for each partial distribution.
    NK = C.N1
    add = {
            "ZA": C.C1,
            "AWR": C.C2,
            }
    out.update(add)
    pdistr = {}
    for j in range(NK):
        tp, i = sandy.read_tab1(df, i)
        # Flag specifying the energy distribution law used for a particular
        # subsection (partial energy distribution)
        LF = tp.L2
        # Fractional part of the particular cross section which can be
        # described by the kth partial energy distribution at the i-th
        # incident energy point
        sub = {
                "NBT_P": tp.NBT,
                "INT_P": tp.INT,
                "E_P": tp.x,
                "P": tp.y,
                "LF": LF,
                }
        # General Evaporation Spectrum (LF=5)
        if LF == 5:
            """
            Found in:
                100-Fm-255g.jeff33 (x6)
                88-Ra-226g.jeff33 (x6)
                91-Pa-233g.jeff33 (x6)
                92-U-239g.jeff33
                92-U-240g.jeff33
            """
            sub["U"] = tp.C1
            T, i = sandy.read_tab1(df, i)
            sub["NBT_THETA"] = T.NBT
            sub["INT_THETA"] = T.INT
            sub["E_THETA"] = T.x
            sub["THETA"] = T.y
            T, i = sandy.read_tab1(df, i)
            sub["NBT_G"] = T.NBT
            sub["INT_G"] = T.INT
            sub["E_G"] = T.x
            sub["G"] = T.y
        # Simple Maxwellian Fission Spectrum (LF=7) /
        # Evaporation Spectrum (LF=9)
        elif LF in (7, 9):
            """
            Found in:
                27-Co-59g.jeff33
            """
            sub["U"] = tp.C1
            T, i = sandy.read_tab1(df, i)
            sub["NBT_THETA"] = T.NBT
            sub["INT_THETA"] = T.INT
            sub["E_THETA"] = T.x
            sub["THETA"] = T.y
        # Energy-Dependent Watt Spectrum (LF=11)
        elif LF == 11:
            sub["U"] = tp.C1
            T, i = sandy.read_tab1(df, i)
            sub["NBT_A"] = T.NBT
            sub["INT_A"] = T.INT
            sub["E_A"] = T.x
            sub["A"] = T.y
            T, i = sandy.read_tab1(df, i)
            sub["NBT_B"] = T.NBT
            sub["INT_B"] = T.INT
            sub["E_B"] = T.x
            sub["B"] = T.y
        # Energy-Dependent Fission Neutron Spectrum (Madland and Nix) (LF=12)
        elif LF == 12:
            TM, i = sandy.read_tab1(df, i)
            sub["EFL"] = T.C1
            sub["EHL"] = T.C2
            sub["NBT_TM"] = T.NBT
            sub["INT_TM"] = T.INT
            sub["E_TM"] = T.x
            sub["TM"] = T.y
        # Arbitrary Tabulated Function (LF=1)
        elif LF == 1:
            T2, i = sandy.read_tab2(df, i)
            NZ = T2.NZ  # number of incident energies for which distr. is given
            sub["NBT_EIN"] = T2.NBT
            sub["INT_EIN"] = T2.INT
            edistr = {}
            for k in range(NZ):
                T1, i = sandy.read_tab1(df, i)
                e_in = T1.C2
                edistr[e_in] = {
                        "EOUT": T1.x,
                        "EDISTR": T1.y,
                        "NBT": T1.NBT,
                        "INT": T1.INT,
                        }
            sub["EIN"] = edistr
        pdistr[j] = sub
    if pdistr:
        out["PDISTR"] = pdistr
    return out


def write_mf5(sec):
    """
    Given the content of MF5 for energy distributions of secondary particles
    as nested dictionaries, write it to string.
    Parameters
    ----------
    sec : 'dict'
        Multiline string reproducing the content of a ENDF-6 section.
    Returns
    -------
    `str`
        Multiline string reproducing the content of a ENDF-6 section.
    Examples
    --------
    >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 922380)
    >>> sec = sandy.sections.read_mf5(endf6, 9237, 18)
    >>> print(write_mf5(sec)[:1000])
    92238.0000 236.005800          0          0          1          09237 5 18    1
    0.00000000 0.00000000          0          1          1          29237 5 18    2
             2          2                                            9237 5 18    3
    1.000000-5 1.00000000 30000000.0 1.00000000                      9237 5 18    4
    0.00000000 0.00000000          0          0          1         969237 5 18    5
            96          2                                            9237 5 18    6
    0.00000000 1.000000-5          0          0          1        3039237 5 18    7
           303          2                                            9237 5 18    8
    1.000000-5 1.95793-12 2.000000-5 2.76893-12 4.000000-5 3.91586-129237 5 18    9
    6.000000-5 4.79593-12 8.000000-5 5.53787-12 1.000000-4 6.19152-129237 5 18   10
    2.000000-4 8.75614-12 4.000000-4 1.23831-11 6.000000-4 1.51661-119237 5 18   11
    8.000000-4 1.75123-11 1.000000-3 1.95793-11 2.000000-3 2.76893-119237 5 18   12
    4.000000-3 3.91586-11 6.000
    """
    NK = len(sec["PDISTR"])
    lines = sandy.write_cont(sec["ZA"], sec["AWR"], 0, 0, NK, 0)
    for k, sub_info in sec["PDISTR"].items():
        LF = sub_info['LF']
        if LF not in [5, 7, 9, 11, 12]:
            lines += sandy.write_tab1(0,
                                      0,
                                      0,
                                      LF,
                                      sub_info["NBT_P"],
                                      sub_info["INT_P"],
                                      sub_info["E_P"],
                                      sub_info["P"])
        else:
            lines += sandy.write_tab1(sub_info['U'],  # C1
                                      0,  # C2
                                      0,  # L1
                                      LF,  # L2
                                      sub_info["NBT_P"],  # NBT
                                      sub_info["INT_P"],  # INT
                                      sub_info["E_P"],  # X
                                      sub_info["P"])  # Y
        if LF == 5:
            lines += sandy.write_tab1(0,  # C1
                                      0,  # C2
                                      0,  # L1
                                      0,  # L2
                                      sub_info["NBT_THETA"],  # NBT
                                      sub_info["INT_THETA"],  # INT
                                      sub_info["E_THETA"],  # X
                                      sub_info["THETA"])  # Y
            lines += sandy.write_tab1(0,  # C1
                                      0,  # C2
                                      0,  # L1
                                      0,  # L2
                                      sub_info["NBT_G"],  # NBT
                                      sub_info["INT_G"],  # INT
                                      sub_info["E_G"],  # X
                                      sub_info["G"])  # Y
        elif LF in (7, 9):
            lines += sandy.write_tab1(0,  # C1
                                      0,  # C2
                                      0,  # L1
                                      0,  # L2
                                      sub_info["NBT_THETA"],  # NBT
                                      sub_info["INT_THETA"],  # INT
                                      sub_info["E_THETA"],  # X
                                      sub_info["THETA"])  # Y
        elif LF == 11:
            lines += sandy.write_tab1(0,  # C1
                                      0,  # C2
                                      0,  # L1
                                      0,  # L2
                                      sub_info["NBT_A"],  # NBT
                                      sub_info["INT_A"],  # INT
                                      sub_info["E_A"],  # X
                                      sub_info["A"])  # Y
            lines += sandy.write_tab1(0,  # C1
                                      0,  # C2
                                      0,  # L1
                                      0,  # L2
                                      sub_info["NBT_B"],  # NBT
                                      sub_info["INT_B"],  # INT
                                      sub_info["E_B"],  # X
                                      sub_info["B"])  # Y
        elif LF == 1:
            NZ = len(sub_info['EIN'])
            lines += sandy.write_tab2(0,  # C1
                                      0,  # C2
                                      0,  # L1
                                      0,  # L2
                                      NZ,  # NZ
                                      sub_info["NBT_EIN"],  # NBT
                                      sub_info["INT_EIN"])  # INT
            for EIN, edistr in sub_info['EIN'].items():
                lines += sandy.write_tab1(0,  # C1
                                          EIN,  # C2
                                          0,  # L1
                                          0,  # L2
                                          edistr["NBT"],  # NBT
                                          edistr["INT"],  # INT
                                          edistr["EOUT"],  # X
                                          edistr["EDISTR"])  # Y
    return "\n".join(sandy.write_eol(lines, sec["MAT"], 5, sec["MT"]))
