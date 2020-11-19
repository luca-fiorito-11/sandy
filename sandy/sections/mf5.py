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
