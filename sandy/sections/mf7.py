# -*- coding: utf-8 -*-
"""
This module contains only two public functions:

    * `read_mf7`
    * `write_mf7`

Function `read` reads a MF7/MT section from a string and produces a content
object with a dictionary-like structure.
The content object can be accessed using most of the keywords specified in
the ENDF6 manual for this specific MF section.

Function `write_mf7` writes a content object for a MF7/MT section into a
string.
MAT, MF, MT and line numbers are also added (each line ends with a `\n`).
"""

import sandy

__author__ = "Aitor Bengoechea"
__all__ = [
        "read_mf7",
        "write_mf7",
        ]

mf = 7


def read_mf7(tape, mat, mt):
    """
    Parse MAT/MF=MAT/7 section from `sandy.Endf6` object and return
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
    if mt == 2:
        out = _read_elastic_scattering(tape, mat, mt)
    elif mt == 4:
        out = _read_incoherent_inelastic(tape, mat, mt)
    return out


def write_mf7(sec):
    """
    Write MT section for MF7

    Parameters
    ----------
    sec : 'dict'
        Multiline string reproducing the content of a ENDF-6 section.

    Returns
    -------
    `str`
    """
    if sec["MT"] == 2:
        return _write_elastic_scattering(sec)
    elif sec["MT"] == 4:
        return _write_inelastic_scattering(sec)


def _read_elastic_scattering(tape, mat, mt):
    """
    Parse MAT/MF=MAT/7 section for tsl from `sandy.Endf6` object
    and return structured content in nested dictionaries for elastic coherent
    and incoherent scattering.

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
    out : `dict`
        Content of the ENDF-6 tape structured as nested `dict`.

    Examples
    --------
    Incoherent
    >>> tls = sandy.get_endf6_file("endfb_80", 'tsl', 10)
    >>> sandy.sections.mf7._read_elastic_scattering(tls, 10, 2)
    {'MAT': 10,
     'MF': 7,
     'MT': 2,
     'ZA': 110.0,
     'AWR': 0.99928,
     'LTHR': 2,
     'SB': 80.31784,
     'NBT': [9],
     'INT': [2],
     'TINT': array([115.  , 188.15, 208.15, 228.15, 233.15, 248.15, 253.15, 268.15,
            273.15]),
     'W': array([14.70372, 19.1224 , 20.37892, 21.65261, 21.97355, 22.94205,
            23.26671, 24.24591, 24.57398])}

    Coherent
    >>> tls = sandy.get_endf6_file("endfb_80", 'tsl', 26)
    >>> sandy.sections.mf7._read_elastic_scattering(tls, 26, 2)['T'].keys()
    dict_keys([296.0, 400.0, 500.0, 600.0, 700.0, 800.0, 1000.0, 1200.0])
    """
    df = tape._get_section_df(mat, mf, mt)
    out = {
            "MAT": mat,
            "MF": mf,
            "MT": mt,
            }
    i = 0
    C, i = sandy.read_cont(df, i)
    LTHR = C.L1
    add = {
        "ZA": C.C1,
        "AWR": C.C2,
        "LTHR": LTHR,
        }
    out.update(add)
    if LTHR == 1:  # Coherent
        add_temp = {}
        C, i = sandy.read_tab1(df, i)
        temp = C.C1  # Temperature
        LT = C.L1  # Temperature flag
        add_2 = {
            "LT": LT,  # Flag for temperature dependence
            "NBT": C.NBT,  # Number of different pairs alpha S
            "INT": C.INT,  # Number of Bragg edges
            "S": C.y,  # S matrix values for a all the energy but only 1 temp
            }
        add["EINT"] = C.x  # Energy array, constant
        add_temp[temp] = add_2
        for j in range(LT):
            C, i = sandy.read_list(df, i)
            temp = C.C1 # Temperature
            add_2 = {
                "LI": C.L1,  # Flag indicating how to interpolate
                "S": C.B,  # S matrix values for energy array but only 1 temp
                }
            add_temp[temp] = add_2
        add['T'] = add_temp
    elif LTHR == 2:  # Incoherent
        C, i = sandy.read_tab1(df, i)
        add = {
            'SB': C.C1,  # characteristic bound cross section (barns)
            'NBT': C.NBT,
            'INT': C.INT,
            'TINT': C.x,  # Temperature array
            'W': C.y  # Debye-Waller integral divided by the atomic mass(eV^-1)
            }
    out.update(add)
    return out


def _read_incoherent_inelastic(tape, mat, mt):
    """
    Parse MAT/MF=MAT/7 section for tsl from `sandy.Endf6` object
    and return structured content in nested dictionaries for incoherent
    inelastic scattering.

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
    out : `dict`
        Content of the ENDF-6 tape structured as nested `dict`.

    Coherent
    >>> tls = sandy.get_endf6_file("endfb_80", 'tsl', 26)
    >>> dict = sandy.sections.mf7._read_incoherent_inelastic(tls, 26, 4)
    >>> dict['BN']
    [6.153875, 197.6285, 8.93478, 5.000001, 0.0, 1.0]

    Temperature for the first beta:
    >>> dict['beta/T'][0.0].keys()
    dict_keys([296.0, 400.0, 500.0, 600.0, 700.0, 800.0, 1000.0, 1200.0])

    Effective temperature:
    >>> dict['effective T'][0]['T_eff']
    array([ 433.3817,  506.3929,  586.9472,  673.3305,  763.3208,  855.6755,
           1044.799 , 1237.451 ])
    """
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
            "LAT": C.L2, # Flag indicating which temperature has been used to compute alpha and beta
            "LASYM": C.N1  # Symmetry or asymmetry of S matrix
            }
    out.update(add)
    C, i = sandy.read_list(df, i)
    NS = C.N2
    B = C.B
    add = {
        'LLN': C.L1,  # S values form: direct or ln(S)
        'NS': NS,  # Number of non-principal scattering atom types
        'BN': C.B,  # List of constants
        }
    out.update(add)
    if B[0] != 0:
        C, i = sandy.read_tab2(df, i)
        add = {
                'NR': C.NZ,  # Number of interpolation ranges for alpha and beta
                'NP': C.NBT,  # Number of alphas
                'BINT': C.INT,  # Interpolation schemes used.
        }
        Num_beta = C.NBT[0]
        add_beta = {}
        for j in range(Num_beta):
            add_temp = {}
            T, i = sandy.read_tab1(df, i)
            temp = T.C1  # Temperature
            b = T.C2  # beta
            LT = T.L1  # Temperature flag
            add_2 = {
                "LT": LT,  # Temperature flag
                "INT": T.INT,  # Number of interpolation ranges
                "NBT": T.NBT,  # Number of different pairs alpha S
                "S": T.y,  # S matrix values for all alpha in a given beta and temp 
                }
            add["alpha"] = T.x  # Alpha values, constant values
            add_temp[temp] = add_2
            for z in range(LT):
                C, i = sandy.read_list(df, i)
                temp = C.C1  # Temperature
                add_2 = {
                    "LI": C.L1,  # Flag indicating how to interpolate
                    "S": C.B,  # S matrix values for all alpha in a given beta and temp
                    }
                add_temp[temp] = add_2
            add_beta[b] = add_temp
        add['beta/T'] = add_beta
        # Effective temperature:
        add_efective_temp = {}
        C, i = sandy.read_tab1(df, i)
        add_2 = ({
                "NR": C.NBT,  # Number of interpolation ranges for alpha and beta
                "NP": C.INT,  # Number of temperatures
                "TINT": C.x,  # temperatures
                "T_eff": C.y,  # Table of effective temperatures (K) for the shortcollision-time approximation
                })
        add_efective_temp[0] = add_2
        add['effective T'] = add_efective_temp
        if NS >= 1 and B[6] == 0:
            C, i = sandy.read_tab1(df, i)
            add_2 = ({
                "NR": C.NBT,  # Number of interpolation ranges for alpha and beta
                "NP": C.INT,  # Number of temperatures
                "TINT": C.x,  # temperatures (K)
                "T_eff": C.y,  # Table of effective temperatures (K) for the shortcollision-time approximation
                })
            add['effective T'][1] = add_2
            if NS >= 2 and B[12] == 0:
                C, i = sandy.read_tab1(df, i)
                add_2 = ({
                        "NR": C.NBT,  # Number of interpolation ranges for alpha and beta
                        "NP": C.INT,  # Number of temperatures
                        "TINT": C.x,  # temperatures
                        "T_eff": C.y,  # Table of effective temperatures (K) for the shortcollision-time approximation
                        })
                add['effective T'][2] = add_2
            if NS == 3 and B[18] == 0:
                C, i = sandy.read_tab1(df, i)
                add_2 = ({
                        "NR": C.NBT,  # Number of interpolation ranges for alpha and beta
                        "NP": C.INT,  # Number of temperatures
                        "TINT": C.x,  # temperatures (K)
                        "T_eff": C.y,  # Table of effective temperatures (K) for the shortcollision-time approximation
                        })
                add['effective T'][3] = add_2
        out.update(add)
    return out


def _write_elastic_scattering(sec):
    """
    Given the content of MF7 for elastic scattering as nested dictionaries,
    write it to string.

    Parameters
    ----------
    sec : 'dict'
        Content of the ENDF-6 tape structured as nested `dict`.

    Returns
    -------
    `str`
        Multiline string reproducing the content of a ENDF-6 section.

    Examples
    --------
    Coherent elastic scattering:
    >>> tls = sandy.get_endf6_file("endfb_80", 'tsl', 26)
    >>> sec = sandy.sections.mf7._read_elastic_scattering(tls, 26, 2)
    >>> assert len(_write_elastic_scattering(sec)) == len(tls.data[(26, 7, 2)])

    Incoherent elastic scattering:
    >>> tls = sandy.get_endf6_file("endfb_80", 'tsl', 10)
    >>> sec = sandy.sections.mf7._read_elastic_scattering(tls, 10, 2)
    >>> print(_write_elastic_scattering(sec))
     110.000000 9.992800-1          2          0          0          0  10 7  2    1
     80.3178400 0.00000000          0          0          1          9  10 7  2    2
              9          2                                              10 7  2    3
     115.000000 14.7037200 188.150000 19.1224000 208.150000 20.3789200  10 7  2    4
     228.150000 21.6526100 233.150000 21.9735500 248.150000 22.9420500  10 7  2    5
     253.150000 23.2667100 268.150000 24.2459100 273.150000 24.5739800  10 7  2    6
    """
    if sec['LTHR'] == 1:
        lines = sandy.write_cont(sec["ZA"],
                                 sec["AWR"],
                                 sec['LTHR'],
                                 0,
                                 0,
                                 0)
        diff = list(sec['T'].keys())[0]
        for T, T_info in sec['T'].items():
            if T == diff:
                lines += sandy.write_tab1(T,
                                          0,
                                          T_info['LT'],
                                          0,
                                          T_info['NBT'],
                                          T_info['INT'],
                                          sec['EINT'],
                                          T_info['S'])
            else:
                lines += sandy.write_list(T,
                                          0,
                                          T_info['LI'],
                                          0,
                                          0,
                                          T_info['S'])
    elif sec['LTHR'] == 2:
        lines = sandy.write_cont(sec["ZA"], sec["AWR"], sec['LTHR'], 0, 0, 0)
        lines += sandy.write_tab1(sec['SB'],
                                  0,
                                  0,
                                  0,
                                  sec['NBT'],
                                  sec['INT'],
                                  sec['TINT'],
                                  sec['W'])
    return "\n".join(sandy.write_eol(lines, sec["MAT"], 7, sec["MT"]))


def _write_inelastic_scattering(sec):
    """
    Given the content of MF7 for inelastic scattering as nested dictionaries,
    write it to string.

    Parameters
    ----------
    sec : 'dict'
        Content of the ENDF-6 tape structured as nested `dict`.

    Returns
    -------
    `str`
        Multiline string reproducing the content of a ENDF-6 section.

    Examples
    --------
    Incoherent inelastic scattering:
    >>> tls = sandy.get_endf6_file("endfb_80", 'tsl', 26)
    >>> sec = sandy.sections.mf7._read_incoherent_inelastic(tls, 26, 4)
    >>> assert len(_write_inelastic_scattering(sec)) == len(tls.data[(26, 7, 4)])
    """
    lines = sandy.write_cont(sec["ZA"],
                             sec["AWR"],
                             0,
                             sec["LAT"],
                             sec["LASYM"],
                             0)
    lines += sandy.write_list(0,
                              0,
                              sec['LLN'],
                              0,
                              sec['NS'],
                              sec['BN'])
    if sec['BN'][0] != 0:
        lines += sandy.write_tab2(0,
                                  0,
                                  0,
                                  0,
                                  sec['NR'],
                                  sec['NP'],
                                  sec['BINT'])
        for beta, beta_info in sec['beta/T'].items():
            diff_T = list(beta_info.keys())[0]
            for T, T_info in beta_info.items():
                if diff_T == T:
                    lines += sandy.write_tab1(T,
                                              beta,
                                              T_info['LT'],
                                              0,
                                              T_info['NBT'],
                                              T_info['INT'],
                                              sec['alpha'],
                                              T_info['S'])
                else:
                    lines += sandy.write_list(T,
                                              beta,
                                              T_info['LI'],
                                              0,
                                              0,
                                              T_info['S'])
        for T_eff_numb, T_eff_info in sec['effective T'].items():
            lines += sandy.write_tab1(0,
                                      0,
                                      0,
                                      0,
                                      T_eff_info['NR'],
                                      T_eff_info['NP'],
                                      T_eff_info['TINT'],
                                      T_eff_info['T_eff'])
    return "\n".join(sandy.write_eol(lines, sec["MAT"], 7, sec["MT"]))
