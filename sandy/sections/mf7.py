# -*- coding: utf-8 -*-
"""
This module contains only two public functions:

    * `read_mf7`
    * `write_mf7`

Function `read` reads a MF7/MT section from a string and produces a content
object with a dictionary-like structure.
The content object can be accessed using most of the keywords specified in
the ENDF6 manual for this specific MF section.

Function `write_mf8` writes a content object for a MF7/MT section into a
string.
MAT, MF, MT and line numbers are also added (each line ends with a `\n`).
"""

import sandy

__author__ = "Aitor Bengoechea"


mf = 7


def read_mf7(tape, mat, mt):
    """
    Parse MAT/MF=MAT/8 section from `sandy.Endf6` object and return
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


def _read_elastic_scattering(tape, mat, mt):
    """
    Parse MAT/MF=MAT/8 section for tsl from `sandy.Endf6` object
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
     'ZA': 110,
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
    >>> sandy.sections.mf7._read_elastic_scattering(tls, 26, 2)['temperature'].keys()
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
        "ZA": int(C.C1),
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
            "LT": LT,
            "NBT": C.NBT,  # Number of different pairs alpha S
            "S": C.y,
            }
        add["EINT"] = C.x  # Energy array, constant
        add_temp[temp] = add_2
        for j in range(LT):
            C, i = sandy.read_list(df, i)
            temp = C.C1
            add_2 = {
                "LI": C.L1,
                "S": C.B,
                }
            add_temp[temp] = add_2
        add['temperature'] = add_temp
    if LTHR == 2:  # Incoherent
        C, i = sandy.read_tab1(df, i)
        add = {
            'SB': C.C1,
            'NBT': C.NBT,
            'INT': C.INT,
            'TINT': C.x,
            'W': C.y
            }
    out.update(add)
    return out


def _read_incoherent_inelastic(tape, mat, mt):
    """
    Parse MAT/MF=MAT/8 section for tsl from `sandy.Endf6` object
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
    >>> dict['beta/temperature'][0.0].keys()
    dict_keys([296.0, 400.0, 500.0, 600.0, 700.0, 800.0, 1000.0, 1200.0])

    Effective temperature:
    >>> dict['effective temp'][0]['T_eff']
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
            "ZA": int(C.C1),
            "AWR": C.C2,
            "LAT": C.L2,
            "LASYM": C.N1
            }
    C, i = sandy.read_list(df, i)
    NS = C.N2
    B = C.B
    add = {
        'LLN': C.L1,
        'NS': NS,
        'BN': C.B,
        }
    out.update(add)
    if B[0] != 0:
        C, i = sandy.read_tab2(df, i)
        add = {
                'NR': C.NZ,
                'NP': C.NBT,
                'BINT': C.INT,
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
                "LT": LT,
                "NBT": T.NBT,  # Number of different pairs alpha S
                "S": T.y,
                }
            add["alpha"] = T.x  # Alpha values, constant values
            add_temp[temp] = add_2
            for z in range(LT):
                C, i = sandy.read_list(df, i)
                temp = C.C1  # Temperature
                add_2 = {
                    "LI": C.L1,
                    "S": C.B,
                    }
                add_temp[temp] = add_2
            add_beta[b] = add_temp
        add['beta/temperature'] = add_beta
        add_efective_temp = {}
        C, i = sandy.read_tab1(df, i)
        add_2.update({
                "NR": C.NBT,
                "NP": C.INT,
                "TINT": C.x,
                "T_eff": C.y,
                })
        add_efective_temp[0] = add_2
        add['effective temp'] = add_efective_temp
        if NS > 1 and B[6] != 0:
            add_2.update({
                "NR": C.NBT,
                "NP": C.INT,
                "TINT": C.x,
                "T_eff": C.y,
                })
            add['effective temp'][1] = add_2
            if NS > 2 and B[12] != 0:
                add_2.update({
                    "NR": C.NBT,
                    "NP": C.INT,
                    "TINT": C.x,
                    "T_eff": C.y,
                    })
                add['effective temp'][2] = add_2
            if NS > 3 and B[18] != 0:
                add_2.update({
                    "NR": C.NBT,
                    "NP": C.INT,
                    "TINT": C.x,
                    "T_eff": C.y,
                    })
                add['effective temp'][3] = add_2
        out.update(add)
    return out
