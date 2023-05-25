import logging
import re

import sandy

__author__ = "Luca Fiorito"
__all__ = [
        "read_mf1",
        "write_mf1",
        ]

allowed_mt = (
        451,
        452,
        455,
        456,
        458,
        )
mf = 1


def read_mf1(tape, mat, mt):
    r"""
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

    Examples
    --------
    Since the outputs are very large and they are much, I am only going to
    check only some information for the test:

    **mt = 451** :
    >>> tape = sandy.get_endf6_file("endfb_71", 'xs', 922350)
    >>> test = sandy.read_mf1(tape, 9228, 451)
    >>> test['SECTIONS'][::5]
    [(1, 451, 934, 7),
    (1, 460, 253745, 0),
    (3, 4, 106, 7),
    (3, 20, 24, 6),
    (3, 52, 77, 6),
    (3, 57, 34, 7),
    (3, 62, 29, 7),
    (3, 67, 37, 7),
    (3, 72, 33, 7),
    (3, 77, 28, 7),
    (3, 82, 25, 7),
    (3, 87, 21, 7),
    (3, 102, 106, 6),
    (4, 53, 178, 6),
    (4, 58, 172, 7),
    (4, 63, 53, 7),
    (4, 68, 69, 7),
    (4, 73, 72, 7),
    (4, 78, 65, 7),
    (4, 83, 56, 7),
    (4, 88, 52, 7),
    (6, 16, 699, 6),
    (12, 18, 5, 6),
    (14, 4, 1, 6),
    (15, 18, 54, 6),
    (33, 2, 113418, 7),
    (33, 102, 23060, 7)]

    **mt = 452** :
    >>> test = sandy.read_mf1(tape, 9228, 452)
    >>> test['E']
    array([1.00e-05, 2.53e-02, 5.00e-02, 1.00e+01, 1.00e+02, 1.00e+03,
        5.50e+03, 7.75e+03, 1.00e+04, 1.50e+04, 2.00e+04, 3.00e+04,
        4.00e+04, 5.00e+04, 6.00e+04, 7.00e+04, 8.00e+04, 9.00e+04,
        1.00e+05, 1.20e+05, 1.30e+05, 1.40e+05, 1.50e+05, 1.70e+05,
        2.00e+05, 2.50e+05, 3.00e+05, 3.50e+05, 4.00e+05, 5.00e+05,
        6.00e+05, 7.00e+05, 8.00e+05, 9.00e+05, 1.00e+06, 1.20e+06,
        1.40e+06, 1.60e+06, 1.80e+06, 2.00e+06, 2.20e+06, 2.40e+06,
        2.60e+06, 2.80e+06, 3.00e+06, 3.50e+06, 4.00e+06, 4.50e+06,
        5.00e+06, 5.50e+06, 6.00e+06, 6.50e+06, 7.00e+06, 7.50e+06,
        8.00e+06, 8.50e+06, 9.00e+06, 9.50e+06, 1.00e+07, 1.05e+07,
        1.10e+07, 1.15e+07, 1.20e+07, 1.25e+07, 1.30e+07, 1.35e+07,
        1.40e+07, 1.45e+07, 1.50e+07, 1.55e+07, 1.60e+07, 1.65e+07,
        1.70e+07, 1.75e+07, 1.80e+07, 1.85e+07, 1.90e+07, 1.95e+07,
        2.00e+07])

    **mt = 455** :
    >>> test = sandy.read_mf1(tape, 9228, 455)
    >>> test['LAMBDA']
    [0.013336, 0.032739, 0.12078, 0.30278, 0.84949, 2.853]
    >>> test['NU']
    array([0.01585, 0.01585, 0.0167 , 0.0167 , 0.009  , 0.009  ])

    **mt = 456** :
    >>> test = sandy.read_mf1(tape, 9228, 456)
    >>> test['NU']
    array([2.42085 , 2.42085 , 2.42085 , 2.42085 , 2.417948, 2.417933,
        2.417857, 2.417818, 2.41778 , 2.414463, 2.412632, 2.409341,
        2.407132, 2.406774, 2.408063, 2.410846, 2.414175, 2.417648,
        2.421063, 2.428448, 2.431798, 2.434909, 2.43778 , 2.444676,
        2.451972, 2.455226, 2.45777 , 2.459884, 2.461767, 2.466965,
        2.472104, 2.481385, 2.491777, 2.502975, 2.516006, 2.540238,
        2.565402, 2.590224, 2.613985, 2.636654, 2.658912, 2.681288,
        2.703809, 2.727271, 2.7515  , 2.811124, 2.87605 , 2.95168 ,
        3.031212, 3.119151, 3.210984, 3.306917, 3.396401, 3.467412,
        3.535635, 3.608747, 3.680833, 3.752163, 3.821917, 3.890283,
        3.957945, 4.02513 , 4.092022, 4.160664, 4.232148, 4.305359,
        4.379344, 4.453294, 4.52846 , 4.605076, 4.681103, 4.754749,
        4.824436, 4.890934, 4.955851, 5.019246, 5.081174, 5.141692,
        5.200845])

    **mt = 458** :
    >>> test = sandy.read_mf1(tape, 9228, 458)
    >>> test['POLYNOMIALS'][1]
    {'EFR': -0.266,
    'DEFR': 0.0266,
    'ENP': 0.3004,
    'DENP': 0.03004,
    'END': 0.0,
    'DEND': 0.0,
    'EGP': 0.0777,
    'DEGP': 0.00777,
    'EGD': -0.075,
    'DEGD': 0.0075,
    'EB': -0.075,
    'DEB': 0.0075,
    'ENU': -0.1,
    'DENU': 0.01,
    'ER': -0.0379,
    'DER': 0.00379,
    'ET': -0.1379,
    'DET': 0.01379}

    >>> tape = sandy.get_endf6_file("endfb_71", 'nfpy', 922350)
    >>> test = read_mf1(tape, 9228, 451)
    >>> test['SECTIONS']
    [(1, 451, 17, 2), (8, 454, 2501, 2), (8, 459, 2501, 2)]

    >>> tape = sandy.get_endf6_file("endfb_71", 'decay', 922350)
    >>> test = sandy.read_mf1(tape, 3515, 451)
    >>> print("\n".join(test['DESCRIPTION']))
     92-U -235  BNL        EVAL-NOV05 Conversion from ENSDF           
     /ENSDF/                                               20111222   
    ----ENDF/B-VII.1      Material 3515                               
    -----RADIOACTIVE DECAY DATA                                       
    ------ENDF-6 FORMAT                                               
    *********************** Begin Description ***********************
    **         ENDF/B-VII.1 RADIOACTIVE DECAY DATA FILE            **
    **         Produced at the NNDC from the ENSDF database        **
    **               Translated into ENDF format by:               **
    **    T.D. Johnson, E.A. McCutchan and A.A. Sonzogni, 2011     **
    *****************************************************************
    ENSDF evaluation authors: E. BROWNE
    Parent Excitation Energy: 0
    Parent Spin & Parity: 7/2-
    Parent half-life: 703.8E+6 Y 5
    Decay Mode: A
    ************************ Energy  Balance ************************
    Mean Gamma Energy:      1.486E2 +- 1.440E0 keV
    Mean X-Ray+511 Energy:  1.553E1 +- 7.609E-1 keV
    Mean CE+Auger Energy:   4.170E1 +- 1.313E0 keV
    Mean B- Energy:         0.000E0 +- 0.000E0 keV
    Mean B+ Energy:         0.000E0 +- 0.000E0 keV
    Mean Neutrino Energy:   0.000E0 +- 0.000E0 keV
    Mean Neutron Energy:    0.000E0 +- 0.000E0 keV
    Mean Proton Energy:     0.000E0 +- 0.000E0 keV
    Mean Alpha Energy:      4.339E3 +- 1.648E2 keV
    Mean Recoil Energy:     7.386E1 +- 2.806E0 keV
    Sum Mean Energies:      4.619E3 +- 1.649E2 keV
    Q effective:            4.679E3 keV
    Missing Energy:         5.951E1 keV
    Deviation:              1.272E0 %
    ************************ End Description ************************
    """
    if mt == 451:
        out = _read_intro(tape, mat)
    elif mt == 452:
        out = _read_nubar(tape, mat)
    elif mt == 455:
        out = _read_dnubar(tape, mat)
    elif mt == 456:
        out = _read_pnubar(tape, mat)
    elif mt == 458:
        out = _read_fission_energy(tape, mat)
    elif mt in allowed_mt:
        raise ValueError("'MF={mf}/MT={mt}' not yet implemented")
    else:
        raise ValueError("'MF={mf}/MT={mt}' not allowed")
    return out


def _read_nubar(tape, mat):
    mt = 452
    df = tape._get_section_df(mat, mf, mt)
    out = {
            "MAT": mat,
            "MF": mf,
            "MT": mt,
            }
    i = 0
    C, i = sandy.read_cont(df, i)
    LNU = C.L2
    add = {
            "ZA": C.C1,
            "AWR": C.C2,
            "LNU": LNU,
            }
    out.update(add)
    if LNU == 1:
        L, i = sandy.read_list(df, i)
        add = {
                "C": L.B,
                }
        out.update(add)
    elif LNU == 2:
        T, i = sandy.read_tab1(df, i)
        add = {
                "NBT": T.NBT,
                "INT": T.INT,
                "E": T.x,
                "NU": T.y,
                }
        out.update(add)
    else:
        raise ValueError(f"'LNU' cannot be '{LNU}'")
    return out


def _read_pnubar(tape, mat):
    mt = 456
    df = tape._get_section_df(mat, mf, mt)
    out = {
            "MAT": mat,
            "MF": mf,
            "MT": mt,
            }
    i = 0
    C, i = sandy.read_cont(df, i)
    LNU = C.L2
    add = {
            "ZA": C.C1,
            "AWR": C.C2,
            "LNU": LNU,
            }
    out.update(add)
    if LNU == 1:
        L, i = sandy.read_list(df, i)
        add = {
                "NU": L.B,
                }
        out.update(add)
    elif LNU == 2:
        T, i = sandy.read_tab1(df, i)
        add = {
                "NBT": T.NBT,
                "INT": T.INT,
                "E": T.x,
                "NU": T.y,
                }
        out.update(add)
    else:
        raise ValueError(f"'LNU' cannot be '{LNU}'")
    return out


def _read_dnubar(tape, mat):
    mt = 455
    df = tape._get_section_df(mat, mf, mt)
    out = {
            "MAT": mat,
            "MF": mf,
            "MT": mt,
            }
    i = 0
    C, i = sandy.read_cont(df, i)
    LDG = C.L1
    LNU = C.L2
    add = {
            "ZA": C.C1,
            "AWR": C.C2,
            "LDG": LDG,
            "LNU": LNU,
            }
    out.update(add)
    if LDG == 0 and LNU == 2:
        L, i = sandy.read_list(df, i)
        add = {
                "LAMBDA": L.B,
                }
        out.update(add)
        T, i = sandy.read_tab1(df, i)
        add = {
                "NBT": T.NBT,
                "INT": T.INT,
                "E": T.x,
                "NU": T.y,
                }
        out.update(add)
    elif LDG == 1 and LNU == 2:
        T2, i = sandy.read_tab2(df, i)
        NE = T2.NZ
        add = {
                "ENBT": T2.NBT,
                "EINT": T2.INT,
                }
        out.update(add)
        add = {}
        for j in range(NE):
            L, i = sandy.read_list(df, i)
            E = L.C2
            LAMBDA = L.B[::2]
            ALPHA = L.B[1::2]
            add[E] = {
                    "LAMBDA": LAMBDA,
                    "ALPHA": ALPHA,
                    }
        out["EGROUPS"] = add
        T, i = sandy.read_tab1(df, i)
        add = {
                "NBT": T.NBT,
                "INT": T.INT,
                "E": T.x,
                "NU": T.y,
                }
        out.update(add)
    elif LDG == 0 and LNU == 1:
        logging.warning(f"""'(LDG, LNU) = ({LDG}, {LNU})' is not validated.
                        Please, report any possible bug/error.""")
        L, i = sandy.read_list(df, i)
        add = {
                "LAMBDA": L.B,
                }
        out.update(add)
        T, i = sandy.read_tab1(df, i)
        add = {
                "NBT": T.NBT,
                "INT": T.INT,
                "E": T.x,
                "NU": T.y,
                }
        out.update(add)
    elif LDG == 1 and LNU == 1:
        logging.warning(f"""'(LDG, LNU) = ({LDG}, {LNU})' is not validated.
                        Please, report any possible bug/error.""")
        T2, i = sandy.read_tab2(df, i)
        NE = T2.NZ
        add = {
                "ENBT": T2.NBT,
                "EINT": T2.INT,
                }
        out.update(add)
        add = {}
        for j in range(NE):
            L, i = sandy.read_list(df, i)
            E = L.C2
            LAMBDA = L.B[::2]
            ALPHA = L.B[1::2]
            add[E] = {
                    "LAMBDA": LAMBDA,
                    "ALPHA": ALPHA,
                    }
        out["EGROUPS"] = add
        T, i = sandy.read_tab1(df, i)
        add = {
                "NBT": T.NBT,
                "INT": T.INT,
                "E": T.x,
                "NU": T.y,
                }
        out.update(add)
    else:
        raise ValueError(f"'(LDG, LNU)' cannot be '({LDG}, {LNU})'")
    return out


def _read_intro(tape, mat):
    mt = 451
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
            "LRP": C.L1,   # Flag indicating whether resolved and/or unresolved resonance parameters are given in File 2
            "LFI": C.L2,   # Flag indicating whether this material fissions
            "NLIB": C.N1,  # Library identifier (e.g. NLIB= 0 for ENDF/B)
            "MOD": C.N2,   # Modification number for this material
            }
    out.update(add)
    C, i = sandy.read_cont(df, i)
    add = {
            "ELIS": C.C1,  # Excitation energy of the target nucleus relative to 0.0 for the ground state.
            "STA": C.C2,   # Target stability flag
            "LIS": C.L1,   # State number of the target nucleus
            "LISO": C.L2,  # Isomeric state number
            "NFOR": C.N2,  # Library format. NFOR=6 for ENDF-6
            }
    out.update(add)
    C, i = sandy.read_cont(df, i)
    add = {
            "AWI": C.C1,   # Mass of the projectile in neutron mass units
            "EMAX": C.C2,  # Upper limit of the energy range for evaluation
            "LREL": C.L1,  # Library release number; for example, LREL=2 for the ENDF/B-VI.2 library
            "NSUB": C.N1,  # Sub-library number
            "NVER": C.N2,  # Library version number; for example, NVER=7 for version ENDF/B-VII
            }
    out.update(add)
    C, i = sandy.read_cont(df, i)
    NWD = C.N1
    NXC = C.N2
    add = {
            "TEMP": C.C1,  # Target temperature (Kelvin) for data that have been generated by Doppler broadening
            "LDRV": C.L1,  # Special derived material flag that distinguishes between different evaluations with the same material keys
#            "NWD": NWD,    # Number of records with descriptive text for this material
#            "NXC": NXC,    # Number of records in the directory for this material
            }
    out.update(add)
    descr = []
    for j in range(NWD):
        T, i = sandy.read_text(df, i)
        descr.append(T[0])
    add = {
        "DESCRIPTION": descr,
        }
    out.update(add)
    # add = {
    #         "ZSYMAM": T.HL[:11],   # Character representation of the material
    #         "ALAB": T.HL[11:22],   # Mnemonic for the originating laboratory(s)
    #         "EDATE": T.HL[22:32],  # Date of evaluation
    #         "AUTH": T.HL[33:],     # Author(s) name(s)
    #         }
    # out.update(add)
    # T, i = sandy.read_text(df, i)
    # add = {
    #         "REF": T.HL[1:22],      # Primary reference for the evaluation
    #         "DDATE": T.HL[22:32],   # Original distribution date
    #         "RDATE": T.HL[33:43],   # Date and number of the last revision to this evaluation
    #         "ENDATE": T.HL[55:63],  # Author(s) name(s)
    #         }
    # out.update(add)
    # H1, i = sandy.read_text(df, i)
    # H2, i = sandy.read_text(df, i)
    # H3, i = sandy.read_text(df, i)
    # add = {
    #         "HSUB": "\n".join([H1.HL, H2.HL, H3.HL])
    #         }
    # out.update(add)
    # lines = []
    # for j in range(NWD - 5):
    #     T, i = sandy.read_text(df, i)
    #     lines.append(T.HL)
    # add = "\n".join(lines)
    # out.update({
    #     "DESCRIPTION": add,
    #     })
    # try:
    #     sections = _get_sections(df.iloc[-NXC:])
    # except Exception as e:
    #     msg = f"reported sections in MAT{mat}/MF1/MT451 are not consistent"
    #     logging.warning(msg)
    #     logging.warning(f"captured error: '{e}'")
    # else:
    #     out.update({
    #         "SECTIONS": sections,
    #         })
    sections = []
    for j in range(NXC):
        T, i = sandy.read_text(df, i)
        s = tuple(map(int, re.findall(".{11}" , T[0])[2:]))
        sections.append(s)
    out.update({
        "SECTIONS": sections,
        })
    return out


def _get_sections(df):
    sections = []
    for ipos, row in df.iterrows():
        MAT = int(row.L1)
        MF = int(row.L2)
        MT = int(row.N1)
        MOD = int(row.N2)
        add = MAT, MF, MT, MOD
        sections.append(add)
    return sections
    

def _read_fission_energy(tape, mat):
    mt = 458
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
            }
    out.update(add)
    L, i = sandy.read_list(df, i)
    NPLY = L.L2
    add = {
            "NPLY": NPLY,   # Order of the polynomial expansion of the energy-components
            "N2": L.N2,
            }
    out.update(add)
    poly = {}
    LIST = L.B[:]
    for pol_order in range(NPLY + 1):
        B = LIST[:18]
        LIST = LIST[18:]
        add = {
            "EFR": B[0],    # Kinetic energy of the fission products (following prompt neutron emission from the fission fragments)
            "DEFR": B[1],
            "ENP": B[2],    # Kinetic energy of the prompt fission neutrons
            "DENP": B[3],
            "END": B[4],    # Kinetic energy of the delayed fission neutrons
            "DEND": B[5],
            "EGP": B[6],    # Total energy released by the emission of prompt g rays
            "DEGP": B[7],
            "EGD": B[8],    # Total energy released by the emission of delayed g rays
            "DEGD": B[9],
            "EB": B[10],    # Total energy released by delayed betas
            "DEB": B[11],
            "ENU": B[12],   # Energy carried away by neutrinos
            "DENU": B[13],
            "ER": B[14],    # Total energy less the energy of the neutrinos (ET - ENU); equal to the pseudo-Q-value in File 3 for MT=18
            "DER": B[15],
            "ET": B[16],    # Sum of all the partial energies
            "DET": B[17],
            }
        poly.update({
            pol_order: add,
            })
    out.update({
        "POLYNOMIALS": poly,
        })
    return out


def write_mf1(sec):
    """
    Given the content of a MF1 section as nested dictionaries, write it
    to string.

    Returns
    -------
    `str`
        Multiline string reproducing the content of a ENDF-6 section.

    Notes
    -----
    .. note:: The end-of-line records MAT, MF, MT and line number are added at
              the end of each line.

    .. important:: The string does not end with a newline symbol `\n`.

    Examples
    --------
    Since the outputs are very large and they are much, I am only going to
    check only some information for the test:

    **mt = 452** :
    >>> tape = sandy.get_endf6_file("endfb_71", 'xs', 922350)
    >>> sec = sandy.read_mf1(tape, 9228, 452)
    >>> text = sandy.write_mf1(sec)
    >>> print(text[:1000])
    92235.0000 233.024800          0          2          0          09228 1452    1
    0.00000000 0.00000000          0          0          1         799228 1452    2
            79          2                                            9228 1452    3
    1.000000-5 2.43670000 2.530000-2 2.43670000 5.000000-2 2.436700009228 1452    4
    10.0000000 2.43670000 100.000000 2.43380000 1000.00000 2.433800009228 1452    5
    5500.00000 2.43380000 7750.00000 2.43380000 10000.0000 2.433800009228 1452    6
    15000.0000 2.43056800 20000.0000 2.42882200 30000.0000 2.425701009228 1452    7
    40000.0000 2.42366200 50000.0000 2.42347400 60000.0000 2.424763009228 1452    8
    70000.0000 2.42754600 80000.0000 2.43087500 90000.0000 2.434348009228 1452    9
    100000.000 2.43776300 120000.000 2.44514800 130000.000 2.448498009228 1452   10
    140000.000 2.45160900 150000.000 2.45448000 170000.000 2.461376009228 1452   11
    200000.000 2.46867200 250000.000 2.47192600 300000.000 2.474470009228 1452   12
    350000.000 2.47658400 40000

    **mt = 455** :
    >>> sec = sandy.read_mf1(tape, 9228, 455)
    >>> text = sandy.write_mf1(sec)
    >>> print(text[:1000])
    92235.0000 233.024800          0          2          0          09228 1455    1
    0.00000000 0.00000000          0          0          6          09228 1455    2
    1.333600-2 3.273900-2 1.207800-1 3.027800-1 8.494900-1 2.853000009228 1455    3
    0.00000000 0.00000000          0          0          1          69228 1455    4
            6          2                                            9228 1455    5
    1.000000-5 1.585000-2 2.530000-2 1.585000-2 50000.0000 1.670000-29228 1455    6
    4000000.00 1.670000-2 7000000.00 9.000000-3 20000000.0 9.000000-39228 1455    7

    **mt = 456** :
    >>> sec = sandy.read_mf1(tape, 9228, 456)
    >>> text = sandy.write_mf1(sec)
    >>> print(text[:1000])
    92235.0000 233.024800          0          2          0          09228 1456    1
    0.00000000 0.00000000          0          0          1         799228 1456    2
            79          2                                            9228 1456    3
    1.000000-5 2.42085000 2.530000-2 2.42085000 5.000000-2 2.420850009228 1456    4
    10.0000000 2.42085000 100.000000 2.41794800 1000.00000 2.417933009228 1456    5
    5500.00000 2.41785700 7750.00000 2.41781800 10000.0000 2.417780009228 1456    6
    15000.0000 2.41446300 20000.0000 2.41263200 30000.0000 2.409341009228 1456    7
    40000.0000 2.40713200 50000.0000 2.40677400 60000.0000 2.408063009228 1456    8
    70000.0000 2.41084600 80000.0000 2.41417500 90000.0000 2.417648009228 1456    9
    100000.000 2.42106300 120000.000 2.42844800 130000.000 2.431798009228 1456   10
    140000.000 2.43490900 150000.000 2.43778000 170000.000 2.444676009228 1456   11
    200000.000 2.45197200 250000.000 2.45522600 300000.000 2.457770009228 1456   12
    350000.000 2.45988400 40000
    """
    mt = sec["MT"]
    if mt == 451:
        out = _write_intro(sec)
    elif mt == 452:
        out = _write_nubar(sec)
    elif mt == 455:
        out = _write_dnubar(sec)
    elif mt == 456:
        out = _write_pnubar(sec)
    elif mt in allowed_mt:
        raise ValueError(f"'MT={mt}' not yet implemented")
    else:
        raise ValueError(f"'MT={mt}' not allowed")
    return out


def _write_nubar(sec):
    mat = sec["MAT"]
    mt = 452
    LNU = sec["LNU"]
    lines = sandy.write_cont(
            sec["ZA"],
            sec["AWR"],
            0,
            LNU,
            0,
            0,
            )
    if LNU == 1:
        lines += sandy.write_list(
                0,
                0,
                0,
                0,
                0,
                sec["C"],
                )
    elif LNU == 2:
        lines += sandy.write_tab1(
                0,
                0,
                0,
                0,
                sec["NBT"],
                sec["INT"],
                sec["E"],
                sec["NU"],
                )
    else:
        raise ValueError(f"'LNU' cannot be '{LNU}'")
    return "\n".join(sandy.write_eol(lines, mat, mf, mt))


def _write_pnubar(sec):
    mat = sec["MAT"]
    mt = 456
    LNU = sec["LNU"]
    lines = sandy.write_cont(
            sec["ZA"],
            sec["AWR"],
            0,
            LNU,
            0,
            0,
            )
    if LNU == 1:
        lines += sandy.write_list(
                0,
                0,
                0,
                0,
                0,
                sec["NU"],
                )
    elif LNU == 2:
        lines += sandy.write_tab1(
                0,
                0,
                0,
                0,
                sec["NBT"],
                sec["INT"],
                sec["E"],
                sec["NU"],
                )
    else:
        raise ValueError(f"'LNU' cannot be '{LNU}'")
    return "\n".join(sandy.write_eol(lines, mat, mf, mt))


def _write_dnubar(sec):
    mat = sec["MAT"]
    mt = 455
    LDG = sec["LDG"]
    LNU = sec["LNU"]
    lines = sandy.write_cont(
            sec["ZA"],
            sec["AWR"],
            LDG,
            LNU,
            0,
            0,
            )
    if LDG == 0 and LNU == 2:
        lines += sandy.write_list(
                0,
                0,
                0,
                0,
                0,
                sec["LAMBDA"],
                )
        lines += sandy.write_tab1(
                0,
                0,
                0,
                0,
                sec["NBT"],
                sec["INT"],
                sec["E"],
                sec["NU"],
                )
    elif LDG == 1 and LNU == 2:
        lines += sandy.write_tab2(
                0,
                0,
                0,
                0,
                len(sec["EGROUPS"]),
                sec["ENBT"],
                sec["EINT"],
                )
        for e, v in sec["EGROUPS"]:
            LAMBDA = v["LAMBDA"]
            ALPHA = v["ALPHA"]
            lines += sandy.write_list(
                    0,
                    e,
                    0,
                    0,
                    0,
                    [item for pair in zip(LAMBDA, ALPHA) for item in pair],
                    )
        lines += sandy.write_tab1(
                0,
                0,
                0,
                0,
                sec["NBT"],
                sec["INT"],
                sec["E"],
                sec["NU"],
                )
    elif LDG == 0 and LNU == 1:
        lines += sandy.write_list(
                0,
                0,
                0,
                0,
                0,
                sec["LAMBDA"],
                )
        lines += sandy.write_tab1(
                0,
                0,
                0,
                0,
                sec["NBT"],
                sec["INT"],
                sec["E"],
                sec["NU"],
                )
    elif LDG == 1 and LNU == 1:
        lines += sandy.write_tab2(
                0,
                0,
                0,
                0,
                len(sec["EGROUPS"]),
                sec["ENBT"],
                sec["EINT"],
                )
        for e, v in sec["EGROUPS"]:
            LAMBDA = v["LAMBDA"]
            ALPHA = v["ALPHA"]
            lines += sandy.write_list(
                    0,
                    e,
                    0,
                    0,
                    0,
                    [item for pair in zip(LAMBDA, ALPHA) for item in pair],
                    )
        lines += sandy.write_tab1(
                0,
                0,
                0,
                0,
                sec["NBT"],
                sec["INT"],
                sec["E"],
                sec["NU"],
                )
    else:
        raise ValueError(f"'(LDG, LNU)' cannot be '({LDG}, {LNU})'")
    return "\n".join(sandy.write_eol(lines, mat, mf, mt))


def _write_intro(sec):
    mat = sec["MAT"]
    mt = 451
    lines = sandy.write_cont(
            sec["ZA"],
            sec["AWR"],
            sec["LRP"],
            sec["LFI"],
            sec["NLIB"],
            sec["MOD"],
            )
    lines += sandy.write_cont(
            sec["ELIS"],
            sec["STA"],
            sec["LIS"],
            sec["LISO"],
            0,
            sec["NFOR"],
            )
    lines += sandy.write_cont(
            sec["AWI"],
            sec["EMAX"],
            sec["LREL"],
            0,
            sec["NSUB"],
            sec["NVER"],
            )
    NWD = len(sec["DESCRIPTION"])
    NXC = len(sec["SECTIONS"])
    lines += sandy.write_cont(
            sec["TEMP"],
            0,
            sec["LDRV"],
            0,
            NWD,
            NXC,
            )
    for t in sec["DESCRIPTION"]:
        lines += sandy.write_text(t)
    for MF, MT, NL, MOD in sec["SECTIONS"]:
        t = " "*22 + f"{MF:>11d}{MT:>11d}{NL:>11d}{MOD:>11d}"
        lines += sandy.write_text(t)
    return "\n".join(sandy.write_eol(lines, mat, mf, mt))
