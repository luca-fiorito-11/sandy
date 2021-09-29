"""
This module contains only two public functions:

    * `read_mf6`
    * `write_mf6`

Function `read` reads a MF6/MT section from a string and produces a content
object with a dictionary-like structure.
The content object can be accessed using most of the keywords specified in
the ENDF6 manual for this specific MF section.

Function `write_mf6` writes a content object for a MF6/MT section into a
string.
MAT, MF, MT and line numbers are also added (each line ends with a `\n`).
"""
__author__ = "Aitor Bengoechea"
__all__ = [
        "read_mf6",
        "write_mf6",
        ]

import sandy
import logging
from sandy import zam


def read_mf6(tape, mat, mt):
    """
    Parse MAT/MF=6/MT section from `sandy.Endf6` object and return
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
    Since the outputs are very large, I am only going to check only some
    information for the test:

    LAW 1:

    >>> tape = sandy.get_endf6_file("endfb_71", 'xs', 70140)
    >>> test = read_mf6(tape, 725, 5)
    >>> test["NK"][10010]["EGROUPS"][1e-05]
     {'ND': 0,
     'NA': 1,
     'NW': 6,
     'NEP': 2,
     'Ep': [0.0, 1e-05],
     'b': [100000.0, 0.0, 0.0, 0.0]}

    LAW 2:

    >>> import pprint
    >>> tape = sandy.get_endf6_file("endfb_71", 'xs', 10010)
    >>> test = read_mf6(tape, 125, 102)
    >>> test["NK"][10020]['AWP'] = round (test["NK"][10020]['AWP'], 5)
    >>> pprint.pprint(test["NK"][10020])
    {'AWP': 1.99626,
     'E': array([1.e-05, 2.e+07]),
     'LAW': 4,
     'NP': [2],
     'NR': [2],
     'Y': array([1., 1.])}

    LAW 6:

    >>> tape = sandy.get_endf6_file("endfb_71", 'xs', 10020)
    >>> test = read_mf6(tape, 128, 16)
    >>> test["NK"][10]['APSX'] = round (test["NK"][10]['APSX'],5)
    >>> test["NK"][10]
    {'AWP': 1.0,
     'LAW': 6,
     'NR': [2],
     'NP': [2],
     'E': array([3.339002e+06, 1.500000e+08]),
     'Y': array([2., 2.]),
     'APSX': 2.99862,
     'NPSX': 3}

    LAW 7:

    >>> tape = sandy.get_endf6_file("endfb_71", 'xs', 40090)
    >>> test = read_mf6(tape, 425, 16)
    >>> test["NK"][10]["EGROUPS"][1748830.0]['COSGROUPS'][-1.0]
    {'NRP': [15],
     'NEP': [2],
     'EP_INT': array([  1092.99,   1093.  ,   3278.9 ,   7650.8 ,  12023.  ,  20766.  ,
             29510.  ,  55741.  ,  71043.  ,  81973.  ,  90716.  ,  95088.  ,
             99460.  , 101650.  , 101651.  ]),
     'E_p': array([0.00000e+00, 1.16614e-06, 1.58588e-06, 1.54624e-06, 7.09710e-07,
            1.79581e-07, 2.86776e-08, 0.00000e+00]),
     'E_distr': array([7.40674e-07, 1.46654e-06, 1.61094e-06, 1.07195e-06, 4.02172e-07,
            9.52648e-08, 4.69275e-09])}
    """
    mf = 6
    df = tape._get_section_df(mat, mf, mt)
    out = {"MAT": mat, "MF": mf, "MT": mt}
    i = 0
    C, i = sandy.read_cont(df, i)
    out.update({
                "ZA": C.C1,
                "AWR": C.C2,
                "LCT": C.L2,  # Reference system for secondary energy and angle
               })
    subsections = {}
    # Each subsection describes one reaction product. There can be more than
    # one subsection for a given particle, but the combination of the product
    # identifier with its final isomeric state create a unique identifier for
    # each reaction.
    for a in range(C.N1):
        T, i = sandy.read_tab1(df, i)
        LAW = T.L2  # Distintion between different distribution function
        ZAP = T.C1  # Product identifier
        LIP = T.L1  # Product isomeric state identifier
        ZAM = zam.za2zam(ZAP, meta=LIP)
        add = {
                "AWP": T.C2,  # Product mass in neutron units
                "LAW": T.L2,
                "NR": T.NBT,
                "NP": T.INT,
                "E": T.x,  # Neutron incident energy
                "Y": T.y,  # The product multiplicity
                }
        # LAW dependent structures:
        if LAW == 1:  # Continuum Energy-Angle Distributions
            L, i = sandy.read_tab2(df, i)
            NE = L.NBT[0]  # How many NE incident energies
            add.update({
                        "LANG": L.L1,  # Angular representation identificator
                        "LEP": L.L2,  # Interpolation for secondary energy
                        "ENR": L.NR,  # I put here ENR insted of NR to do not overwrite NR of the tab1 section.
                        "ENE": L.NBT,  # Number of different product energy
                        "EINT": L.INT,  # product energy diferent values
                        })
            add_e = {}
            # To repeat for all the NE incident energies
            for j in range(NE):
                T, i = sandy.read_list(df, i)
                E = T.C2  # Incident energy
                if int(T.L2) == 0:
                    Ep = T.B[::2]
                    b = T.B[1::2]
                else:
                    Ep = T.B[::T.L2+2]
                    b = T.B[::1]
                    del b[::T.L2+2]  # To delete from b the energy values
                add_2 = {
                         "ND": T.L1,  # Number of discrete energies
                         "NA": T.L2,  # Number of angular parameters
                         "NW": T.NPL,  # Total number of words
                         "NEP": T.N2,  # Secondary energy points in distribution
                         "Ep": Ep,  # The energy of the product emitted
                         "b": b,  # Coefficients for the angular representation
                         # the contents of the b depend on LANG
                       }
                add_e[E] = add_2
            add["EGROUPS"] = add_e

        elif LAW == 2:  # Discrete Two-Body Scattering
            L, i = sandy.read_tab2(df, i)
            NE = L.NBT[0]
            add.update({
                        "ENR": L.NZ,
                        "ENE": L.NBT,  # Number of different product energy
                        "EINT": L.INT,  # Product energy values
                        })
            add_e = {}
            # To repeat list records for all the incident energies
            for j in range(NE):
                T, i = sandy.read_list(df, i)
                E = T.C2  # Incident energy
                add_2 = {
                         "LANG": T.L1,  # Angular representation identificator
                         "NW": T.NPL,  # Number of parameters
                         "NL": T.N2,  # Highest coefficient identificator
                         "Al": T.B,  # Coeficcient of angular representation
                         }
                add_e[E] = add_2
            add["EGROUPS"] = add_e

        elif LAW == 5:  # Charged-Particle Elastic Scattering
            logging.warning(f"""'(LAW) = ({LAW})' is not validated.
                            Please report any posible error/bug.""")
            L, i = sandy.read_tab2(df, i)
            NE = L.NBT[0]  # How many NE incident energies
            LIDP = L.L1
            add.update({
                        "SPI": L.C1,
                        "LIDP": LIDP,
                        "ENR": L.NR,
                        "ENE": L.NBT,
                        "EINT": L.INT,
            })
            add_e = {}
            for j in range(NE):  # To repeat the LIST records for all the NE
                T, i = sandy.read_list(df, i)
                E = T.C2
                LTP = T.L1
                add_2 = {
                         "LTP": T.L1,
                         "NW": T.NPL,
                         "NL": T.N2,
                         }
                # We have to do the distintion between the different Ai organizations.
                if LTP == 1 and LIDP == 0 or LIDP == 1:
                    mark = T.N2+1
                    b = T.B[0:mark]
                    Ra = T.B[mark::2]
                    Ia = T.B[mark+1::2]
                    add_2["A"] = {
                                             "B": b,
                                             "Ra": Ra,
                                             "Ia": Ia,
                                             }
                elif LTP == 2:
                    add_2["A"] = {
                                             "C": T.B,
                                             }
                elif LTP > 2:
                    nu = T.B[::2]
                    p = T.B[1::2]
                    add_2["A"] = {
                                             "nu": nu,
                                             "p": p,
                                             }
                add_e[E] = add_2
            add["EGROUPS"] = add_e

        elif LAW == 6:  # N-Body Phase-Space Distributions
            T, i = sandy.read_cont(df, i)
            add.update({
                        "APSX": T.C1,  # Total mass(neutron uni) of N particles
                        "NPSX": T.N2,  # Number of particles distributed
                             })
        elif LAW == 7:  # Laboratory Angle-Energy Law
            L, i = sandy.read_tab2(df, i)
            NE = L.NBT[0]  # How many interpolation range we have
            # Interpolation parameters for incident energy E
            add.update({
                                "ENR": L.NZ,
                                "ENE": L.NBT,  # The incident energies number
                                "EINT": L.INT,  # Incident energy values
                             })
            add_e = {}
            # To repeat for all NE incident energies
            for j in range(NE):
                T, i = sandy.read_tab2(df, i)
                # Interpolation parameters for emission cosine
                E = T.C2  # Incident energy
                NMU = T.NBT[0]  # Number of possible emission cosines
                add_2 = {
                            "NRM": T.NZ,
                            "NMU": T.NBT,  # Number of different cosines
                            "NU_INT": T.INT,  # Value of emission cosines
                         }
                add_e[E] = add_2
                # To repeat for all the NMU emission cosines
                add_nu = {}
                for hz in range(NMU):
                    Z, i = sandy.read_tab1(df, i)
                    # Interpolation parameters for secondary energy E′
                    nu = Z.C2  # Value for emission cosine
                    E_p = Z.y[::2]
                    E_distr = Z.y[1::2]
                    add_3 = {
                                "NRP": Z.NBT,
                                "NEP": Z.INT,
                                "EP_INT": Z.x,  # Energies for E'
                                "E_p": E_p,  # Secondary energy value
                                "E_distr": E_distr,  # Distribution according to nu, E and E'
                                }
                    add_nu[nu] = add_3
                add_e[E]["COSGROUPS"] = add_nu
            add["EGROUPS"] = add_e
        subsections[ZAM] = add
    out["NK"] = subsections
    return out


def write_mf6(sec):
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
    As we have a law-dependent structure, I will develop a test for each law.
    LAW 1:
    >>> tape = sandy.get_endf6_file("endfb_71", 'xs', 70140)
    >>> sec = read_mf6(tape, 725, 5)
    >>> text = write_mf6(sec)
    >>> print(text[:1000])
     7014.00000 13.8828000          0          3         22          0 725 6  5    1
     1.00000000 1.00000000          0          1          1         19 725 6  5    2
             19          2                                             725 6  5    3
     1.000000-5 7.880151-1 20000000.0 7.880151-1 20000010.0 7.880151-1 725 6  5    4
     23000000.0 8.736206-1 27000000.0 1.04544100 30000000.0 1.09352800 725 6  5    5
     35000000.0 1.15156600 40000000.0 1.20729400 50000000.0 1.29613300 725 6  5    6
     60000000.0 1.36702400 70000000.0 1.42382400 80000000.0 1.51414400 725 6  5    7
     90000000.0 1.61200000  100000000 1.70471900  110000000 1.93624400 725 6  5    8
      120000000 2.02639000  130000000 2.09783300  140000000 2.17252600 725 6  5    9
      150000000 2.23843500                                             725 6  5   10
     0.00000000 0.00000000          2          1          1          1 725 6  5   11
             19          2                                             725 6  5   12
     0.00000000 1.000000-5

    LAW 2:
    >>> tape = sandy.get_endf6_file("endfb_71", 'xs', 10010)
    >>> sec = read_mf6(tape, 125, 102)
    >>> text = write_mf6(sec)
    >>> print(text[:1000])
     1001.00000 9.991673-1          0          2          2          0 125 6102    1
     0.00000000 2223300.00          0          2          1          2 125 6102    2
              2          2                                             125 6102    3
     1.000000-5 1.00000000 20000000.0 1.00000000                       125 6102    4
     0.00000000 0.00000000          0          0          1         96 125 6102    5
             96          2                                             125 6102    6
     0.00000000 1.000000-5          0          0          2          2 125 6102    7
    -6.160950-8-8.86237-13                                             125 6102    8
     0.00000000 2.000000-5          0          0          2          2 125 6102    9
    -8.712894-8-1.77244-12                                             125 6102   10
     0.00000000 5.000000-5          0          0          2          2 125 6102   11
    -1.377628-7-4.43100-12                                             125 6102   12
     0.00000000 1.000000-4

    LAW 6:
    >>> tape = sandy.get_endf6_file("endfb_71", 'xs', 10020)
    >>> sec = read_mf6(tape, 128, 16)
    >>> text = write_mf6(sec)
    >>> print (text)
     1002.00000 1.99680000          0          1          2          0 128 6 16    1
     1.00000000 1.00000000          0          6          1          2 128 6 16    2
              2          2                                             128 6 16    3
     3339002.00 2.00000000  150000000 2.00000000                       128 6 16    4
     2.99862000 0.00000000          0          0          0          3 128 6 16    5
     1001.00000 9.986200-1          0          6          1          2 128 6 16    6
              2          2                                             128 6 16    7
     3339002.00 1.00000000  150000000 1.00000000                       128 6 16    8
     2.99862000 0.00000000          0          0          0          3 128 6 16    9

    LAW 7:
    >>> tape = sandy.get_endf6_file("endfb_71", 'xs', 40090)
    >>> sec = read_mf6(tape, 425, 16)
    >>> text = write_mf6(sec)
    >>> print(text[:1000])
     4009.00000 8.93478000          0          1          2          0 425 6 16    1
     1.00000000 1.00000000          0          7          1         17 425 6 16    2
             17          1                                             425 6 16    3
     1748830.00 2.00000000 2700000.00 2.00000000 3250000.00 2.00000000 425 6 16    4
     3900000.00 2.00000000 4500000.00 2.00000000 5900000.00 2.00000000 425 6 16    5
     6400000.00 2.00000000 7050000.00 2.00000000 8400000.00 2.00000000 425 6 16    6
     10100000.0 2.00000000 11000000.0 2.00000000 12000000.0 2.00000000 425 6 16    7
     13000000.0 2.00000000 14200000.0 2.00000000 15400000.0 2.00000000 425 6 16    8
     17000000.0 2.00000000 20000000.0 2.00000000                       425 6 16    9
     0.00000000 0.00000000          0          0          1         18 425 6 16   10
             18          2                                             425 6 16   11
     0.00000000 1748830.00          0          0          1         21 425 6 16   12
             21          2

    """
    lines = sandy.write_cont(
                sec["ZA"],
                sec["AWR"],
                0,
                sec['LCT'],
                len(sec["NK"]),
                0,
                )
    for product, NK in sec["NK"].items():  # For the rest of the NK subsections
        # [MAT, 6, MT/ ZAP, AWP, LIP, LAW,NR,NP/Eint/yi(E)]TAB1
        ZAP, LIP = zam.zam2za(product)
        lines += sandy.write_tab1(
                    ZAP,
                    NK["AWP"],
                    LIP,
                    NK["LAW"],
                    NK["NR"],
                    NK["NP"],
                    NK["E"],
                    NK["Y"],
                    )
        # LAW dependent structures:
        if NK["LAW"] == 1:
            # [MAT, 6, MT/ 0.0, 0.0, LANG, LEP, NR, NE/ Eint]TAB2 for each NK
            lines += sandy.write_tab2(
                            0,
                            0,
                            NK["LANG"],
                            NK["LEP"],
                            NK["ENR"],
                            NK["ENE"],
                            NK["EINT"],
                    )
            # repeat the LIST structure for all the NE incident energies
            for key, NK_E in NK["EGROUPS"].items():  # NE incident energies
                # Because the dividision that i made between energies and b
                    size_E = len(NK_E["Ep"])
                    size_b = len(NK_E["b"])
                    add = [0]*(size_E + size_b)
                    if NK_E["NA"] == 0:
                        add[::2] = NK_E["Ep"]
                        add[1::2] = NK_E["b"]
                    else:
                        n = NK_E["NA"]+1
                        x = [NK_E["b"][i:i + n] for i in range(0, size_b, n)]
                        for i in range(size_E):
                            add.append(NK_E["Ep"][i])
                            for j in range(len(x[i])):
                                add.append(x[i][j])
                # [MAT, 6, MT/ 0.0, E1, ND, NA, NW, NEP/
                # E′1, b0(E1,E′1), b1(E1,E′1 ), -------- bNA(E1,E′1),
                # E′2, b0(E1,E′2), b1(E1,E′2 ), -------- bNA(E1,E′2),
                # --------------------------------------------
                # E′N EP , b0(E1,E′NEP), b1(E1,E′NEP), ---- bNA(E1,E′NEP)]LIST
                        lines += sandy.write_list(
                               0,
                               key,
                               NK_E["ND"],
                               NK_E["NA"],
                               NK_E["NEP"],
                               add,
                                )
        elif NK["LAW"] == 2:
            # [MAT, 6, MT/ 0.0, 0.0, 0, 0, NR, NE/ Eint]TAB2 for each NK
            lines += sandy.write_tab2(
                    0,
                    0,
                    0,
                    0,
                    NK["ENR"],
                    NK["ENE"],
                    NK["EINT"],
            )
            # [MAT, 6, MT/ 0.0,E1,LANG,0,NW,NL/Al(E)]LIST for each Ne
            for key, NK_E in NK["EGROUPS"].items():
                lines += sandy.write_list(
                    0,
                    key,
                    NK_E["LANG"],
                    0,
                    NK_E["NL"],
                    NK_E["Al"],
                    )
        elif NK["LAW"] == 5:
            LAW = NK["LAW"]
            logging.warning(f"""'(LAW) = ({LAW})' is not validated.
                            Please report any posible error/bug.""")
            LIDP = NK["LIDP"]
            lines += sandy.write_tab2(
                                        NK["SPI"],
                                        LIDP,
                                        NK["ENR"],
                                        NK["ENE"],
                                        NK["EINT"],
                                        )
            for key, NK_E in NK['EGROUPS'].items():
                LTP = NK_E["LTP"]
                if LTP == 1 and LIDP == 0 or LIDP == 1:
                    add = NK_E["A"]["B"]
                    Ra = NK_E["A"]["Ra"]
                    Ia = NK_E["A"]["Ia"]
                    add_2 = [0]*(len(Ra)+len(Ia))
                    add_2[::2] = Ra
                    add_2[1::2] = Ia
                    add.append(add_2)
                elif LTP == 2:
                    add = NK_E["A"]["C"]
                elif LTP > 2:
                    nu = NK_E["A"]["nu"]
                    p = NK_E["A"]["p"]
                    add = [0]*(len(nu) + len(p))
                    add[::2] = nu
                    add[1::2] = p
                lines += sandy.write_list(
                                            0,
                                            key,
                                            LTP,
                                            NK_E["NW"],
                                            NK_E["NL"],
                                            add,
                                        )
        elif NK["LAW"] == 6:
            lines += sandy.write_cont(
                NK["APSX"],
                0,
                0,
                0,
                0,
                NK["NPSX"],
                )
        elif NK["LAW"] == 7:
            # [MAT, 6, MT/ 0.0,0.0,0,0,NR,NE/Eint]TAB2
            lines += sandy.write_tab2(
                    0,
                    0,
                    0,
                    0,
                    NK["ENR"],
                    NK["ENE"],
                    NK["EINT"],
                    )
            # [MAT, 6, MT/ 0.0, E1, 0, 0, NRM, NMU/μint]TAB2 for all the NE
            for key, NK_E in NK["EGROUPS"].items():
                lines += sandy.write_tab2(
                    0,
                    key,
                    0,
                    0,
                    NK_E["NRM"],
                    NK_E["NMU"],
                    NK_E["NU_INT"],
                    )
                # [MAT, 6, MT/ 0.0, μ1, 0, 0, NRP, NEP/E′int/
                #     E′1,f(μ1,E1,E′1 ), E′2,f(μ1,E1,E′2 ), -----
                #     -------------------------
                # E′N EP ,f(μ1,E1,E′N EP )]TAB1 for all the emission cosines
                for key_nu, NK_E_NU in NK_E["COSGROUPS"].items():
                    add = [0]*(len(NK_E_NU["E_p"])+len(NK_E_NU["E_distr"]))
                    add[::2] = NK_E_NU["E_p"]
                    add[1::2] = NK_E_NU["E_distr"]
                    lines += sandy.write_tab1(
                        0,
                        key_nu,
                        0,
                        0,
                        NK_E_NU["NRP"],
                        NK_E_NU["NEP"],
                        NK_E_NU["EP_INT"],
                        add,
                        )
    return "\n".join(sandy.write_eol(lines, sec["MAT"], 6, sec["MT"]))
