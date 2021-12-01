# -*- coding: utf-8 -*-
"""
This module contains only two public functions:

    * `read_mf8`
    * `write_mf8`

Function `read` reads a MF8/MT section from a string and produces a content
object with a dictionary-like structure.
The content object can be accessed using most of the keywords specified in
the ENDF6 manual for this specific MF section.

Function `write_mf8` writes a content object for a MF8/MT section into a
string.
MAT, MF, MT and line numbers are also added (each line ends with a `\n`).
"""
import math

import sandy

__author__ = "Luca Fiorito"
__all__ = [
        "read_mf8",
        "write_mf8",
        ]

mf = 8


def read_mf8(tape, mat, mt):
    """
    Parse MAT/MF=8/MT section from `sandy.Endf6` object and return
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
    if mt in (454, 459):
        out = _read_fy(tape, mat, mt)
    elif mt == 457:
        out = _read_rdd(tape, mat)
    else:
        out = _read_nucl_prod(tape, mat, mt)
#    else:
#        raise ValueError(f"'MF={mf}/MT={mt}' not yet implemented")
    return out


def write_mf8(sec):
    """
    Write MT section for MF8

    Parameters
    ----------
    sec : 'dict'
        Multiline string reproducing the content of a ENDF-6 section.

    Returns
    -------
    `str`
    """
    if sec["MT"] in (454, 459):
        return _write_fy(sec)
    elif sec["MT"] == 457:
        return _write_rdd(sec)


def _read_nucl_prod(tape, mat, mt):
    """
    Parse MAT/MF=8/MT section for radioactive nuclide production from
    `sandy.Endf6` object and return structured content in nested dcitionaries.

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
    NS = C.N1
    add = {
            "ZAM": int(C.C1*10),
            "AWR": C.C2,
            "LIS": C.L1,   # State number (including ground and all levels) of the target (ZA)
            "LISO": C.L2,  # Isomeric state number of the target
            "NS": NS,      # Total number of states (LFS) of the radioactive reaction product for which decay data are given
            "NO": C.N2,    # Flag denoting where the decay information is to be given for an important radioactive end product.
            }
    out.update(add)
    products = {}
    for j in range(NS):
        L, i = sandy.read_list(df, i)
        LFS = L.L2
        ND = int(L.NPL / 6)
        LIST = L.B[:]
        decay_modes = []
        for k in range(ND):
            B = LIST[:6]
            LIST = LIST[6:]
            add = {
                "HL": B[0],    # half-life of the nuclide ZAP in seconds
                "RTYP": B[1],  # Mode of decay using the same definitions specified in MT=457
                "ZAN": B[2],   # Z and mass identifier of the next nuclide produced along the chain
                "BR": B[3],    # Branching ratio for the production of that particular ZAN and level
                "END": B[4],   # Endpoint energy of the particle or quantum emitted
                "CT": B[5],    # Chain terminator that gives minimal information about the formation and decay of ZAN
                }
            decay_modes.append(add)
        ZAP = int(L.C1)
        add = {
            "ZAP": ZAP,
            "ELFS": L.C2,
            "LMF": L.L1,
            "LFS": LFS,
            }
        if ND > 0:
            add["DECAY"] = decay_modes
        products[ZAP] = add
    out.update({
        "PRODUCTS": products,
        })
    return out


def _read_fy(tape, mat, mt):
    """
    Parse MAT/MF=8/MT section for fission yields from `sandy.Endf6` object
    and return structured content in nested dcitionaries.

    Parameters
    ----------
    tape : `sandy.Endf6`
        endf6 object containing requested section
    mat : `int`
        MAT number

    Returns
    -------
    `dict`
        Content of the ENDF-6 tape structured as nested `dict`.

    Notes
    -----
    .. note:: Fission yields are only contained in sections with `mt=454` (IFY)
              or `mt=459` (CFY).

    Examples
    --------
    >>> nfpy = sandy.get_endf6_file("jeff_33", "nfpy", 922350)
    >>> IFY = sandy.sections.mf8.read_mf8(nfpy, 9228, 454)
    >>> IFY["E"][0.0253]['ZAP'][10010]
    {'FY': 1.711e-05, 'DFY': 2.9483e-06}

    >>> nfpy = sandy.get_endf6_file("jeff_33", "nfpy", 922350)
    >>> IFY = sandy.sections.mf8.read_mf8(nfpy, 9228, 459)
    >>> IFY["E"][0.0253]['ZAP'][10010]
    {'FY': 1.711e-05, 'DFY': 1.8479e-06}
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
            "ZAM": int(C.C1*10),
            "AWR": C.C2,
            }
    out.update(add)
    nranges = C.L1
    eranges = {}
    for j in range(nranges):
        L, i = sandy.read_list(df, i)
        data = {}
        for zafp, fps, fy, dfy in zip(*[iter(L.B)]*4):
            zap = int(zafp*10 + fps)
            data[zap] = {
                    "FY": fy,
                    "DFY": dfy,
                    }
        section = {"ZAP": data}
        section["INTERP"] = L.L1
        eranges[L.C1] = section
    out["E"] = eranges
    return out


def _read_rdd(tape, mat):
    """
    Parse MAT/MF=8/MT=457 section from `sandy.Endf6` object and return
    structured content in nested dcitionaries.

    Parameters
    ----------
    tape : `sandy.Endf6`
        endf6 object containing requested section
    mat : `int`
        MAT number

    Returns
    -------
    `dict`
        Content of the ENDF-6 tape structured as nested `dict`.

    Examples
    --------
    >>> decay = sandy.get_endf6_file("jeff_33", "decay", 922350)
    >>> rdd = sandy.sections.mf8.read_mf8(decay, 3542, 457)
    >>> rdd['SPECTRA'][0]['ER'][19595.0]
    {'DER': 4.0,
    'RTYP': 4.0,
    'TYPE': 0.0,
    'RI': 0.00011,
    'DRI': 2e-05,
    'RIS': 0.0,
    'DRIS': 0.0,
    'RICC': 9670.0,
    'DRICC': 967.0,
    'RICK': 0.0,
    'DRICK': 0.0,
    'RICL': 0.0,
    'DRICL': 0.0}

    >>> decay = sandy.get_endf6_file("jeff_33", "decay", 922350)
    >>> rdd = sandy.sections.mf8.read_mf8(decay, 3542, 457)
    >>> rdd["DLAMBDA"]
    2.2171470275223715e-20
    """
    mt = 457
    df = tape._get_section_df(mat, mf, mt)
    out = {
            "MAT": mat,
            "MF": mf,
            "MT": mt,
            }
    i = 0
    C, i = sandy.read_cont(df, i)
    add = {
            # Designation of the original (radioactive) nuclide (Z*1000 + A)
            "ZA": C.C1,
            # Ratio of the LIS state nuclide mass to that of neutron
            "AWR": C.C2,
            # State of the original nuclide (LIS=0, ground state,
            # LIS=1, first excited state, etc.)
            "LIS": C.L1,
            # Isomeric state number for the original nuclide
            # (LISO=0, ground state; LISO=1, first isomeric state; etc.)
            "LISO": C.L2,
            # Nucleus stability flag (NST=0, radioactive; NST=1, stable)
            "NST": C.N1,
            }
    out.update(add)
    # Total number of radiation types (STYP) for which spectral information is
    # given (NSP may be zero)
    NSP = C.N2
    L, i = sandy.read_list(df, i)
    add = {
            # half-life of the original nuclide (seconds)
            "HL": L.C1,
            # uncertainty on half-life
            "DHL": L.C2,
            # list of average decay energies (eV) for different radiation types
            "E": L.B[::2],
            # list of uncertainties on average decay energy (eV) for different
            # radiation types
            "DE": L.B[1::2],
            # decay constant in 1/s, 0 if stable
            "LAMBDA": math.log(2.0) / L.C1 if L.C1 else 0,
            # uncertainty on decay constant
            "DLAMBDA": math.log(2.0)*L.C2 / L.C1**2 if L.C1 else 0
            }

    out.update(add)
    L, i = sandy.read_list(df, i)
    add = {
            # Spin of the nuclide in its LIS state
            "SPI": L.C1,
            # Parity of the nuclide in its LIS state
            "PAR": L.C2,
            }
    out.update(add)
    # Total number of decay modes given
    NDK = L.N2
    ###########################
    # READ DECAY MODES
    ###########################
    if NDK > 0:
        dk = {}
        # Update list of decay modes when nuclide is radioactive
        for idk, data in enumerate(zip(*[iter(L.B)]*6)):
            # Decay Mode (Multiple particle decay is also allowed using
            # combination of RTYP variables)
            RTYP = str(data[0]).replace(".", "")
            # Isomeric state flag for daughter nuclide
            RFS = data[1]
            residual_state = int(RFS)
            key = f"{RTYP}x{residual_state}"
            decay = {
                "RTYP": RTYP,
                "RFS": RFS,
                # Total decay energy (eV) available in the corresponding
                # decay process
                "Q": data[2],
                # Uncertainty on the total decay heat available
                "DQ": data[3],
                # Branching ratio
                "BR": data[4],
                # Uncertainty on branching ratio
                "DBR": data[5],
                }
            dk[key] = decay
        out["DK"] = dk
    ###########################
    # READ SPECTRA
    ###########################
    if NSP > 0:
        spectra = {}
        # Update list of spectra
        for ist in range(NSP):
            L, i = sandy.read_list(df, i)
            # Decay spectrum type
            STYP = int(L.C2)
            spectra[STYP] = {}
            # Continuum spectrum flag (0=no continuous spectrum, 1=only
            # continuous spectrum, 2=both discrete and continuum spectra)
            spectra[STYP]["LCON"] = LCON = L.L1
            # Total number of tabulated discrete energies for a given spectral
            # type (STYP)
            NER = L.N2
            # Discrete spectrum normalization factor
            spectra[STYP]["FD"] = L.B[0]
            spectra[STYP]["DFD"] = L.B[1]
            # Average decay energy of radiation produced
            spectra[STYP]["ERAV"] = L.B[2]
            spectra[STYP]["DERAV"] = L.B[3]
            # Continuum spectrum normalization factor
            spectra[STYP]["FC"] = L.B[4]
            spectra[STYP]["DFC"] = L.B[5]
            if LCON != 1:
                discrete_spectrum = {}
                for ier in range(NER):
                    discr = {}
                    L, i = sandy.read_list(df, i)
                    # Discrete energy (eV) of radiation produced
                    ER = L.C1
                    # Uncertainty on discrete energy
                    discr["DER"] = L.C2
                    # Number of entries given for each discrete energy (ER)
                    NT = len(L.B)
                    if NT > 0:
                        # Decay mode
                        discr["RTYP"] = L.B[0]
                    if NT > 1:
                        # Type of transition for beta and electron capture
                        discr["TYPE"] = L.B[1]
                    if NT > 2:
                        # intensity of discrete radiation produced
                        # (relative units)
                        discr["RI"] = L.B[2]
                    if NT > 3:
                        # uncertainty on the intensity of discrete radiation
                        # produced
                        discr["DRI"] = L.B[3]
                    if NT > 4:
                        # Internal pair formation coefficient
                        # (STYP=2.0 positron intensity, STYP=0.0 otherwise)
                        discr["RIS"] = L.B[4]
                    if NT > 5:
                        # Uncertainty on internal pair formation coefficient
                        discr["DRIS"] = L.B[5]
                    if NT > 6:
                        # Total internal conversion coefficient
                        # (STYP=0.0 only)
                        discr["RICC"] = L.B[6]
                    if NT > 7:
                        # Uncertainty on RICC1
                        discr["DRICC"] = L.B[7]
                    if NT > 8:
                        # K-shell internal conversion coefficient
                        # (STYP=0.0 only)
                        discr["RICK"] = L.B[8]
                    if NT > 9:
                        # Uncertainty on RICC1
                        discr["DRICK"] = L.B[9]
                    if NT > 10:
                        # L-shell internal conversion coefficient
                        # (STYP=0.0 only)
                        discr["RICL"] = L.B[10]
                    if NT > 11:
                        # Uncertainty on RICL1
                        discr["DRICL"] = L.B[11]
                    discrete_spectrum[ER] = discr
                if discrete_spectrum:
                    spectra[STYP]["ER"] = discrete_spectrum
            if LCON != 0:
                spectra[STYP]["CONT"] = {}
                cont = {}
                T, i = sandy.read_tab1(df, i)
                # Decay mode
                cont["RTYP"] = T.C1
                # Flag indicating whether covariance data are given
                # (0=no, 1=yes)
                cont["LCOV"] = T.L2
                cont["NBT"] = T.NBT
                cont["INT"] = T.INT
                cont["E"] = T.x,
                # Normalized spectrum of the continuum component of the
                # radiation (1/eV)
                cont["RP"] = T.y
                spectra[STYP]["CONT"][RTYP] = cont
        if spectra:
            out["SPECTRA"] = spectra
    return out


def _write_fy(sec):
    """
    Given the content of MF8 for fission yields as nested dictionaries,
    write it to string.

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
    Independent fission yield:
    >>> nfpy = sandy.get_endf6_file("jeff_33", "nfpy", 922350)
    >>> sec = sandy.sections.mf8.read_mf8(nfpy, 9228, 454)
    >>> text = _write_fy(sec)
    >>> print(text[:1000])
     922350.000 233.025000          3          0          0          09228 8454    1
     2.530000-2 0.00000000          2          0       3932        9839228 8454    2
     1001.00000 0.00000000 1.711000-5 2.948300-6 1002.00000 0.000000009228 8454    3
     8.400000-6 2.438900-6 1003.00000 0.00000000 1.080000-4 5.863500-69228 8454    4
     2003.00000 0.00000000 0.00000000 0.00000000 2004.00000 0.000000009228 8454    5
     1.700000-3 8.134700-5 2006.00000 0.00000000 2.668000-5 4.878400-69228 8454    6
     3006.00000 0.00000000 0.00000000 0.00000000 3008.00000 0.000000009228 8454    7
     7.292000-7 2.368000-7 3009.00000 0.00000000 4.071000-7 1.049100-79228 8454    8
     4008.00000 0.00000000 7.292000-7 1.650700-7 4009.00000 0.000000009228 8454    9
     4.071000-7 1.653500-7 4010.00000 0.00000000 5.201000-6 9.085300-79228 8454   10
     4012.00000 0.00000000 1.261000-7 4.498200-8 5009.00000 0.000000009228 8454   11
     4.071000-7 1.101800-7 5010.00000 0.00000000 5.201000-6 8.994800-79228 8454   12
     5012.00000 0.00000000 1.261

    Cumulative fission yield:
    >>> nfpy = sandy.get_endf6_file("jeff_33", "nfpy", 922350)
    >>> sec = sandy.sections.mf8.read_mf8(nfpy, 9228, 459)
    >>> text = _write_fy(sec)
    >>> print(text[:1000])
     922350.000 233.025000          3          0          0          09228 8459    1
     2.530000-2 0.00000000          2          0       3932        9839228 8459    2
     1001.00000 0.00000000 1.711000-5 1.847900-6 1002.00000 0.000000009228 8459    3
     8.400000-6 1.503600-6 1003.00000 0.00000000 1.080000-4 3.996000-69228 8459    4
     2003.00000 0.00000000 1.080000-4 3.996000-6 2004.00000 0.000000009228 8459    5
     1.702100-3 4.930000-5 2006.00000 0.00000000 2.668000-5 1.840900-69228 8459    6
     3006.00000 0.00000000 2.668000-5 1.840900-6 3008.00000 0.000000009228 8459    7
     7.292000-7 1.181300-7 3009.00000 0.00000000 4.071000-7 2.849700-89228 8459    8
     4008.00000 0.00000000 1.341800-6 1.181300-7 4009.00000 0.000000009228 8459    9
     6.126800-7 2.849700-8 4010.00000 0.00000000 5.201000-6 2.548500-79228 8459   10
     4012.00000 0.00000000 1.261000-7 3.013800-8 5009.00000 0.000000009228 8459   11
     4.071000-7 2.849700-8 5010.00000 0.00000000 5.201000-6 2.548500-79228 8459   12
     5012.00000 0.00000000 2.522
    """
    LE = len(sec["E"])
    lines = sandy.write_cont(sec["ZAM"], sec["AWR"], LE, 0, 0, 0)
    for E, FY_information in sec['E'].items():
        interp = FY_information['INTERP']
        energy_data = FY_information['ZAP']
        NFP = len(energy_data)
        B = []
        for zap, fy_dat in energy_data.items():
            z, a, m = sandy.zam.expand_zam(zap)
            zaps, fps = sandy.zam.get_za(z, a, m, method='none')
            fy = fy_dat['FY']
            dfy = fy_dat['DFY']
            B.extend([zaps, fps, fy, dfy])
        lines += sandy.write_list(
                    E,
                    0,
                    interp,
                    0,
                    NFP,
                    B,
                    )
    return "\n".join(sandy.write_eol(lines, sec["MAT"], 8, sec["MT"]))


def _write_rdd(sec):
    """
    Given the content of MF8 for radioactive decay as nested dictionaries,
    write it to string.

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
    Stable nuclide:
    >>> decay = sandy.get_endf6_file("jeff_33", "decay", 551340)
    >>> sec = sandy.sections.mf8.read_mf8(decay,1803,457)
    >>> text = _write_rdd(sec)
    >>> print(text[:1000])
     55134.0000 132.756000          0          0          0          51803 8457    1
     65146100.0 44179.7000          0          0          6          01803 8457    2
     163850.000 403.088000 1554660.00 895.993000 0.00000000 0.000000001803 8457    3
     4.00000000 1.00000000          0          0         12          21803 8457    4
     20.0000000 0.00000000 1233300.00 800.000000 3.000000-6 1.000000-61803 8457    5
     10.0000000 0.00000000 2058980.00 330.000000 9.999970-1 1.000000-61803 8457    6
     0.00000000 0.00000000          0          0          6         121803 8457    7
     1.000000-2 0.00000000 1554380.00 895.965000 0.00000000 0.000000001803 8457    8
     242760.000 50.0000000          0          0         12          01803 8457    9
     1.00000000 0.00000000 2.409990-2 3.099990-3 0.00000000 0.000000001803 8457   10
     8.700000-2 3.000000-3 7.220000-2 3.289850-3 1.200000-2 2.551310-31803 8457   11
     326585.000 14.0000000          0          0         12          01803 8457   12
     1.00000000 0.00000000 1.709

    Unstable nuclide:
    >>> decay = sandy.get_endf6_file("jeff_33", "decay", 922350)
    >>> sec = sandy.sections.mf8.read_mf8(decay,3542,457)
    >>> text = _write_rdd(sec)
    >>> print(text[:1000])
     92235.0000 233.025000          0          0          0          63542 8457    1
     2.22102+16 1.57788+13          0          0          6          03542 8457    2
     50671.7000 4291.63000 163616.000 1708.01000 4464600.00 163255.0003542 8457    3
     3.50000000-1.00000000          0          0         12          23542 8457    4
     40.0000000 0.00000000 4678700.00 700.000000 1.00000000 1.000000-43542 8457    5
     60.0000000 0.00000000  176400000 4600000.00 7.20000-11 2.10000-113542 8457    6
     0.00000000 0.00000000          2          0          6         503542 8457    7
     5.710000-1 3.000000-3 148550.000 804.564000 5.82048-10 1.79455-103542 8457    8
     19595.0000 4.00000000          0          0         12          03542 8457    9
     4.00000000 0.00000000 1.100000-4 2.000000-5 0.00000000 0.000000003542 8457   10
     9670.00000 967.000000 0.00000000 0.00000000 0.00000000 0.000000003542 8457   11
     31580.0000 10.0000000          0          0         12          03542 8457   12
     4.00000000 0.00000000 3.000
    """
    if sec['NST'] == 1:
        lines = sandy.write_cont(sec["ZA"],
                                 sec["AWR"],
                                 sec['LIS'],
                                 sec['LISO'],
                                 1,
                                 0,
                                 )
        lines += sandy.write_list(0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  [0]*6,
                                  )
        lines += sandy.write_list(sec['SPI'],
                                  sec['PAR'],
                                  0,
                                  0,
                                  0,
                                  [0]*6,
                                  )
    elif sec['NST'] == 0:
        lines = sandy.write_cont(sec["ZA"],
                                 sec["AWR"],
                                 sec['LIS'],
                                 sec['LISO'],
                                 0,
                                 len(sec['SPECTRA']),
                                 )
        size_E = len(sec['E'])
        size_DE = len(sec['DE'])
        add = [0]*(size_E + size_DE)
        add[::2] = sec['E']
        add[1::2] = sec['DE']
        lines += sandy.write_list(
                            sec['HL'],
                            sec['DHL'],
                            0,
                            0,
                            0,
                            add,
                            )
        NDK = len(sec['DK'])
        add = []
        for key, decay_data in sec['DK'].items():
            add.extend([int(decay_data['RTYP']),
                        decay_data['RFS'],
                        decay_data['Q'],
                        decay_data['DQ'],
                        decay_data['BR'],
                        decay_data['DBR'],
                        ])
        lines += sandy.write_list(
            sec['SPI'],
            sec['PAR'],
            0,
            0,
            NDK,
            add,
            )
        for STYP, spectra in sec['SPECTRA'].items():
            if spectra['LCON'] != 1 and 'ER' in list(sec['SPECTRA'][STYP].keys()):
                add = [spectra['FD'],
                       spectra['DFD'],
                       spectra['ERAV'],
                       spectra['DERAV'],
                       spectra['FC'],
                       spectra['DFC'],
                       ]
                lines += sandy.write_list(
                    0,
                    STYP,
                    spectra['LCON'],
                    0,
                    len(spectra['ER']),
                    add
                    )
                for ER, discr in spectra['ER'].items():
                    lines += sandy.write_list(
                        ER,
                        discr['DER'],
                        0,
                        0,
                        0,
                        list(discr.values())[1:],
                        )
                if spectra['LCON'] != 0 and 'CONT' in list(sec['SPECTRA'][STYP].keys()):
                    for RTYP, cont in spectra['CONT'].items():
                        lines += sandy.write_tab1(
                            cont['RTYP'],
                            0.0,
                            0,
                            cont['LCOV'],
                            cont['NBT'],
                            cont['INT'],
                            cont['E'][0],
                            cont["RP"],
                            )
        return "\n".join(sandy.write_eol(lines, sec["MAT"], 8, sec["MT"]))
