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
        # "write_mf8",
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
    return out


def write(sec):
    """
    Write MT section for MF8

    Parameters
    ----------
    sec : `sandy.utils.Section`
        dictionary with MT section for MF8

    Returns
    -------
    `str`
    """
    if sec["MT"] in (454, 459):
        return _write_fy(sec)
    elif sec["MT"] == 457:
        return _write_rdd(sec)


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
            "LAMBDA": math.log(2.0)/L.C1 if L.C1 else 0
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
            decay = {
                    # Isomeric state flag for daughter nuclide
                    "RFS": data[1],
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
            dk[RTYP] = decay
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
    LE = len(sec["E"])
    text = write_cont(sec["ZA"], sec["AWR"], LE, 0, 0, 0)
    for i,(e,esec) in enumerate(sorted(sec["E"].items())):
        tab = [ksec[j] for k,ksec in sorted(esec["FY"].items()) for j in ("ZAFP","FPS","YI","DYI")]
        NFP = len(esec["FY"])
        I = LE-1 if i == 0 else esec["I"]
        text += write_list(e, 0, I, 0, NFP, tab)
    TextOut = []; iline = 1
    for line in text:
        if iline > 99999:
            iline = 1
        TextOut.append("{:<66}{:4}{:2}{:3}{:5}\n".format(line, sec["MAT"], sec["MF"], sec["MT"], iline))
        iline += 1
    return "".join(TextOut)


def _write_rdd(*args, **kwargs):
    pass
