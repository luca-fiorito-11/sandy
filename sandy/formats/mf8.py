# -*- coding: utf-8 -*-
"""
This module contains only two public functions:

    * `read`
    * `write`

Function `read` reads a MF8/MT section from a string and produces a content object with a dictionary-like 
structure.
The content object can be accessed using most of the keywords specified in the ENDF6 manual for this specific 
MF section.

Function `write` writes a content object for a MF8/MT section into a string.
MAT, MF, MT and line numbers are also added (each line ends with a `\n`).
"""

import numpy as np
import pdb

from sandy.formats.records import *
from sandy.formats.utils import Section

__author__ = "Luca Fiorito"
__all__ = ["read", "write"]

def read(text):
    """Read MT section for MF8
    
    Parameters
    ----------
    text : `str`
        one string containing the whole section
    
    Returns
    -------
    `sandy.formats.utils.Section`
        MF8 content sructured as a collection of dictionaries 
    """
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    if MT in (454, 459):
        return _read_fy(text)
    elif MT == 457:
        return _read_rdd(text)



def write(sec):
    """Write MT section for MF8
    
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



def _read_fy(text):
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    out = {"MAT" : MAT, "MF" : MF, "MT" : MT}
    i = 0
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2, "E" : {}})
    for j in range(C.L1):
        L, i = read_list(str_list, i)
        FY = { int(ZAFP*10+FPS) : {"ZAFP" : int(ZAFP), "FPS" : int(FPS), "YI" : YI, "DYI" : DYI} for ZAFP,FPS,YI,DYI in  zip(*[iter(L.B)]*4)}
        out["E"].update({ L.C1 : { "FY" : FY } })
        if j > 0:
            out["E"][L.C1].update({ "I" : L.L1 })
    return Section(out)



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



def _read_rdd(text):
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    out = {"MAT" : MAT, "MF" : MF, "MT" : MT}
    i = 0
    C, i = read_cont(str_list, i)
    out["ZA"] = ZA = C.C1       # Designation of the original (radioactive) nuclide (ZA = Z*1000 + A)
    out["AWR"] = AWR = C.C2     # Ratio of the LIS state nuclide mass to that of neutron
    out["LIS"] = LIS = C.L1     # State of the original nuclide (LIS=0, ground state, LIS=1, first excited state, etc.)
    out["LISO"] = LISO = C.L2   # Isomeric state number for the original nuclide (LISO=0, ground state; LISO=1, first isomeric state; etc.)
    out["NST"] = NST = C.N1     # Nucleus stability flag (NST=0, radioactive; NST=1, stable)
    NSP = C.N2     # Total number of radiation types (STYP) for which spectral information is given (NSP may be zero)
    L, i = read_list(str_list, i)
    out["HL"] = HL = L.C1       # half-life of the original nuclide (seconds)
    out["DHL"] = DHL = L.C2     # uncertainty on half-life
    out["E"] = E = L.B[::2]     # list of average decay energies (eV) for different radiation types.
    out["DE"] = DE = L.B[1::2]  # list of uncertainties on average decay energy (eV) for different radiation types.
    out["LAMBDA"] = np.log(2)/HL if HL else 0
    L, i = read_list(str_list, i)
    out["SPI"] = SPI = L.C1 # Spin of the nuclide in its LIS state
    out["PAR"] = PAR = L.C2 # Parity of the nuclide in its LIS state
    NDK = L.N2              # Total number of decay modes given
    dk = {}
    # Update list of decay modes when nuclide is radioactive
    for idk,(RTYP,RFS,Q,DQ,BR,DBR) in enumerate(zip(*[iter(L.B)]*6)):
        RTYP = str(RTYP).replace(".", "") # Decay Mode (Multiple particle decay is also allowed using combination of RTYP variables)
        decay = {
                 "RFS" : RFS,   # Isomeric state flag for daughter nuclide
                 "Q" : Q,       # Total decay energy (eV) available in the corresponding decay process
                 "DQ" : DQ,     # Uncertainty on the total decay heat available
                 "BR" : BR,     # Branching ratio
                 "DBR" : DBR    # Uncertainty on branching ratio
                }
        dk[RTYP] = decay
    if NDK > 0:
        out["DK"] = dk
    spectra = {}
    for ist in range(NSP):
        L, i = read_list(str_list, i)
        STYP = int(L.C2)  # Decay spectrum type
        spectra[STYP] = {}
        spectra[STYP]["LCON"] = LCON = L.L1 # Continuum spectrum flag (0=no continuous spectrum, 1=only continuous spectrum, 2=both discrete and continuum spectra)
        NER = L.N2   # Total number of tabulated discrete energies for a given spectral type (STYP)
        spectra[STYP]["FD"] = FD = L.B[0]   # Discrete spectrum normalization factor
        spectra[STYP]["DFD"] = L.B[1]
        spectra[STYP]["ERAV"] = L.B[2]      # Average decay energy of radiation produced
        spectra[STYP]["DERAV"] = L.B[3]
        spectra[STYP]["FC"] = L.B[4]        # Continuum spectrum normalization factor
        spectra[STYP]["DFC"] = L.B[5]
        if LCON != 1:
            discrete_spectrum = {}
            for ier in range(NER):
                discr = {}
                L, i = read_list(str_list, i)
                ER = L.C1                 # Discrete energy (eV) of radiation produced
                discr["DER"] = DER = L.C2 # Uncertainty on discrete energy
                NT = len(L.B)             # Number of entries given for each discrete energy (ER)
                discr["RTYP"] = L.B[0]    # Decay mode
                discr["TYPE"] = L.B[1]    # Type of transition for beta and electron capture
                discr["RI"] = L.B[2]      # intensity of discrete radiation produced (relative units).
                discr["DRI"] = L.B[3]     # uncertainty on the intensity of discrete radiation produced
                discr["RIS"] = L.B[4]     # Internal pair formation coefficient (STYP=2.0 positron intensity, STYP=0.0 otherwise)
                discr["DRIS"] = L.B[5]    # Uncertainty on internal pair formation coefficient 
                if NT == 12:
                    discr["RICC"] = L.B[6]   # Total internal conversion coefficient (STYP=0.0 only)
                    discr["DRICC"] = L.B[7]  # Uncertainty on RICC1
                    discr["RICK"] = L.B[8]   # K-shell internal conversion coefficient (STYP=0.0 only)
                    discr["DRICK"] = L.B[9]  # Uncertainty on RICC1
                    discr["RICL"] = L.B[10]  # L-shell internal conversion coefficient (STYP=0.0 only)
                    discr["DRICL"] = L.B[11] # Uncertainty on RICL1
                discrete_spectrum[ER] = discr
            if discrete_spectrum:
                spectra[STYP]["ER"] = discrete_spectrum
        if LCON != 0:
            spectra[STYP]["CONT"] = {}
            cont = {}
            T, i = read_tab1(str_list, i)
            cont["RTYP"] = T.C1        # Decay mode 
            cont["LCOV"] = LCOV = T.L2 # Flag indicating whether covariance data are given (0=no, 1=yes).
            cont["NBT"] = T.NBT
            cont["INT"] = T.INT
            cont["E"] = T.x,
            cont["RP"] = T.y           # Normalized spectrum of the continuum component of the radiation (1/eV)
            spectra[STYP]["CONT"][RTYP] = cont
    if spectra:
        out["SPECTRA"] = spectra
    return Section(out)