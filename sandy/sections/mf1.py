import pdb
import logging

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
    """
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
        raise sandy.Error("'MF={mf}/MT={mt}' not yet implemented")
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
            "NWD": NWD,    # Number of records with descriptive text for this material
            "NXC": NXC,    # Number of records in the directory for this material
            }
    out.update(add)
    T, i = sandy.read_text(df, i)
    add = {
            "ZSYMAM": T.HL[:11],   # Character representation of the material
            "ALAB": T.HL[11:22],   # Mnemonic for the originating laboratory(s)
            "EDATE": T.HL[22:32],  # Date of evaluation
            "AUTH": T.HL[33:],     # Author(s) name(s)
            }
    out.update(add)
    T, i = sandy.read_text(df, i)
    add = {
            "REF": T.HL[1:22],      # Primary reference for the evaluation
            "DDATE": T.HL[22:32],   # Original distribution date
            "RDATE": T.HL[33:43],   # Date and number of the last revision to this evaluation
            "ENDATE": T.HL[55:63],  # Author(s) name(s)
            }
    out.update(add)
    H1, i = sandy.read_text(df, i)
    H2, i = sandy.read_text(df, i)
    H3, i = sandy.read_text(df, i)
    add = {
            "HSUB": "\n".join([H1.HL, H2.HL, H3.HL])
            }
    out.update(add)
    lines = []
    for j in range(NWD):
        T, i = sandy.read_text(df, i)
        lines.append(T.HL)
    add = "\n".join(lines)
    out.update({
        "DESCRIPTIION": add,
        })
    try:
        sections = _get_sections(df.iloc[-NXC:])
    except Exception as e:
        msg = f"reported sections in MAT{mat}/MF1/MT451 are not cosistent"
        logging.warning(msg)
        logging.warning(f"captured error: '{e}'")
    else:
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
    """
    mt = sec["MT"]
    if mt == 452:
        out = _write_nubar(sec)
    elif mt == 455:
        out = _write_dnubar(sec)
    elif mt == 456:
        out = _write_pnubar(sec)
    elif mt in allowed_mt:
        raise sandy.Error(f"'MT={mt}' not yet implemented")
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
                len(sec["EGRPOUPS"]),
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
                len(sec["EGRPOUPS"]),
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
