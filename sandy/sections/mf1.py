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
    if mt == 452:
        out = _read_nubar(tape, mat)
    elif mt == 455:
        out = _read_dnubar(tape, mat)
    elif mt == 456:
        out = _read_pnubar(tape, mat)
    elif mt in allowed_mt:
        raise sandy.Error(f"'MT={mt}' not yet implemented")
    else:
        raise ValueError(f"'MT={mt}' not allowed")
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
