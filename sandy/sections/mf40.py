# -*- coding: utf-8 -*-
"""
"""

import sandy

__author__ = "Jan Malec"
__all__ = [
        "read_mf40",
        "write_mf40"
        ]


def read_mf40(tape, mat, mt, mf=40):
    """
    Parse MAT/MF=8/MT section for MF40 covariance matrix from
    `sandy.Endf6` object and return structured content in nested dcitionaries.

    Parameters
    ----------
    tape : `sandy.Endf6`
        endf6 object containing requested section
    mat : `int`
        MAT number
    mtt : `int`
        MT number
    mf : `int`, optional, default is `40`
        MF number, it also allows parsing nubar covariances with MF=31

    Returns
    -------
    `dict`
        Content of the ENDF-6 tape structured as nested `dict`.

    Notes
    -----
    """
    df = tape._get_section_df(mat, mf, mt)
    out = {
            "MAT": mat,
            "MF": mf,
            "MT": mt,
            }
    i = 0
    C, i = sandy.read_cont(df, i)
    out["ZA"] = C.C1
    out["AWR"] = C.C2
    out["NS"] = C.N1
    nsub = C.N1
    subs = {}
    for j in range(nsub):
        C, i = sandy.read_cont(df, i)
        qm = C.C1
        qi = C.C2
        izap = C.L1
        lfs = C.L2
        nl = C.N2
        lfs_key = (izap, lfs)
        sub_list = []
        for _ in range(nl):
            C, i = sandy.read_cont(df, i)
            xmf1 = C.C1
            xlfs1 = C.C2
            mat1 = C.L1
            mt1 = C.L2
            nc = C.N1
            ni = C.N2
            sub = {
                "QM": qm,
                "QI": qi,
                "IZAP": izap,
                "LFS": lfs,
                "XMF1": xmf1,
                "XLFS1": xlfs1,
                "MAT1": mat1,
                "MT1": mt1,
            }
            # Read NC blocks
            ncdict = {}
            for k in range(nc):
                C, i = sandy.read_cont(df, i)
                lty = C.L2
                subsub = {"LTY": lty}
                L, i = sandy.read_list(df, i)
                if lty == 0:
                    subsub["E1"] = L.C1
                    subsub["E2"] = L.C2
                    subsub["CI"] = L.B[:L.N2]
                    subsub["XMTI"] = L.B[L.N2:]
                elif lty in (1, 2, 3):
                    subsub["E1"] = L.C1
                    subsub["E2"] = L.C2
                    subsub["MATS"] = L.L1
                    subsub["MTS"] = L.L2
                    subsub["XMFS"] = L.B[0]
                    subsub["XLFSS"] = L.B[1]
                    subsub["EI"] = L.B[2:2 + L.N2]
                    subsub["WEI"] = L.B[2 + L.N2:]
                else:
                    raise AssertionError
                ncdict[k] = subsub
            sub["NC"] = ncdict
            # Read NI blocks
            nidict = {}
            for k in range(ni):
                L, i = sandy.read_list(df, i)
                lb = L.L2
                subsub = {"LB": lb}
                if lb in [0, 1, 2, 3, 4]:
                    subsub["LT"] = L.L1
                    subsub["NT"] = L.NPL
                    subsub["NP"] = L.N2
                    if subsub["LT"] == 0:
                        subsub["EK"] = L.B[::2]
                        subsub["FK"] = L.B[1::2]
                    else:
                        nk = subsub["NP"] - subsub["LT"]
                        arrk = L.B[:nk * 2]
                        subsub["EK"] = arrk[::2]
                        subsub["FK"] = arrk[1::2]
                        arrl = L.B[nk * 2:]
                        subsub["EL"] = arrl[::2]
                        subsub["FL"] = arrl[1::2]
                elif lb == 5:
                    subsub["LS"] = L.L1
                    subsub["NT"] = L.NPL
                    subsub["NE"] = L.N2
                    subsub["EK"] = L.B[:L.N2]
                    subsub["FKK"] = L.B[L.N2:]
                elif lb == 6:
                    subsub["NT"] = L.NPL
                    subsub["NER"] = L.N2
                    subsub["NEC"] = (L.NPL - 1) // L.N2
                    ner = subsub["NER"]
                    nec = subsub["NEC"]
                    subsub["EK"] = L.B[:ner]
                    subsub["EL"] = L.B[ner:ner + nec]
                    subsub["FKL"] = L.B[ner + nec:]
                elif lb in (8, 9):
                    subsub["LT"] = L.L1
                    subsub["NT"] = L.NPL
                    subsub["NP"] = L.N2
                    subsub["EK"] = L.B[:L.N2]
                    subsub["FK"] = L.B[L.N2:]
                else:
                    raise AssertionError("LB not found")
                nidict[k] = subsub
            sub["NI"] = nidict
            sub_list.append(sub)
        subs[lfs_key] = sub_list
    if subs:
        out["SUB"] = subs
    return out

def write_mf40(sec):
    """
    Given the content of a MF40 section as nested dictionaries, write it
    to string.

    Parameters
    ----------
    sec : `dict`
        Content of the ENDF-6 tape structured as nested `dict`, typically
        obtained from `read_mf40`.

    Returns
    -------
    `str`
        Multiline string reproducing the content of an ENDF-6 section.

    Notes
    -----
    The string does not end with a newline symbol `\\n`.
    """
    lines = sandy.write_cont(
        sec["ZA"],
        sec["AWR"],
        0,
        0,
        sec.get("NS", len(sec.get("SUB", {}))),  # NS might not be present if SUB is empty
        0,
    )
    if "SUB" in sec:
        # SUB keys are (IZAP, LFS) and values are lists of subsections
        for (izap, lfs), sublist in sorted(sec["SUB"].items()):
            nl = len(sublist)
            # First CONT record for the (IZAP, LFS) group
            lines += sandy.write_cont(
                sublist[0]["QM"],
                sublist[0]["QI"],
                izap,
                lfs,
                0,
                nl
            )
            for sub in sublist:
                # Second CONT for each subsection
                lines += sandy.write_cont(
                    sub["XMF1"],
                    sub["XLFS1"],
                    sub["MAT1"],
                    sub["MT1"],
                    len(sub.get("NC", {})),
                    len(sub.get("NI", {})),
                )
                # NC subsections
                for k in sorted(sub.get("NC", {}).keys()):
                    nc_sub = sub["NC"][k]
                    lines += sandy.write_cont(
                        0, 0, 0, nc_sub["LTY"], 0, 0
                    )
                    if nc_sub["LTY"] == 0:
                        b = list(nc_sub["CI"]) + list(nc_sub["XMTI"])
                        nci = len(nc_sub["CI"])
                        lines += sandy.write_list(
                            nc_sub["E1"],
                            nc_sub["E2"],
                            0,
                            0,
                            nci,
                            b,
                        )
                    elif nc_sub["LTY"] in (1, 2, 3):
                        nei = len(nc_sub["EI"])
                        b = [nc_sub["XMFS"], nc_sub["XLFSS"]] + list(nc_sub["EI"]) + list(nc_sub["WEI"])
                        lines += sandy.write_list(
                            nc_sub["E1"],
                            nc_sub["E2"],
                            nc_sub["MATS"],
                            nc_sub["MTS"],
                            nei,
                            b,
                        )
                    else:
                        raise ValueError(f"Cannot write unsupported LTY={nc_sub['LTY']} in NC section")
                # NI subsections
                for k in sorted(sub.get("NI", {}).keys()):
                    ni_sub = sub["NI"][k]
                    lb = ni_sub["LB"]
                    b = []
                    l1 = 0
                    npl = ni_sub.get("NT", 0)
                    n2 = ni_sub.get("NP", 0)

                    if lb in [0, 1, 2, 3, 4]:
                        l1 = ni_sub.get("LT", 0)
                        if l1 == 0:
                            for ek, fk in zip(ni_sub.get("EK", []), ni_sub.get("FK", [])):
                                b.extend([ek, fk])
                        else:
                            for ek, fk in zip(ni_sub.get("EK", []), ni_sub.get("FK", [])):
                                b.extend([ek, fk])
                            for el, fl in zip(ni_sub.get("EL", []), ni_sub.get("FL", [])):
                                b.extend([el, fl])
                    elif lb == 5:
                        l1 = ni_sub.get("LS", 0)
                        n2 = ni_sub.get("NE", 0)
                        b = list(ni_sub.get("EK", [])) + list(ni_sub.get("FKK", []))
                    elif lb == 6:
                        l1 = 0
                        n2 = ni_sub.get("NER", 0)
                        b = list(ni_sub.get("EK", [])) + list(ni_sub.get("EL", [])) + list(ni_sub.get("FKL", []))
                    elif lb in [8, 9]:
                        l1 = ni_sub.get("LT", 0)
                        b = list(ni_sub.get("EK", [])) + list(ni_sub.get("FK", []))
                    else:
                        raise ValueError(f"Cannot write unsupported LB={lb} in NI section")

                    lines += sandy.write_list(
                        0,
                        0,
                        l1,
                        lb,
                        n2,
                        b
                    )
    return "\n".join(sandy.write_eol(lines, sec["MAT"], sec["MF"], sec["MT"]))
