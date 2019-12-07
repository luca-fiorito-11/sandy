import pdb

import sandy

__author__ = "Luca Fiorito"
__all__ = [
        "read_mf3",
        "write_mf3"
        ]

def read_mf3(tape, mat, mt):
    mf = 3
    df = tape._get_section_df(mat, mf, mt)
    out = {"MAT" : mat, "MF" : mf, "MT" : mt}
    i = 0
    C, i = sandy.read_cont(df, i)
    out.update({
            "ZA" : C.C1,
            "AWR" : C.C2,
            "PFLAG" : C.L2,
            })
    T, i = sandy.read_tab1(df, i)
    out.update({"QM" : T.C1, "QI" : T.C2, "LR" : T.L2, "NBT" : T.NBT, "INT" : T.INT, "E" : T.x, "XS" : T.y})
    return out

def write_mf3(sec):
    lines = sandy.write_cont(sec["ZA"], sec["AWR"], 0, sec["PFLAG"], 0, 0)
    lines += sandy.write_tab1(sec["QM"], sec["QI"], 0, sec["LR"], sec["NBT"], sec["INT"], sec["E"], sec["XS"])
    return "\n".join(sandy.write_eol(lines, sec["MAT"], 3, sec["MT"]))

def _read_errorr(text):
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    out = {"MAT" : MAT, "MF" : MF, "MT" : MT}
    i = 0
    L, i = read_list(str_list, i)
    out.update({"XS" : L.B})
    return out

def _read_groupr(text):
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    out = {"MAT" : MAT, "MF" : MF, "MT" : MT}
    i = 0
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2, "NL" : C.L1, "NZ" : C.L2, "LRFLAG" : C.N1, "NGN" : C.N2})
    groups = {}
    for ig in range(out["NGN"]):
        L, i = read_list(str_list, i)
        group = {"TEMPIN" : L.C1, "NG2" : L.L1, "IG2LO" : L.L2, "IG" : L.N2, "DATA" : L.B}
        groups[L.N2] = group
    out["GROUPS"] = groups
    return out