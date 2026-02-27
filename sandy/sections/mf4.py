__author__ = "Luca Fiorito"

def read_mf4(tape, mat, mt):
    from ..records import read_cont, read_list, read_tab1, read_tab2

    mf = 4
    df = tape._get_section_df(mat, mf, mt)
    out = {"MAT" : mat, "MF" : mf, "MT" : mt}
    i = 0
    C, i = read_cont(df, i)
    out.update({
            "ZA" : C.C1,
            "AWR" : C.C2,
            "LTT" : C.L2,
            })
    C, i = read_cont(df, i)
    out.update({
            "LI" : C.L1,
            "LCT" : C.L2,
            })
    # polynomial distributions
    if out["LTT"] in (1, 3):
        T2, i = read_tab2(df, i)
        energy_points = T2.NZ
        edistr = {}
        for j in range(energy_points):
            L, i = read_list(df, i)
            edistr[L.C2] = {
                    "COEFF" : L.B,
                    "T" : L.C1,
                    "LT" : L.L1
                    }
        out["LPC"] = {
                "NE" : T2.NZ,
                "NBT" : T2.NBT,
                "INT" : T2.INT, 
                "E" : edistr}
    # tabulated distributions
    if out["LTT"] in (2, 3):
        T2, i = read_tab2(df, i)
        energy_points = T2.NZ
        edistr = {}
        for j in range(energy_points):
            T, i = read_tab1(df, i)
            edistr[T.C2] = {
                    "T" : T.C1,
                    "LT" : T.L1,
                    "NBT" : T.NBT,
                    "INT" : T.INT,
                    "MU" : T.x,
                    "ADISTR" : T.y,
                    }
        out["TAB"] = {
                "NE" : T2.NZ,
                "NBT" : T2.NBT,
                "INT" : T2.INT, 
                "E" : edistr}
    return out



def write_mf4(sec):
    from ..records import write_cont, write_list, write_tab1, write_tab2, write_eol

    if "LCP" in sec and "TAB" in sec:
        sec["LTT"] == 3
        sec["LI"] == 0
    elif "LCP" in sec:
        sec["LTT"] == 1
        sec["LI"] == 0
    elif "TAB" in sec:
        sec["LTT"] == 2
        sec["LI"] == 0
    else:
        sec["LTT"] == 0
        sec["LI"] == 1
    lines = write_cont(sec["ZA"], sec["AWR"], 0, sec["LTT"], 0, 0)
    lines += write_cont(0, sec["AWR"], sec["LI"], sec["LCT"], 0, 0)
    if sec["LTT"] in (1,3):
        lpc = sec["LPC"]
        lines += write_tab2(0, 0, 0, 0, len(lpc["E"]), lpc["NBT"], lpc["INT"])
        for e,sub in sorted(lpc["E"].items()):
            lines += write_list(sub["T"], e, sub["LT"], 0, 0, sub["COEFF"])
    if sec["LTT"] in (2,3):
        tab = sec["TAB"]
        lines += write_tab2(0, 0, 0, 0, len(tab["E"]), tab["NBT"], tab["INT"])
        for e,sub in sorted(tab["E"].items()):
            lines += write_tab1(sub["T"], e, sub["LT"], 0, sub["NBT"], sub["INT"], sub["MU"], sub["ADISTR"])
    return "\n".join(write_eol(lines, sec["MAT"], 4, sec["MT"]))
