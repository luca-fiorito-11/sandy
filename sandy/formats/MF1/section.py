# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 09:23:33 2018

@author: fiorito_l
"""
from ..records import read_cont, read_tab1, read_control, read_text, read_list, write_cont, write_tab1, write_list
from ..utils import Section
import sys

def read(text):
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    if MT == 451:
        return read_info(text)
    elif MT in (452,455,456):
        return read_nubar(text)

def write(sec):
    if sec["MT"] == 451:
        return write_info(sec)
    elif sec["MT"] in (452,455,456):
        return write_nubar(sec)

def read_errorr(text):
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    out = {"MAT" : MAT, "MF" : MF, "MT" : MT}
    i = 0
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2, "ERRFLAG" :C.N1})
    L, i = read_list(str_list, i)
    out.update({"EG" : L.B})
    return out

def read_info(text):
    from sandy.csvq import elements, metastates
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    out = {"MAT" : MAT, "MF" : MF, "MT" : MT}
    i = 0
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2, "LRP" : C.L1, "LFI" : C.L2, "NLIB" :C.N1, "NMOD" : C.N2})
    C, i = read_cont(str_list, i)
    out.update({"ELIS" : C.C1, "STA" : C.C2, "LIS" : C.L1, "LISO" : C.L2, "NFOR" : C.N2})
    C, i = read_cont(str_list, i)
    out.update({"AWI" : C.C1, "EMAX" : C.C2, "LREL" : C.L1, "NSUB" : C.N1, "NVER" : C.N2})
    C, i = read_cont(str_list, i)
    out.update({"TEMP" : C.C2, "LDRV" : C.L1, "NWD" : C.N1, "NXC" : C.N2})
    out["Z"] = int(out["ZA"]//1000)
    out["A"] = int(out["ZA"]-out["Z"]*1000)
    out["SYM"] = elements["SYM"].to_dict()[out["Z"]]
    out["M"] = metastates["META"].to_dict()[out["LISO"]]
    out["TAG"] = "{}-{}-{}{}".format(out["Z"], out["SYM"], out["A"], out["M"])
    TEXT = []
    for j in range(out["NWD"]):
        T, i = read_text(str_list, i)
        TEXT.append(T)
    out.update({ "TEXT" : TEXT })
    # This part is not given in PENDF files
    if out["LRP"] != 2:
#        groups = TEXT[0][:11].split("-")
#        out["Z"] = int(groups[0])
#        out["SYM"] = groups[1].strip()
#        out["A"] = re.sub(r"\D", "", groups[2])
#        out["M"] = re.sub(r"[0-9\s]", "", groups[2].lower())
#        if not out["M"]: out["M"] = 'g'
#            out["A"] = int(TEXT[0][7:10])
#            out["M"] =  'g' if TEXT[0][10:11] is ' ' else TEXT[0][10:11].lower()
#            out["Z"] = int(TEXT[0][:3])
#            out["SYM"] = TEXT[0][4:6].rstrip()
#            out["A"] = int(TEXT[0][7:10])
#            out["M"] =  'g' if TEXT[0][10:11] is ' ' else TEXT[0][10:11].lower()
        out['ALAB'] = TEXT[0][11:22]
        out['EDATE'] = TEXT[0][22:32]
        out['AUTH'] = TEXT[0][33:66]
        out['REF'] = TEXT[1][1:22]
        out['DDATE'] = TEXT[1][22:32]
        out['RDATE'] = TEXT[1][33:43]
        out['ENDATE'] = TEXT[1][55:63]
        out['LIBVER'] = TEXT[2][:22].strip('- ')
        out['SUB'] = TEXT[3].strip('- ')
        out['FOR'] = TEXT[4].strip('- ')
    out.update({ "RECORDS" : [] })
    for j in range(out["NXC"]):
        C, i = read_cont(str_list, i)
        out["RECORDS"].append((C.L1,C.L2,C.N1,C.N2))
    return Section(out)

def write_info(sec):
    text = write_cont(sec["ZA"], sec["AWR"], sec["LRP"], sec["LFI"], sec["NLIB"], sec["NMOD"])
    text += write_cont(sec["ELIS"], sec["STA"], sec["LIS"], sec["LISO"], 0, sec["NFOR"])
    text += write_cont(sec["AWI"], sec["EMAX"], sec["LREL"], 0, sec["NSUB"], sec["NVER"])
    text += write_cont(0, sec["TEMP"], sec["LDRV"], 0, len(sec["TEXT"]), len(sec["RECORDS"]))
    text += sec["TEXT"]
    text += [ " "*22 + "{:>11}{:>11}{:>11}{:>11}".format(*x) for x in sec["RECORDS"]]
    TextOut = []; iline = 1
    for line in text:
        if iline > 99999:
            iline = 1
        TextOut.append("{:<66}{:4}{:2}{:3}{:5}\n".format(line, sec["MAT"], sec["MF"], sec["MT"], iline))
        iline += 1
    return "".join(TextOut)

def read_nubar(text):
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    out = {"MAT" : MAT, "MF" : MF, "MT" : MT}
    i = 0
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2, "LDG" : C.L1, "LNU" : C.L2})
    if out["MT"] == 455:
        if out["LDG"] == 0:
            L, i = read_list(str_list, i)
            out.update({ "NNF" : L.NPL, "LAMBDAS" : L.B })
        elif out["LDG"] == 1:
            # None found in JEFF33 and ENDFB8, hence not implemented
            sys.exit("ERROR: Not implemented format")
            pass
    if out["LNU"] == 1:
        # None found in JEFF33 and ENDFB8 neither for MT455 nor for MT456
        L, i = read_list(str_list, i)
        out.update({ "NC" : L.NPL, "C" : L.B})
    else:
        # JEFF33 and ENDFB8 only have lin-lin interpolation schemes
        T, i = read_tab1(str_list, i)
        out.update({"NBT" : T.NBT, "INT" : T.INT, "E" : T.x, "NUBAR" : T.y})
    return out


def write_nubar(sec):
    text = write_cont(sec["ZA"], sec["AWR"], sec["LDG"], sec["LNU"], 0, 0)
    if sec["MT"] == 455:
        if sec["LDG"] == 0:
            text += write_list(0, 0, 0, 0, 0, sec["LAMBDAS"])
        elif sec["LDG"] == 1:
            sys.exit("ERROR: Not found in JEFF33 and ENDFB8, hence not implemented")
    if sec["LNU"] == 1:
        text += write_list(0, 0, 0, 0, 0, sec["C"])
    else:
        text += write_tab1(0, 0, 0, 0, sec["NBT"], sec["INT"], sec["E"], sec["NUBAR"])
    TextOut = []; iline = 1
    for line in text:
        if iline > 99999:
            iline = 1
        TextOut.append("{:<66}{:4}{:2}{:3}{:5}\n".format(line, sec["MAT"], sec["MF"], sec["MT"], iline))
        iline += 1
    return "".join(TextOut)