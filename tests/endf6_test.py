# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 09:17:18 2018

@author: Fiorito_L
"""
import pytest
import numpy as np

import pandas as pd

import sandy
from sandy.formats.endf6 import Endf6
from sandy.formats.mf1 import write as write_mf1
from sandy.formats.mf3 import write as write_mf3
from sandy.formats.mf4 import write as write_mf4
from sandy.formats.mf5 import write as write_mf5
from sandy.formats.mf8 import write as write_mf8
from sandy.data import U5, U8, Fe56, Pu9, H1, RDD

@pytest.fixture(scope="module")
def testPu9():
    tape = Endf6.from_text("\n".join(Pu9.endf6))
    assert (tape.index.get_level_values("MAT").unique() == 9437).all()
    return tape

@pytest.fixture(scope="module")
def testH1():
    tape = Endf6.from_text("\n".join(H1.pendf))
    assert (tape.index.get_level_values("MAT").unique() == 125).all()
    return tape

@pytest.fixture(scope="module")
def testFe56():
    tape = Endf6.from_text("\n".join(Fe56.endf6))
    assert (tape.index.get_level_values("MAT").unique() == 2631).all()
    return tape

@pytest.fixture(scope="module")
def testU5():
    tape = Endf6.from_text("\n".join(U5.nfy))
    assert (tape.index.get_level_values("MAT").unique() == 9228).all()
    return tape

@pytest.fixture(scope="module")
def testU8():
    tape = Endf6.from_text("\n".join(U8.endf6))
    assert (tape.index.get_level_values("MAT").unique() == 9237).all()
    return tape

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.info
def test_read_info(testPu9):
    S = testPu9.read_section(9437, 1, 451)
    text = write_mf1(S)
#    assert testPu9.TEXT.loc[9437,1,451] == text

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.info
def test_update_info(testPu9):
    testPu9.loc[9437,3,1].TEXT = "\n".join(testPu9.loc[9437,3,1].TEXT.splitlines()[:10]) + "\n"
    testPu9 = Endf6(testPu9.drop([(9437,3,102)]))
    new = testPu9.update_info()
    recordsold = testPu9.read_section(9437,1,451)["RECORDS"]
    recordsnew = new.read_section(9437,1,451)["RECORDS"]
    assert (3,102,147,1) in recordsold
    assert (3,102,147,1) not in recordsnew
    assert (3,1,188,1) in recordsold
    assert (3,1,10,1) in recordsnew

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.xs
def test_read_xs(testPu9):
    S = testPu9.read_section(9437, 3, 102)
    text = write_mf3(S)
#    assert testPu9.TEXT.loc[9437,3,102] == text

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.xs
def test_extract_xs(testH1):
    testH1.get_xs(listmat=[125], listmt=[1,2,102,4])
    xs = testH1.get_xs(listmat=[9437])
    assert xs.empty

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.xs
def test_reconstruct_xs(testH1):
    xs = testH1.get_xs()
    xs[(125,51)] = 1
    xs[(125,2)] = xs[(125,102)] = xs[(125,1)]
    SUM = xs[(125,1)].values*2 + 1
    rec_xs = xs.reconstruct_sums(drop=False)
    assert (rec_xs[(125,4)].values == 1).all()
    assert np.isclose(rec_xs[(125,1)].values,SUM).all()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.xs
def test_update_xs(testH1):
    xs = testH1.get_xs()
    xs[(125,2)] = 1
    new = testH1.update_xs(xs)
    assert (new.read_section(125,3,2)["XS"] == 1).all()
    assert (testH1.read_section(125,3,2)["XS"] != new.read_section(125,3,2)["XS"]).all()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.xs
@pytest.mark.cov
def test_read_xs_cov(testPu9):
    testPu9.read_section(9437, 33, 1)
    testPu9.read_section(9437, 33, 2)
    testPu9.read_section(9437, 33, 18)
    testPu9.read_section(9437, 31, 456)
    testPu9.read_section(9437, 33, 102)

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.xs
@pytest.mark.cov
def test_extract_xs_cov(testPu9):
    testPu9.get_xs_cov()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.nubar
def test_read_nubar(testPu9):
    for mt in (452,455,456):
        S = testPu9.read_section(9437, 1, mt)
        text = write_mf1(S)
#        assert testPu9.TEXT.loc[9437,1,mt] == text

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.nubar
def test_extract_nubar(testPu9):
    testPu9.get_nubar(listmt=[452,456])

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.nubar
def test_reconstruct_nubar(testPu9):
    nubar = testPu9.get_nubar()
    nubar[(9437,455)] = 1
    SUM = nubar[(9437,456)].values + 1
    rec_nubar = nubar.reconstruct_sums()
    assert np.isclose(rec_nubar[(9437,452)].values,SUM).all()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.nubar
def test_update_nubar(testPu9):
    nubar = testPu9.get_nubar()
    nubar[(9437,452)] = 1
    new = testPu9.update_nubar(nubar)
    assert (new.read_section(9437,1,452)["NUBAR"] == 1).all()
    assert (testPu9.read_section(9437,1,452)["NUBAR"] != new.read_section(9437,1,452)["NUBAR"]).all()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.lpc
def test_read_lpc(testFe56):
    S = testFe56.read_section(2631, 4, 2)
    assert S["LTT"] == 3
    assert "LPC" in S
    assert S["LPC"]["NE"] == len(S["LPC"]["E"]) == 1782
    assert "TAB" in S
    assert S["TAB"]["NE"] == len(S["TAB"]["E"]) == 28
    text = write_mf4(S)
#    assert testFe56.TEXT.loc[2631,4,2] == text
    for mt in range(51,83):
        S = testFe56.read_section(2631, 4, mt)
        text = write_mf4(S)
#        assert testFe56.TEXT.loc[2631,4,mt] == text

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.lpc
def test_extract_lpc(testFe56):
    testFe56.get_lpc()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.lpc
def test_convert_lpc(testFe56):
    Lpc = testFe56.get_lpc()
    C = Lpc.to_tab(2631, 2, 1e-5)
    assert (C == 0.5).all()
    C = Lpc.to_tab(2631, 2, 1e4)
    assert (C >= 0).all()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.lpc
def test_addpoints_lpc(testFe56):
    lpcold = testFe56.get_lpc(listmt=[2])
    lpcnew = lpcold.add_points([1e4, 1e5, 1e6])
    assert 1e4 in lpcnew.loc[2631,2].index
    assert 1e4 not in lpcold.loc[2631,2].index
    assert 1e5 in lpcnew.loc[2631,2].index
    assert 1e5 not in lpcold.loc[2631,2].index
    assert 1e6 in lpcnew.loc[2631,2].index
    assert 1e6 in lpcold.loc[2631,2].index
    lpcnew = lpcold.add_points([])
    assert (lpcnew.index == lpcold.index).all()
    lpcold = testFe56.get_lpc(listmt=[810])
    lpcnew = lpcold.add_points([1e4])
    assert 1e4 not in lpcnew.loc[2631,810].index
    assert 1e4 not in lpcold.loc[2631,810].index

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.lpc
def test_update_lpc(testFe56):
    lpc = testFe56.get_lpc()
    new = testFe56.update_lpc(lpc)
#    assert new.TEXT[2631,4,2] == testFe56.TEXT[2631,4,2]
    lpc.loc[2631,2,1e5:2e7] = 0
    new = testFe56.update_lpc(lpc)
#    assert new.TEXT[2631,4,2] != testFe56.TEXT[2631,4,2]
    new_sec = new.read_section(2631,4,2)
    assert (np.array(new_sec["LPC"]["E"][2e7]["COEFF"]) == 0).all()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.lpc
@pytest.mark.cov
def test_read_lpc_cov(testFe56):
    S = testFe56.read_section(2631, 34, 2)
    assert S["NMT1"] == len(S["REAC"]) == 1
    assert len(S["REAC"][0,2]["P"]) == S["REAC"][0,2]["NL"]*(S["REAC"][0,2]["NL"]+1)//2
    assert len(S["REAC"][0,2]["P"][1,1]["NI"]) == 3

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.lpc
@pytest.mark.cov
def test_extract_lpc_cov(testFe56):
    C = testFe56.get_lpc_cov()
    assert C.index.names == C.columns.names == ['MAT', 'MT', 'L', 'E']
    assert (C.values == C.values.T).all()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.lpc
def test_perturb_lpc(testU8):
    energies =  [1e-05, 0.0253, 2e2, 1e3, 1e4, 6e4, 4e5, 7e5, 1.2e6,  1.8e6, 2.6e6, 3.2e6, 4e6, 4.6e6, 5.4e6, 6e6, 3e7]
    le = len(energies)
    lpc = testU8.get_lpc()
    smp = np.random.normal(1.0, 0.5, le*6)
    ls = [1]*le + [2]*le + [3]*le + [4]*le + [5]*le + [6]*le
    perts = pd.DataFrame.from_dict(dict(E=energies*6, SMP=smp.tolist(), L=ls, MT=2, MAT=9237)).set_index(["MAT","MT","L","E"]).SMP
    out = lpc.perturb(perts)
    assert ((out/lpc).stack() >= 0).all()
    assert ((out/lpc).stack() <= 2).all()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.lpc
@pytest.mark.cov
@pytest.mark.plot
def test_plot_std_lpc(testU8):
    lpccov = testU8.get_lpc_cov()
    ax = lpccov.plot_std(display=False)    
    assert ax is not None

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.lpc
@pytest.mark.tpd
def test_lpc_to_tpd(testU8):
    lpc = testU8.get_lpc()
    tpd = lpc.to_tpd()
    assert (tpd.index == lpc.index).all()
    assert (tpd.loc[9237,2,1e-5] == 0.5).all()
    assert np.isclose(tpd.loc[9237,51,1e7][-0.94], 0.210904)
    assert np.isclose(tpd.loc[9237,51,5e4][0.99], 0.412180)
#    A=testU8.update_tpd(tpd)
#    pytest.set_trace()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.chi
def test_read_chi(testPu9):
    S = testPu9.read_section(9437, 5, 18)
    text = write_mf5(S)
#    assert testPu9.TEXT.loc[9437,5,18] == text

@pytest.fixture(scope="module")
@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.chi
def testPu9chi(testPu9):
    return testPu9.get_edistr()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.chi
def test_extract_chi(testPu9chi):
    testPu9chi.add_points([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5])

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.chi
def test_normalize_chi(testPu9chi):
    chi = testPu9chi.normalize()
    for i,v in chi.iterrows():
        dx = v.index.values[1:] - v.index.values[:-1]
        y = (v.values[1:]+v.values[:-1])/2
        assert np.isclose(y.dot(dx), 1, rtol=1e-10)

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.chi
def test_update_chi(testPu9, testPu9chi):
    testPu9chi.loc[(9437,18,1,1e-5)] = 1
    new = testPu9.update_edistr(testPu9chi)
    new_sec = new.read_section(9437,5,18)
    assert (new_sec["PDISTR"][1]["EIN"][1e-5]["EDISTR"] == 1).all()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.chi
def test_add_points_chi(testPu9, testPu9chi):
    extra_points = [1e-5, 1e-4, 1e-3]
    edistr = testPu9chi.add_points(extra_points)
    new = testPu9.update_edistr(edistr)
    new_sec = new.read_section(9437,5,18)
    for e in extra_points:
        assert e in new_sec["PDISTR"][1]["EIN"]

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.chi
@pytest.mark.cov
def test_read_chi_cov(testPu9):
    S = testPu9.read_section(9437, 35, 18)
    assert len(S["SUB"]) == S["NK"] == 8
    for k,v in S["SUB"].items():
        assert v["LS"] == 1
        assert v["LB"] == 7

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.chi
@pytest.mark.cov
def test_extract_chi_cov(testPu9):
    C = testPu9.get_edistr_cov()
    assert C.index.names == ['MAT', 'MT', 'ELO', 'EHI', 'EOUT']
    assert C.columns.names == ['MAT', 'MT', 'ELO', 'EHI', 'EOUT']
    assert (C.values == C.values.T).all()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.chi
def test_perturb_chi(testPu9):
    edistr = sandy.Edistr(testPu9.get_edistr().query("EIN==1e5 | EIN==2e6"))
    mat = 9437
    mt = 18
    std = edistr.loc[mat,mt,1,2e6]*0.6
    cov = np.diag(std**2)
    smp = np.random.multivariate_normal([0]*len(std), cov, 1).flatten()
    perts = pd.DataFrame.from_dict(dict(EOUT=std.index, SMP=smp.tolist(), ELO=1e-5, EHI=3e7, MT=mt, MAT=mat)).set_index(["MAT","MT","ELO","EHI", "EOUT"]).SMP
    out = edistr.perturb(perts, normalize=False)
    assert (out>=0).all().all()
    assert (out.iloc[0] <= out.iloc[0]*2).all()
    assert (out.iloc[-1] <= out.iloc[-1]*2).all()
#    N0 = (out.iloc[0] == 0).sum()
#    N2 = (out.iloc[0] == 2*edistr.iloc[0]).sum()
#    assert N0 > 0
#    assert N0 in range(N2-10, N2+11)
#    N0 = (out.iloc[1] == 0).sum()
#    N2 = (out.iloc[1] == 2*edistr.iloc[1]).sum()
#    assert N0 > 0
#    assert N0 in range(N2-10, N2+11)

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.fy
def test_read_fy(testU5):
    S = testU5.read_section(9228, 8, 454)
    text = write_mf8(S)
#    assert testU5.TEXT.loc[9228,8,454] == text

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.fy
def test_extract_fy(testU5):
    fy = testU5.get_fy(listmt=[454])

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.write
def test_write_to_string(testH1):
    string = testH1.write_string()
    newtape = Endf6.from_text(string)
    assert testH1.equals(newtape)
    
