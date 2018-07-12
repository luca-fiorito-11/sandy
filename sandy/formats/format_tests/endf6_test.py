# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 09:17:18 2018

@author: Fiorito_L
"""
import pytest
import numpy as np

from .. import Endf6
from ..MF1 import write as write_mf1
from ..MF3 import write as write_mf3
from ..MF4 import write as write_mf4
from ..MF5 import write as write_mf5
from ...data import Pu9
from ...data import H1
from ...data import Fe56

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

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.info
def test_read_info(testPu9):
    S = testPu9.read_section(9437, 1, 451)
    text = write_mf1(S)
    assert testPu9.TEXT.loc[9437,1,451] == text

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
    assert testPu9.TEXT.loc[9437,3,102] == text

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
        assert testPu9.TEXT.loc[9437,1,mt] == text

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
    assert testFe56.TEXT.loc[2631,4,2] == text
    for mt in range(51,83):
        S = testFe56.read_section(2631, 4, mt)
        text = write_mf4(S)
        assert testFe56.TEXT.loc[2631,4,mt] == text

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
    assert new.TEXT[2631,4,2] == testFe56.TEXT[2631,4,2]
    lpc.loc[2631,2,1e5:2e7] = 0
    new = testFe56.update_lpc(lpc)
    assert new.TEXT[2631,4,2] != testFe56.TEXT[2631,4,2]
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
@pytest.mark.chi
def test_read_chi(testPu9):
    S = testPu9.read_section(9437, 5, 18)
    text = write_mf5(S)
    assert testPu9.TEXT.loc[9437,5,18] == text

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
@pytest.mark.write
def test_write_to_string(testH1):
    string = testH1.write_string()
    newtape = Endf6.from_text(string)
    assert testH1.equals(newtape)