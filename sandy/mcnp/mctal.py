"""
This module contains all classes and functions to provide a API for a MCNP
*mctal* file.
"""

import pdb
import re
import logging
from io import StringIO

import numpy as np
import pandas as pd


__author__ = "Luca Fiorito"

__all__ = [
        "MctalTally",
        ]


class MctalTally():
    """
    Container for a MCTAL tally.

    Attributes
    ----------
    vals : `pandas.DataFrame`
        dataframe of tally values

        6 columns are provided:

           * f[c|t] : cells, surfaces or detector bins

                      fc = "cumulative" value is given
                      ft = "total" value is given

           * d : total vs. direct or flagged vs. unflagged bins

                 For detectors, n=2 unless there is an ND on the F5 card;
                 for cell and surface tallies, n=1 unless there is an SF or
                 CF card

           * u[c|t] : user bins

                      uc = "cumulative" value is given
                      ut = "total" value is given

           * s[c|t] : segment bins

                      sc = "cumulative" value is given
                      st = "total" value is given

           * m[c,t] : multiplier bins

                      mc = "cumulative" value is given
                      mt = "total" value is given

           * c[c,t] : cosine bins

                      cc = "cumulative" value is given
                      ct = "total" value is given

                      .. note:: `f` flag for point values is not yet
                                implemented

           * e[c,t] : energy bins

                      ec = "cumulative" value is given
                      et = "total" value is given

                      .. note:: `f` flag for point values is not yet
                                implemented

           * t[c,t] : time bins

                      tc = "cumulative" value is given
                      tt = "total" value is given

                      .. note:: `f` flag for point values is not yet
                                implemented

           * vals : tally values

                    The order is what a 9-dimensional Fortran array would have
                    if it were dimensioned (NF,...,NE,NT), where NF is the #
                    of cell, surface, or detector bins, NE is the # of energy
                    bins, and NT is the # of time bins

           * err : tally statistical error (relative)
    """

    def __init__(self, tally, vals, *args, **kwargs):
        self.data = vals
        self.tally = tally

    @property
    def _fkey(self):
        return "ft" if "ft" in self.data else "fc" if "fc" in self.data else "f"

    @property
    def _dkey(self):
        return "dt" if "dt" in self.data else "dc" if "dc" in self.data else "d"

    @property
    def _ukey(self):
        return "ut" if "ut" in self.data else "uc" if "uc" in self.data else "u"

    @property
    def _skey(self):
        return "st" if "st" in self.data else "sc" if "sc" in self.data else "s"

    @property
    def _mkey(self):
        return "mt" if "mt" in self.data else "mc" if "mc" in self.data else "m"

    @property
    def _ckey(self):
        return "ct" if "ct" in self.data else "cc" if "cc" in self.data else "c"

    @property
    def _ekey(self):
        return "et" if "et" in self.data else "ec" if "ec" in self.data else "e"

    @property
    def _tkey(self):
        return "tt" if "tt" in self.data else "tc" if "tc" in self.data else "t"

    def what_bins(self):
        msg = """
        .cell_bins       : "{0._fkey}"
        .detector_bins   : "{0._dkey}"
        .user_bins       : "{0._ukey}"
        .segment_bins    : "{0._skey}"
        .multiplier_bins : "{0._mkey}"
        .cosine_bins     : "{0._ckey}"
        .energy_bins     : "{0._ekey}"
        .time_bins       : "{0._tkey}"
        """.format(self)
        print(msg)

    @property
    def cell_bins(self):
        step = int(self.data.shape[0]/self.nf)
        bins = self.data.f.iloc[::step] \
                        .reset_index(drop=True) \
                        .rename("CELL_BINS")
        bins.index = bins.index.rename("POS")
        return bins

    @property
    def detector_bins(self):
        return self.data[self._dkey].unique().tolist()

    @property
    def user_bins(self):
        return self.data[self._ukey].unique().tolist()

    @property
    def segment_bins(self):
        return self.data[self._skey].unique().tolist()

    @property
    def multiplier_bins(self):
        return self.data[self._mkey].unique().tolist()

    @property
    def cosine_bins(self):
        return self.data[self._ckey].unique().tolist()

    @property
    def energy_bins(self):
        return self.data[self._ekey].unique().tolist()

    @property
    def time_bins(self):
        return self.data[self._tkey].unique().tolist()

    @classmethod
    def from_file(cls, filename="mctal", encoding="utf8", errors='ignore'):
        """
        Parse a mctal output file and return tallies.

        Parameters
        ----------
        filename : `str`, optional, default is `'mctal'`
            name of the mctal file
        encoding : `str`, optional, default is `'utf8'`
            option passed to python built-in function `open`
        errors : `str`, optional, default is `'ignore'`
            option passed to python built-in function `open`

        Returns
        -------
        `dict`
            dictionary of (tally numbers, `MctalTally` objects)
        """
        with open(filename, encoding=encoding, errors=errors) as f:
            string = f.read()
        blocks = re.split("\ntally", string)
        out = {}
        for block in blocks[1:]:
            tal = _from_block(block)
            out.update({tal.tally: tal})
        return out


def from_file(*args, **kwargs):
    """
    Same as `MctalTally.from_file`.
    """
    logging.warning("DEPRECATED! Use 'MctalTally.from_file'")
    return MctalTally.from_file(*args, **kwargs)


def _from_block(text):
    """
    Parse a MCTAL block and return tally.

    TO BE EXTENSIVELY RESTRUCTURED.
    """
    lines = text.splitlines()
    ipos = 0
    tally_number, i, j = list(map(int, lines[ipos].split()))[:3]
    ipos += 1
    if i < 0:
        particles = list(map(int, lines[ipos].split()))
        ipos += 1
    if lines[ipos][0] != "f":  # skip fc comment
        ipos += 1
    # F bins
    fdata = pd.read_csv(StringIO(lines[ipos]), header=None, sep="\s+", engine='python').iloc[0]
    nf = fdata[1]
    fbins = []
    while len(fbins) != nf:
        ipos += 1
        fbins += list(map(int, lines[ipos].split()))
    ipos += 1
    # D bins
    ddata = pd.read_csv(StringIO(lines[ipos]), header=None, sep="\s+", engine='python').iloc[0]
    nd = ddata[1]
    dbins = list(range(nd))
    ipos += 1
    # user bins
    udata = pd.read_csv(StringIO(lines[ipos]), header=None, sep="\s+", engine='python').iloc[0]
    nu = udata[1]
    if nu > 0:
        if udata[0].lower() == "ut":
            ubins = list(range(nu-1)) + ["total"]
        elif udata[0].lower() == "uc":
            ubins = list(range(nu-1)) + ["cumul"]
        else:
            ubins = list(range(nu))
    else:
        ubins = list(range(1))
    ipos += 1
    # segment bin
    sdata = pd.read_csv(StringIO(lines[ipos]), header=None, sep="\s+", engine='python').iloc[0]
    ns = sdata[1]
    if ns > 0:
        if sdata[0].lower() == "st":
            sbins = list(range(ns-1)) + ["total"]
        elif sdata[0].lower() == "sc":
            sbins = list(range(ns-1)) + ["cumul"]
        else:
            sbins = list(range(ns))
    else:
        sbins = list(range(1))
    ipos += 1
    # cosine bins
    mdata = pd.read_csv(StringIO(lines[ipos]), header=None, sep="\s+", engine='python').iloc[0]
    nm = mdata[1]
    if nm > 0:
        if mdata[0].lower() == "mt":
            mbins = list(range(nm-1)) + ["total"]
        elif mdata[0].lower() == "mc":
            mbins = list(range(nm-1)) + ["cumul"]
        else:
            mbins = list(range(nm))
    else:
        mbins = list(range(1))
    ipos += 1
    # cosine bins
    cdata = pd.read_csv(StringIO(lines[ipos]), header=None, sep="\s+", engine='python').iloc[0]
    nc = cdata[1]
    limit = int(np.ceil((nc-1)/6)) if cdata[0].lower() == "ct" or cdata[0].lower() == "cc" else int(np.ceil((nc)/6))
    cbins = []
    for jpos in range(limit):
        ipos += 1
        cbins += list(map(float, lines[ipos].split()))
    if nc > 0:
        if cdata[0].lower() == "ct":
            cbins += ["total"]
        elif cdata[0].lower() == "cc":
            cbins += ["cumul"]
    else:
        cbins = list(range(1))
    ipos += 1
    # energy bins
    edata = pd.read_csv(StringIO(lines[ipos]), header=None, sep="\s+", engine='python').iloc[0]
    ne = edata[1]
    limit = int(np.ceil((ne-1)/6)) if edata[0].lower() == "et" or edata[0].lower() == "ec" else int(np.ceil((ne)/6))
    ebins = []
    for jpos in range(limit):
        ipos += 1
        ebins += list(map(float, lines[ipos].split()))
    if ne > 0:
        if edata[0].lower() == "et":
            ebins += ["total"]
        elif edata[0].lower() == "ec":
            ebins += ["cumul"]
    else:
        ebins = list(range(1))
    ipos += 1
    # time bins
    tdata = pd.read_csv(StringIO(lines[ipos]), header=None, sep="\s+", engine='python').iloc[0]
    nt = tdata[1]
    limit = int(np.ceil((nt-1)/6)) if tdata[0].lower() == "tt" or tdata[0].lower() == "tc" else int(np.ceil((nt)/6))
    tbins = []
    for jpos in range(int(np.ceil((nt-1)/6))):
        ipos += jpos
        tbins += list(map(float, lines[ipos].split()))
    if nt > 0:
        if tdata[0].lower() == "tt":
            tbins += ["total"]
        elif tdata[0].lower() == "tc":
            tbins += ["cumul"]
    else:
        tbins = list(range(1))
    ipos += 1
    # values
    dim = 2 * max(1, nt) * max(1, ne) * max(1, nc) * max(1, nm) * max(1, ns) * max(1, nu) * max(1, nd) * max(1, nf)
    vals = []
    for jpos in range(int(np.ceil(dim/8))):
        ipos += 1
        try:
            A= list(map(float, lines[ipos].split()))
        except:
            pdb.set_trace()
        vals += list(map(float, lines[ipos].split()))
    # vals = np.array(vals).reshape(max(1, nt), max(1, ne), max(1, nc), max(1, nm), max(1, ns), max(1, nu), max(1, nd), max(1, nf), 2)
    df = pd.MultiIndex.from_product(
               [fbins, dbins, ubins, sbins, mbins, cbins, ebins, tbins],
               names=[fdata[0], ddata[0], udata[0], sdata[0], mdata[0], cdata[0], edata[0], tdata[0]]
                                   ).to_frame(index=False)
    df["vals"] = vals[::2]
    df["err"] = vals[1::2]
    tal = MctalTally(tally_number, df)
    tal.nf = int(nf)
    tal.nd = int(nd)
    tal.nu = int(nu)
    tal.ns = int(ns)
    tal.nm = int(nm)
    tal.nc = int(nc)
    tal.ne = int(ne)
    tal.nt = int(nt)
    return tal
