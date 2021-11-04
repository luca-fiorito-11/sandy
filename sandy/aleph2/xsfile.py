import pytest
import re
import h5py
import logging
from io import StringIO
import itertools
import os
from copy import deepcopy

import numpy as np
import pandas as pd

import sandy

__author__ = "Luca Fiorito"

__all__ = [
        "AlephFile",
        ]


test_aleph2xs = """              952410                   4                   4  200.00000000000000
  5.000000000000E+00
  1.000000000000E-11  1.000000000000E-05  1.000000000000E-01  1.000000000000E+01
              952420                                       0  5.000000000000E+00
  1.0000000000000000  2.0000000000000000  3.0000000000000000  4.0000000000000000
                  18                   0  2.000000000000E+02
  1.0000000000000000  2.0000000000000000  3.0000000000000000  4.0000000000000000
              942410                   2  5.000000000000E+00
  6.0000000000000000  5.0000000000000000
               10010                   2  7.000000000000E-01
  6.0000000000000000  5.0000000000000000
"""

pattern_alephxsfile = "^.*_(?P<lib>.*)\\.(?P<tmp>.*)al"


def _lib_from_filename(file):
    fname = os.path.split(file)[1]
    match = re.match(pattern_alephxsfile, fname)
    msg = f"cannot automatically determine 'library' from filename '{fname}'"
    if not match:
        raise sandy.Error(msg)
    ext = match.group("lib")
    if ext.lower() in ext2lib:
        library = ext2lib[ext.lower()]
    else:
        raise sandy.Error(msg)
    return library


def _tmp_from_filename(file):
    fname = os.path.split(file)[1]
    match = re.match(pattern_alephxsfile, fname)
    msg = f"cannot automatically determine 'temperature' from filename '{fname}'"
    if not match:
        raise sandy.Error(msg)
    ext = match.group("tmp")
    if ext in ext2tmp:
        temperature = ext2tmp[ext]
    elif ext in ext2tmp_meta:
        temperature = ext2tmp_meta[ext]
    else:
        raise sandy.Error(msg)
    return temperature


class AlephFile():
    """
    Examples
    --------
    >>> inp = AlephFile.from_text(test_aleph2xs)
    >>> inp.data
    {'nuclide': 952410,
     'energies': array([1.e-11, 1.e-05, 1.e-01, 1.e+01]),
     'fission_qvalue': 200.0,
     'capture_qvalue': 5.0,
     'reactions': {952420: {'qval': 5.0, 'xs': array([1., 2., 3., 4.])},
      18: {'qval': 200.0, 'xs': array([1., 2., 3., 4.])},
      942410: {'qval': 5.0, 'xs': array([6., 5.])},
      10010: {'qval': 0.7, 'xs': array([6., 5.])},
      4581: {'qval': 0.0, 'xs': array([200., 400., 600., 800.])},
      4582: {'qval': 0.0, 'xs': array([ 5., 10., 15., 20.])}}}
    """

    def __repr__(self):
        return f"ALEPH file for nuclide: {self.nuclide}"

    def __init__(self, data):
        self.data = data

    # def __init__(self, text):
    #     self.data = read(text)
    #     self.add_fission_xsenergy()
    #     self.add_capture_xsenergy()

    @property
    def nuclide(self):
        return self.data["nuclide"]

    @property
    def npoints(self):
        return self.data["energies"].size

    def copy(self):
        data = deepcopy(self.data)
        return self.__class__(data)

    def add_capture_xsenergy(self):
        """
        Add

        Returns
        -------
        None.
        """
        mt = 4582
        zap = int(np.floor(self.nuclide / 10) * 10 + 10)
        reacts = self.data["reactions"]
        qval = self.data["capture_qvalue"]
        capture_products = [m for m in range(zap, zap + 10) if m in reacts]
        if capture_products:
            capture_xs = [reacts[i]["xs"].tolist() for i in capture_products]
            if np.unique([len(x) for x in capture_xs]).size > 1:
                msg = "capture cross sections have different " +\
                     f"lengths for '{self.nuclide}'. " +\
                      "Will use only ground state"
                logging.warning(msg)
                capture_xs = [capture_xs[0]]
            xsc = np.array(capture_xs).sum(axis=0)
            xsce = xsc * qval
        reaction = {
            "qval": 0.,
            "xs": xsce,
            }
        out = self.copy()
        out.data["reactions"][mt] = reaction
        return out

    def add_fission_xsenergy(self):
        """
        Add

        Returns
        -------
        None.

        """
        mt = 4581
        mtf = 18
        mtef = 1458
        reacts = self.data["reactions"]
        qval = self.data["fission_qvalue"]
        pad = sandy.utils.pad_from_beginning_fast
        if mtf in reacts:
            xs = pad([reacts[mtf]["xs"]], self.npoints)[0]
        else:
            xs = np.zeros((self.npoints,))
        if mtef in reacts:
            xe = pad([reacts[mtef]["xs"]], self.npoints)[0]
            xsfe = xs * xe
        else:
            xsfe = xs * qval
        reaction = {
            "qval": 0.,
            "xs": xsfe,
            }
        out = self.copy()
        out.data["reactions"][mt] = reaction
        return out

    def add_reaction(self, key, array, qval=0.):
        """

        Parameters
        ----------
        key : `int`
            reaction product ZAM identifier (daughter).
        array : `numpy.ndarray`
            cross section values over original energy grid, or a selection for
            threshold reactions.
        qval : `float`, optional
            reaction qvalue. The default is 0..

        Raises
        ------
        ValueError
            If cross section array contains more entries than the energy
            points.

        Returns
        -------
        None.

        Examples
        --------
        >>> inp = AlephFile.from_text(test_aleph2xs)
        >>> new = [3, 4, 5]
        >>> inp2 = inp.add_reaction(100, new)
        >>> inp2.data
        {'nuclide': 952410,
         'energies': array([1.e-11, 1.e-05, 1.e-01, 1.e+01]),
         'fission_qvalue': 200.0,
         'capture_qvalue': 5.0,
         'reactions': {952420: {'qval': 5.0, 'xs': array([1., 2., 3., 4.])},
          18: {'qval': 200.0, 'xs': array([1., 2., 3., 4.])},
          942410: {'qval': 5.0, 'xs': array([6., 5.])},
          10010: {'qval': 0.7, 'xs': array([6., 5.])},
          4581: {'qval': 0.0, 'xs': array([200., 400., 600., 800.])},
          4582: {'qval': 0.0, 'xs': array([ 5., 10., 15., 20.])},
          100: {'qval': 0.0, 'xs': array([3., 4., 5.])}}}
        """
        if len(array) > self.npoints:
            raise ValueError("reaction over too many points")
        reaction = {
            "qval": qval,
            "xs": np.array(array, dtype=float),
            }
        out = self.copy()
        out.data["reactions"][key] = reaction
        return out

    def to_dataframe(self):
        """
        Order cross sections into a dataframe with energy as index and
        parent/daughter ZAM numbers as columns.

        Returns
        -------
        df : `pandas.DataFrame`
            dtaframe of energy-dependent cross sections.

        Notes
        -----
        .. note:: threshold reactions are padded with zero from the first
                  energy point to the last before the threshold.

        Examples
        --------
        >>> AlephFile.from_text(test_aleph2xs).to_dataframe()
        PARENT           952410
        DAUGHTER         952420      18          942410      10010       4581        4582
        E
        1.00000e-11 1.00000e+00 1.00000e+00 0.00000e+00 0.00000e+00 2.00000e+02 5.00000e+00
        1.00000e-05 2.00000e+00 2.00000e+00 0.00000e+00 0.00000e+00 4.00000e+02 1.00000e+01
        1.00000e-01 3.00000e+00 3.00000e+00 6.00000e+00 6.00000e+00 6.00000e+02 1.50000e+01
        1.00000e+01 4.00000e+00 4.00000e+00 5.00000e+00 5.00000e+00 8.00000e+02 2.00000e+01

        PARENT           952410              ...                        
        DAUGHTER         952420      18      ...      4581        4582  
        E                                    ...                        
        1.00000e-11 1.00000e+00 1.00000e+00  ... 2.00000e+02 5.00000e+00
        1.00000e-05 2.00000e+00 2.00000e+00  ... 4.00000e+02 1.00000e+01
        1.00000e-01 3.00000e+00 3.00000e+00  ... 6.00000e+02 1.50000e+01
        1.00000e+01 4.00000e+00 4.00000e+00  ... 8.00000e+02 2.00000e+01
        """
        keys = []
        vals = []
        parent = self.data['nuclide']
        index = pd.Index(self.data['energies'], name="E")
        for key, item in self.data['reactions'].items():
            keys.append((parent, key))
            vals.append(item["xs"])
        matrix = sandy.utils.pad_from_beginning_fast(vals, index.size)
        columns = pd.MultiIndex.from_tuples(keys, names=["PARENT", "DAUGHTER"])
        df = pd.DataFrame(matrix.T, index=index, columns=columns)
        return df

    def to_string(self):
        """
        Write file content into sting.

        Returns
        -------
        `str`
            string with file content.

        Examples
        --------
        >>> print(AlephFile.from_text(test_aleph2xs).to_string())
                      952410                   4                   4  2.000000000000E+02
          5.000000000000E+00                                                            
          1.000000000000E-11  1.000000000000E-05  1.000000000000E-01  1.000000000000E+01
                          18                   0  2.000000000000E+02                    
          1.000000000000E+00  2.000000000000E+00  3.000000000000E+00  4.000000000000E+00
                       10010                   2  7.000000000000E-01                    
          6.000000000000E+00  5.000000000000E+00                                        
                      942410                   2  5.000000000000E+00                    
          6.000000000000E+00  5.000000000000E+00                                        
                      952420                   0  5.000000000000E+00                    
          1.000000000000E+00  2.000000000000E+00  3.000000000000E+00  4.000000000000E+00
        """
        data = self.data
        lines = []
        reacts = data["reactions"].copy()
        if 4581 in reacts:
            del reacts[4581]
        if 4582 in reacts:
            del reacts[4582]
        nreacts = len(reacts)
        line = write_line(
            data["nuclide"],
            nreacts,
            self.npoints,
            data['fission_qvalue'],
            )
        lines.append(line)
        line = write_line(
            data['capture_qvalue'],
            )
        lines.append(line)
        lines += write_array(data["energies"])
        for k, v in sorted(reacts.items()):
            qval = v["qval"]
            xs = v["xs"]
            lines += write_xs(k, xs, self.npoints, qval)
        return "\n".join(lines)

    def to_hdf5(self, h5file, library, temperature, dtype=None):
        """
        Read ALEPH xs file in ascii format and append content
        to existing/new hdf5 file.

        Datasets of xs are appended to the hdf5 file hierarchical
        sructure with library name, xs temperature and isotope zam
        as group keys.

        Parameters
        ----------
        h5file : `str`
            hdf5 filename.
            If the file does not exist it is created.
        library : `str`
            library identifier
        temperature : `float`
            temperature in K
        dtype : optional, default is None

        Examples
        --------
        >>> AlephFile.from_text(test_aleph2xs).to_hdf5("test.h5",
        ...                                  library="test",
        ...                                  temperature=900)
        Adding 'test/xs/900.0/952410' to file'test.h5'...
        """
        # filename was initially kept to extract temperature and library
        # in `to_hdf5`.
        # Now these data have to be supplied manually
        dct = self.data
        key = f"{library:s}/xs/{temperature:.1f}/{self.nuclide:d}"
        print(f"Adding '{key}' to file'{h5file}'...")
        hdf = sandy.H5File(h5file)   # use aleph h5file module
        hdf.open(mode="a")
        if key in hdf.data:
            logging.warning(f"hdf5 dataset '{key}' already exists "
                            "and will be replaced")
            del hdf.data[key]
        group = hdf.data.create_group(key)
        # some attributes are redundant as they are already in group key
        group.attrs["nuclide"] = dct["nuclide"]   # this is redundant
        group.attrs["temperature"] = float(temperature)
        group.attrs["library"] = library          # this is redundant
        group.attrs["fission_qvalue"] = dct["fission_qvalue"]
        group.attrs["capture_qvalue"] = dct["capture_qvalue"]
        # add energy to the reaction list
        lst = [dct["energies"]]
        products = []
        for k, v in dct["reactions"].items():
            xs = v["xs"]
            lst.append(xs)
            products.append(k)
        lengths = np.array(list(itertools.accumulate([v.size for v in lst])))
        group.create_dataset(
                "reactions",
                data=np.concatenate(lst),
                dtype=dtype,
                )
        group.create_dataset(
                "products",
                data=[0] + products,
                dtype=int,
                )
        group.create_dataset(
                "lengths",
                data=lengths,
                dtype=int,
                )
        hdf.close()

    @classmethod
    def from_file(cls, file):
        """
        Read ALEPH-2 cross sections from file.

        Parameters
        ----------
        file : `str`
            name of the file containing ALEPH-2 cross section data.
        """
        # reading from file is removed from init to be compatible with
        # StringIO
        with open(file) as f:
            text = f.read()
        return cls.from_text(text)

    @classmethod
    def from_text(cls, text):
        """
        

        Parameters
        ----------
        cls : TYPE
            DESCRIPTION.
        text : TYPE
            DESCRIPTION.

        Returns
        -------
        data : TYPE
            DESCRIPTION.

        """
        data = cls(read(text)).add_fission_xsenergy().add_capture_xsenergy()
        return data


def read_line(line, **kwargs):
    """
    Read one line of an aleph2 xs file.

    Parameters
    ----------
    line : `str`
        line of file

    Returns
    -------
    `list` or scalar
        list of values in line. Return scalar if only one value is found

    Examples
    --------
    Read file line by line.
    >>> line = "              952410                   3                   4  200.00000000000000"
    >>> read_line(line)
    [952410.0, 3.0, 4.0, 200.0]

    #>>> line = "  5.000000000000E+00"
    #>>> read_line(line, dtype=int)
    #5
    """
    stream = StringIO(line)
    vals = np.genfromtxt(stream, **kwargs).tolist()
    return vals


def read_array(iterable, size, per_line=4):
    """
    

    Parameters
    ----------
    iterable : TYPE
        DESCRIPTION.
    size : TYPE
        DESCRIPTION.
    per_line : TYPE, optional
        DESCRIPTION. The default is 4.

    Returns
    -------
    TYPE
        DESCRIPTION.

    Examples
    --------
    >>> it = iter(test_aleph2xs.splitlines())
    >>> skipped = next(it)
    >>> skipped = next(it)
    >>> read_array(it, 4)
    array([1.e-11, 1.e-05, 1.e-01, 1.e+01])
    >>> next(it)
    '              952420                                       0  5.000000000000E+00'
    """
    nlines = int(np.ceil(size / per_line))
    array = []
    for x in range(nlines):
        vals = read_line(next(iterable))
        if isinstance(vals, list):
            array += vals
        else:
            array.append(vals)
    return np.array(array)


def read_lines(lines):
    """
    Read ALEPH xs file in ascii format and append content
    to existing/new hdf5 file.

    Datasets of xs are appended to the hdf5 file hierarchical
    sructure with library name, xs temperature and isotope zam
    as group keys.

    Parameters
    ----------
    lines : TYPE
        DESCRIPTION.

    Examples
    --------
    Read test case.
    >>> lines = test_aleph2xs.splitlines()
    >>> read_lines(lines)
    {'nuclide': 952410,
     'energies': array([1.e-11, 1.e-05, 1.e-01, 1.e+01]),
     'fission_qvalue': 200.0,
     'capture_qvalue': 5.0,
     'reactions': {952420: {'qval': 5.0, 'xs': array([1., 2., 3., 4.])},
      18: {'qval': 200.0, 'xs': array([1., 2., 3., 4.])},
      942410: {'qval': 5.0, 'xs': array([6., 5.])},
      10010: {'qval': 0.7, 'xs': array([6., 5.])}}}
    """
    it = iter(lines)
    # READ FIRST LINE
    vals = read_line(next(it))
    zam = int(vals[0])        # ZAM identifier
    nr = int(vals[1])         # number of reactions
    ne = int(vals[2])         # number of energy points
    fiss_qval = vals[3]       # fission q-value
    # READ SECOND LINE
    capt_qval = read_line(next(it))  # radiative capture q-value
    # READ ENERGY BLOCK
    energies = read_array(it, ne)
    # READ REACTIONS BLOCK
    reactions = {}
    for ir in range(nr):
        vals = read_line(next(it))
        zap = int(vals[0])     # identifier of reaction/product
        jskip = int(vals[1])   # number of points to skip before threshold
        qval = 0.              # reaction qvalue
        if len(vals) == 3:
            qval = vals[2]
        xs = read_array(it, ne - jskip)
        if np.isnan(xs).any():
            msg = f"found 'nan' in cross section ZAM={zam} ZAP={zap}"
            logging.warning(msg)
            xs = np.nan_to_num(xs)
        reaction = {
            "qval": qval,
            "xs": xs,
            }
        reactions[zap] = reaction
    dct = {
        "nuclide": zam,
        "energies": energies,
        "fission_qvalue": fiss_qval,
        "capture_qvalue": capt_qval,
        "reactions": reactions,
        }
    return dct


def read(text):
    """

    Parameters
    ----------
    text : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    Examples
    --------
    Read test case.
    >>> read(test_aleph2xs)
    {'nuclide': 952410,
     'energies': array([1.e-11, 1.e-05, 1.e-01, 1.e+01]),
     'fission_qvalue': 200.0,
     'capture_qvalue': 5.0,
     'reactions': {952420: {'qval': 5.0, 'xs': array([1., 2., 3., 4.])},
      18: {'qval': 200.0, 'xs': array([1., 2., 3., 4.])},
      942410: {'qval': 5.0, 'xs': array([6., 5.])},
      10010: {'qval': 0.7, 'xs': array([6., 5.])}}}
    """
    lines = text.splitlines()
    return read_lines(lines)


def write_line(a="", b="", c="", d=""):
    a_ = f"{a:>.12E}" if isinstance(a, float) else a
    b_ = f"{b:>.12E}" if isinstance(b, float) else b
    c_ = f"{c:>.12E}" if isinstance(c, float) else c
    d_ = f"{d:>.12E}" if isinstance(d, float) else d
    return " ".join([f"{a_:>20}", f"{b_:>19}", f"{c_:>19}", f"{d_:>19}"])


def write_array(array, per_line=4):
    foo = sandy.utils.grouper
    lines = [
        write_line(*block) for block in foo(array, per_line, fillvalue="")
        ]
    return lines


def write_xs(idx, array, npoints, qval):
    skip_points = npoints - len(array)
    lines = [write_line(idx, skip_points, qval)]
    lines += write_array(array)
    return lines
