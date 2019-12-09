# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 14:50:33 2019

@author: lfiorito
"""
import pdb
import logging
import io

import pandas as pd

import sandy

__author__ = "Luca Fiorito"
__all__ = [
        "Endf6",
        ]

class _FormattedFile():
    """
    Base class to store ENDF-6 content grouped by `(MAT, MF, MT)`
    
    Attributes
    ----------
    data
    
    keys
    
    mat : `int`
        MAT number.
    mf : `int`
        MF number.
    mt : `int`
        MT number

    Methods
    -------
    add_sections
        Add text section for given `(MAT, MF, MT)`.
    filter_by
        Filter dataframe based on `(MAT, MF, MT)` lists.
    from_file
        Create dataframe by reading a ENDF-6 file.
    from_text
        Create dataframe from endf6 text in string.
    to_series
        Covert content into `pandas.Series`.
    
    Notes
    -----
    This class supports ENDF-6 content from ENDF-6 files, ERRORR files and 
    GROUPR files.
    """
    def __repr__(self):
        return self.to_series().__repr__()
        
    def __init__(self, data, file=None):
        self.data = data
        self.file = file

    @property
    def data(self):
        """
        Dataframe of energy-dependent tabulated cross sections.
        
        Attributes
        ----------
        index : `pandas.Index`
            energy grid in eV
        columns : `pandas.MultiIndex`
            MAT/MT indices
        values : `numpy.array`
            cross sections in barns
        
        Returns
        -------
        `pandas.DataFrame`
            tabulated xs
        
        Raises
        ------
        `sandy.Error`
            if energy grid is not monotonically increasing
        """
        return self._data
    
    @data.setter
    def data(self, data):
        if not isinstance(data, dict):
            raise sandy.Error("'data' is not a 'dict'")
        self._data = data

    @property
    def keys(self):
        return self.data.keys()

    @property
    def _keys(self):
        mat, mf, mt = zip(*self.data.keys())
        return {"MAT" : mat, "MF" : mf, "MT" : mt}

    @property
    def mat(self):
        return sorted(set(self._keys["MAT"]))

    @property
    def mf(self):
        return sorted(set(self._keys["MF"]))

    @property
    def mt(self):
        return sorted(set(self._keys["MT"]))
    
    def to_series(self):
        series = pd.Series(self.data, name=self.file).sort_index(ascending=True)
        series.index.names = ["MAT", "MF", "MT"]
        return series
        
    @classmethod
    def from_file(cls, file):
        """Create dataframe by reading a file.
        
        Parameters
        ----------
        file : `str`
            file name

        Returns
        -------
        `sandy.formats.endf6.BaseFile` or derived instance
            Dataframe containing ENDF6 data grouped by MAT/MF/MT
        """
        with open(file) as f:
            text = f.read()
        out = cls.from_text(text)
        out.file = file
        return out

    @classmethod
    def from_text(cls, text):
        """Create dataframe from endf6 text in string.
        
        Parameters
        ----------
        text : `str`
            string containing the evaluated data

        Returns
        -------
        `sandy.formats.endf6.BaseFile` or derived instance
            Dataframe containing ENDF6 data grouped by MAT/MF/MT
        """
        df = pd.read_fwf(
            io.StringIO(text),
            widths=[66, 4, 2, 3],
            names=["TEXT", "MAT", "MF", "MT"],
            dtype={"TEXT" : str, "MAT" : int, "MF" : int, "MT" : int},
            na_filter=False, # speeds up and does not add NaN in empty lines
            usecols = ("MAT", "MF", "MT") # Do not use TEXT beacuse  the parser does not preserve the whitespaces
        )
        df["TEXT"] = text.splitlines() # use splitlines instead of readlines to remove ""\n"
        data = df[(df.MT>0) & (df.MF>0) & (df.MAT>0)].groupby(["MAT", "MF", "MT"]).agg({"TEXT" : "\n".join}).TEXT
        return cls(data.to_dict())

    def add_section(self, mat, mf, mt, text, inplace=False):
        """Collapse two tapes into a single one.
        If MAT/MF/MT index is present in both tapes, take it from the second.
        
        Parameters
        ----------
        tape : `sandy.formats.endf6.BaseFile` or derived instance
            dataframe for ENDF-6 formatted file
        
        Returns
        -------
        `sandy.formats.endf6.BaseFile` or derived instance
            dataframe with merged content
        """
        key = (mat,mf,mt)
        if inplace:
            self.data[key] = text
        else:
            data = self.data.copy()
            data[key] = text
            return self.__class__(data, file=self.file)

    def filter_by(self, listmat=range(1,10000), listmf=range(1,10000), listmt=range(1,10000), inplace=False):
        """Filter dataframe based on MAT, MF, MT lists.
        
        Parameters
        ----------
        listmat : `list` or `None`
            list of requested MAT values (default is `None`: use all MAT)
        listmf : `list` or `None`
            list of requested MF values (default is `None`: use all MF)
        listmt : `list` or `None`
            list of requested MT values (default is `None`: use all MT)
        
        Returns
        -------
        `sandy.formats.endf6.BaseFile` or derived instance
            Copy of the original instance with filtered MAT, MF and MT sections
        """
        series = self.to_series()
        cond_mat = series.index.get_level_values("MAT").isin(listmat)
        cond_mf = series.index.get_level_values("MF").isin(listmf)
        cond_mt = series.index.get_level_values("MT").isin(listmt)
        d = series.loc[cond_mat & cond_mf & cond_mt].to_dict()
        if inplace:
            self.data = d
        else:
            return self.__class__(d, file=self.file)
    
    def _get_file_format(self):
        """Determine ENDF-6 format type by reading flags "NLIB" and "LRP" of first MAT in file:
            
            * `NLIB = -11 | NLIB = -12` : errorr
            * `NLIB = -1` : gendf
            * `LRP = 2` : pendf
            * `LRP != 2` : endf6
        
        Returns
        -------
        `str`
            type of ENDF-6 format
        """
        lines = self.TEXT.loc[self.mat[0], 1, 451].splitlines()
        C, i = read_cont(lines, 0)
        if C.N1 == -11 or C.N1 == -12:
            ftype = "errorr"
        elif C.N1 == -1:
            ftype = "gendf"
        else:
            if C.L1 == 2:
                ftype = "pendf"
            else:
                ftype = "endf6"
        return ftype



class Endf6(_FormattedFile):
    """
    Container for ENDF-6 file text grouped by MAT, MF and MT numbers.
    
    Methods
    -------
    read_section
        Parse MAT/MF/MT section.
    write_string
        Write ENDF-6 content to string.
    """

    def _get_nsub(self):
        """
        Determine ENDF-6 sub-library type by reading flag "NSUB" of first MAT in file:
            
            * `NSUB = 10` : Incident-Neutron Data
            * `NSUB = 11` : Neutron-Induced Fission Product Yields
        
        Returns
        -------
        `int`
            NSUB value
        """
        return self.read_section(self.mat[0], 1, 451)["NSUB"]

    def _get_section_df(self, mat, mf, mt, delimiter="?"):
        """
        """
        text = self.data[(mat,mf,mt)]
        foo = lambda x : sandy.shared.add_delimiter_every_n_characters(x[:66], 11, delimiter=delimiter)
        newtext = "\n".join(map(foo, text.splitlines()))
        df = pd.read_csv(
            io.StringIO(sandy.shared.add_exp_in_endf6_text(newtext)),
            delimiter=delimiter,
            na_filter=True,
            names=["C1", "C2", "L1", "L2", "N1", "N2"],
        )
        return df
    
    def read_section(self, mat, mf, mt):
        """
        Parse MAT/MF/MT section.
        
        
        Parameters
        ----------
        `mat` : int
            MAT number
        `mf` : int
            MF number
        `mt` : int
            MT number

        Returns
        -------
        `dict`
        """
        foo = eval("sandy.read_mf{}".format(mf))
        return foo(self, mat, mt)

    def write_string(self, title="", skip_title=False, skip_fend=False):
        """
        Write ENDF-6 content to string.
        
        Parameters
        ----------
        title : `str`, optional, default is an empty string
            first line of the file

        Returns
        -------
        `str`
            string containing the ENDF-6 information stored in this instance.
        """
        string = sandy.write_line(title, 1, 0, 0, 0) + "\n"
        for mat,dfmat in self.to_series().groupby('MAT', sort=True):
            for mf,dfmf in dfmat.groupby('MF', sort=True):
                for mt,text in dfmf.groupby('MT', sort=True):
                    string += text.squeeze().encode('ascii', 'replace').decode('ascii')  + "\n"
                    string += sandy.write_line("", mat, mf, 0, 99999) + "\n"
                string += sandy.write_line("", mat, 0, 0, 0) + "\n"
            string += sandy.write_line("", 0, 0, 0, 0) + "\n"
        string += sandy.write_line("", -1, 0, 0, 0)
        return string

    def _update_info(self, descr=None):
        """Update RECORDS item (in DATA column) for MF1/MT451 of each MAT based on the content of the TEXT column.
        """
        from .mf1 import write
        tape = self.copy()
        for mat in sorted(tape.index.get_level_values('MAT').unique()):
            sec = self.read_section(mat,1,451)
            records = pd.DataFrame(sec["RECORDS"], columns=["MF","MT","NC","MOD"]).set_index(["MF","MT"])
            new_records = []
            dfmat=tape.loc[mat]
#            for (mf,mt),text in sorted(tape.loc[mat].query('MT!=451'.format(mat)).TEXT.items()):
            for (mf,mt),text in sorted(dfmat[dfmat.index.get_level_values("MT")!=451].TEXT.items()):
                nc = len(text.splitlines())
                # when copying PENDF sections (MF2/MT152) mod is not present in the dictionary
                try:
                    mod = records.MOD.loc[mf,mt]
                except:
                    mod = 0
                new_records.append((mf,mt,nc,mod))
            if descr is not None:
                sec["TEXT"] = descr
            nc = 4 + len(sec["TEXT"]) + len(new_records) + 1
            mod = records.MOD.loc[1,451]
            new_records = [(1,451,nc,mod)] + new_records
            sec["RECORDS"] = new_records
            text = write(sec)
            tape.loc[mat,1,451].TEXT = text
        return Endf6(tape)

    def parse(self):
        mats = self.index.get_level_values("MAT").unique()
        if len(mats) > 1:
            raise NotImplementedError("file contains more than 1 MAT")
        self.mat = self.endf = mats[0]
        if hasattr(self, "tape"):
            self.filename = os.path.basename(self.tape)
        INFO = self.read_section(mats[0], 1 ,451)
        del INFO["TEXT"], INFO["RECORDS"]
        self.__dict__.update(**INFO)
        self.SECTIONS = self.loc[INFO["MAT"]].reset_index()["MF"].unique()
        self.EHRES = 0
        self.THNMAX = - self.EHRES if self.EHRES != 0 else 1.0E6

