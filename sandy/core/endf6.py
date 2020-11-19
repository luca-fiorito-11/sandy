# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 14:50:33 2019

@author: lfiorito
"""
import pdb
import io
import os
import tempfile

import pandas as pd

import sandy

__author__ = "Luca Fiorito"
__all__ = [
        "Endf6",
        ]

pd.options.display.float_format = '{:.5e}'.format


class _FormattedFile():
    """
    Base class to store ENDF-6 content grouped by `(MAT, MF, MT)`

    Attributes
    ----------
    data

    keys

    kind : `str`
        Kind of ENDF-6 formatted file (`'endf6'`, `'pendf'`, `'gendf'`,
        `'errorr'`) .

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
        Create dataframe by reading a ENDF-6-formatted file.
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
        return {"MAT": mat, "MF": mf, "MT": mt}

    @property
    def mat(self):
        return sorted(set(self._keys["MAT"]))

    @property
    def mf(self):
        return sorted(set(self._keys["MF"]))

    @property
    def mt(self):
        return sorted(set(self._keys["MT"]))

    def to_series(self, **kwargs):
        series = pd.Series(self.data, **kwargs).sort_index(ascending=True)
        series.index.names = ["MAT", "MF", "MT"]
        return series

    @property
    def is_empty(self):
        return False if self.data else True

    @property
    def kind(self):
        """
        Kind of ENDF-6 formatted file (`'endf6'`, `'pendf'`, `'gendf'`,
        `'errorr'`) .

        Returns
        -------
        `str`
            kind of ENDF-6 formatted file

        Eaxmples
        --------
        >>> file = os.path.join(sandy.data.__path__[0], "h1.endf")
        >>> _FormattedFile.from_file(file).kind
        'endf6'

        >>> file = os.path.join(sandy.data.__path__[0], "h1.pendf")
        >>> _FormattedFile.from_file(file).kind
        'pendf'
        """
        if len(self.mat) > 1:
            msg = "Attribute 'kind' does not work if more than 1 MAT number is"
            "found"
            raise AttributeError(msg)
        mat = self.mat[0]
        text = self.data[(mat, 1, 451)]
        lrp = int(text[22:33])
        nlib = int(text[44:55])
        if nlib == -11 or nlib == -12:
            kind = "errorr"
        elif nlib == -1:
            kind = "gendf"
        else:
            if lrp == 2:
                kind = "pendf"
            elif lrp == 1 or lrp == 0:
                kind = "endf6"
            else:
                kind == "unkwown"
        return kind

    @classmethod
    def from_file(cls, file):
        """
        Create dataframe by reading a file.

        Parameters
        ----------
        file : `str`
            filename

        Returns
        -------
        `sandy.formats.endf6.BaseFile` or derived instance
            Dataframe containing ENDF6 data grouped by MAT/MF/MT

        Examples
        --------
        Read hydrogen tape from endf-6 formatted file.
        >>> file = os.path.join(sandy.data.__path__[0], "h1.endf")
        >>> _FormattedFile.from_file(file)
        MAT  MF  MT 
        125  1   451     1.001000+3 9.991673-1          0          0  ...
             2   151     1.001000+3 9.991673-1          0          0  ...
             3   1       1.001000+3 9.991673-1          0          0  ...
                 2       1.001000+3 9.991673-1          0          0  ...
                 102     1.001000+3 9.991673-1          0          0  ...
             4   2       1.001000+3 9.991673-1          0          1  ...
             6   102     1.001000+3 9.991673-1          0          2  ...
             33  1       1.001000+3 9.991673-1          0          0  ...
                 2       1.001000+3 9.991673-1          0          0  ...
                 102     1.001000+3 9.991673-1          0          0  ...
        dtype: object

        Read hydrogen tape from text stream.
        >>> stream = io.StringIO(open(file).read())
        >>> _FormattedFile.from_file(stream)
        MAT  MF  MT 
        125  1   451     1.001000+3 9.991673-1          0          0  ...
             2   151     1.001000+3 9.991673-1          0          0  ...
             3   1       1.001000+3 9.991673-1          0          0  ...
                 2       1.001000+3 9.991673-1          0          0  ...
                 102     1.001000+3 9.991673-1          0          0  ...
             4   2       1.001000+3 9.991673-1          0          1  ...
             6   102     1.001000+3 9.991673-1          0          2  ...
             33  1       1.001000+3 9.991673-1          0          0  ...
                 2       1.001000+3 9.991673-1          0          0  ...
                 102     1.001000+3 9.991673-1          0          0  ...
        dtype: object
        """
        if isinstance(file, io.StringIO):
            text = file.read()
        else:
            with open(file) as f:
                text = f.read()
        return cls.from_text(text)

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

        Examples
        --------
        Read hydrogen tape from text.
        >>> file = os.path.join(sandy.data.__path__[0], "h1.endf")
        >>> text = open(file).read()
        >>> _FormattedFile.from_text(text)
        MAT  MF  MT 
        125  1   451     1.001000+3 9.991673-1          0          0  ...
             2   151     1.001000+3 9.991673-1          0          0  ...
             3   1       1.001000+3 9.991673-1          0          0  ...
                 2       1.001000+3 9.991673-1          0          0  ...
                 102     1.001000+3 9.991673-1          0          0  ...
             4   2       1.001000+3 9.991673-1          0          1  ...
             6   102     1.001000+3 9.991673-1          0          2  ...
             33  1       1.001000+3 9.991673-1          0          0  ...
                 2       1.001000+3 9.991673-1          0          0  ...
                 102     1.001000+3 9.991673-1          0          0  ...
        dtype: object
        """
        df = pd.read_fwf(
            io.StringIO(text),
            widths=[66, 4, 2, 3],
            names=["TEXT", "MAT", "MF", "MT"],
            dtype={"TEXT": str, "MAT": int, "MF": int, "MT": int},
            na_filter=False,  # speeds up and does not add NaN in empty lines
            # Do not use TEXT because  the parser does not preserve the
            # whitespaces
            usecols=("MAT", "MF", "MT"),
            )
        # use splitlines instead of readlines to remove "\n"
        df["TEXT"] = text.splitlines()
        condition = (df.MT > 0) & (df.MF > 0) & (df.MAT > 0)
        data = df[condition].groupby(["MAT", "MF", "MT"])\
                            .agg({"TEXT": "\n".join})\
                            .TEXT\
                            .to_dict()
        return cls(data)

    def _get_section_df(self, mat, mf, mt, delimiter="?"):
        """
        """
        text = self.data[(mat, mf, mt)]
        foo = lambda x : sandy.shared.add_delimiter_every_n_characters(x[:66], 11, delimiter=delimiter)
        newtext = "\n".join(map(foo, text.splitlines()))
        df = pd.read_csv(
            io.StringIO(sandy.shared.add_exp_in_endf6_text(newtext)),
            delimiter=delimiter,
            na_filter=True,
            names=["C1", "C2", "L1", "L2", "N1", "N2"],
        )
        return df

    def add_section(self, mat, mf, mt, text, inplace=False):
        """
        Given MAT, MF and MT add/replace the corresponding section in the
        `Endf6.data`.

        Parameters
        ----------
        mat : `int`
            MAT number
        mf : `int`
            MF number
        mt : `int`
            MT number
        inplace : `bool`, optional, default is `False`
            flag to operate **inplace**

        Returns
        -------
        `sandy._FormattedFile` or derived instance
            object with new section

        Examples
        --------
        >>> tape = sandy.Endf6({(9437, 3, 102) : "lorem ipsum"})
        >>> tape.add_section(9999, 1, 1, "dolor sit amet")
        MAT   MF  MT 
        9437  3   102       lorem ipsum
        9999  1   1      dolor sit amet
        dtype: object

        >>> tape.add_section(9437, 3, 102, "new text", inplace=True)
        >>> tape
        MAT   MF  MT 
        9437  3   102    new text
        dtype: object
        """
        d = self.data.copy()
        key = (mat, mf, mt)
        d[key] = text
        if inplace:
            self.data = d
        else:
            return self.__class__(d)

    def add_sections(self, sections, inplace=False):
        d = self.data.copy()
        for (mat, mf, mt), text in sections.items():
            key = (mat, mf, mt)
            d[key] = text
        if inplace:
            self.data = d
        else:
            return self.__class__(d)

    def delete_section(self, mat, mf, mt, inplace=False, raise_error=True):
        """
        Given MAT, MF and MT delete the corresponding section from the
        `Endf6.data`.

        Parameters
        ----------
        mat : `int`
            MAT number
        mf : `int`
            MF number
        mt : `int`
            MT number
        inplace : `bool`, optional, default is `False`
            flag to operate **inplace**

        Returns
        -------
        `sandy._FormattedFile` or derived instance
            object without given section

        Examples
        --------
        Delete capture cross section from hydrogen file.
        >>> file = os.path.join(sandy.data.__path__[0], "h1.endf")
        >>> tape = _FormattedFile.from_file(file)
        >>> new = tape.delete_section(125, 3, 102)
        >>> new
        MAT  MF  MT 
        125  1   451     1.001000+3 9.991673-1          0          0  ...
             2   151     1.001000+3 9.991673-1          0          0  ...
             3   1       1.001000+3 9.991673-1          0          0  ...
                 2       1.001000+3 9.991673-1          0          0  ...
             4   2       1.001000+3 9.991673-1          0          1  ...
             6   102     1.001000+3 9.991673-1          0          2  ...
             33  1       1.001000+3 9.991673-1          0          0  ...
                 2       1.001000+3 9.991673-1          0          0  ...
                 102     1.001000+3 9.991673-1          0          0  ...
        dtype: object

        This method can also work **inplace**.
        >>> tape.delete_section(125, 3, 102, inplace=True)
        >>> assert tape.data == new.data
        """
        d = self.data.copy()
        key = (mat, mf, mt)
        if key not in d and raise_error is False:
            pass
        else:
            del d[key]
        if inplace:
            self.data = d
        else:
            return self.__class__(d)

    def delete_sections(self, sections, inplace=False, raise_error=True):
        d = self.data.copy()
        for mat, mf, mt in sections:
            key = (mat, mf, mt)
            if key not in d and raise_error is False:
                pass
            else:
                del d[key]
        if inplace:
            self.data = d
        else:
            return self.__class__(d)

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

    def get_value(self, mat, mf, mt, line_number, pos):
        return self._get_section_df(mat, mf, mt)[pos.upper()] \
                   .iloc[line_number - 1]

    def change_value(self, val, mat, mf, mt, line_number, pos, inplace=False,
                     dtype=float):
        items = ("C1", "C2", "L1", "L2", "N1", "N2")
        positions = {y: x for x, y in enumerate(items)}
        step = positions[pos]
        length = 81
        ibeg = length * (line_number - 1) + 11 * step
        iend = length * (line_number - 1) + 11 * (step + 1)
        text = self.data[(mat, mf, mt)]
        new = sandy.write_int(val) if dtype is int else sandy.write_float(val)
        new_text = "".join((text[:ibeg], new, text[iend:]))
        new_tape = self.add_section(mat, mf, mt, new_text, inplace=False)
        if inplace:
            self.data = new_tape.data
        else:
            return new_tape
        print(new_text)


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
        text = self.data[(mat, mf, mt)]
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
        foo = eval(f"sandy.read_mf{mf}")
        return foo(self, mat, mt)

    def write_string(self, title="", skip_title=False, skip_fend=False):
        """
        .. warning:: OUTDATED, use `to_file`.
        """
        string = sandy.write_line(title, 1, 0, 0, 0)
        string += "\n"
        for mat, dfmat in self.to_series().groupby('MAT', sort=True):
            for mf, dfmf in dfmat.groupby('MF', sort=True):
                for mt, text in dfmf.groupby('MT', sort=True):
                    string += text.squeeze()\
                                  .encode('ascii', 'replace')\
                                  .decode('ascii')
                    string += "\n"
                    string += sandy.write_line("", mat, mf, 0, 99999)
                    string += "\n"
                string += sandy.write_line("", mat, 0, 0, 0)
                string += "\n"
            string += sandy.write_line("", 0, 0, 0, 0)
            string += "\n"
        string += sandy.write_line("", -1, 0, 0, 0)
        return string

    def to_string(self, title="", skip_title=False, skip_fend=False):
        """
        Write `Endf6.data` content to string according to the ENDF-6 file
        rules.

        Parameters
        ----------
        title : `str`, optional, default is an empty string
            first line of the file

        Returns
        -------
        `str`
            string containing the ENDF-6 information stored in this instance.

        Notes
        -----
        ..note:: no modification is implemented to the actual content of
                 the `Endf6.data` object.

        Examples
        --------
        >>> file = os.path.join(sandy.data.__path__[0], "h1.endf")
        >>> string = sandy.Endf6.from_file(file).write_string()
        >>> print(string[:81 * 4 - 1])
                                                                             1 0  0    0
         1.001000+3 9.991673-1          0          0          2          5 125 1451    1
         0.000000+0 0.000000+0          0          0          0          6 125 1451    2
         1.000000+0 2.000000+7          3          0         10          3 125 1451    3

        if no modification is applied to the `Endf6` content, the
        `write_string` returns an output identical to the file ASCII content.
        >>> assert string == open(file).read()

        ..note:: differences might appear from the way zeros were handled at
                 the end of ENDF-6 section, or if a different fiel title is
                 given
        """
        string = sandy.write_line(title, 1, 0, 0, 0)
        string += "\n"
        for mat, dfmat in self.to_series().groupby('MAT', sort=True):
            for mf, dfmf in dfmat.groupby('MF', sort=True):
                for mt, text in dfmf.groupby('MT', sort=True):
                    string += text.squeeze()\
                                  .encode('ascii', 'replace')\
                                  .decode('ascii')
                    string += "\n"
                    string += sandy.write_line("", mat, mf, 0, 99999)
                    string += "\n"
                string += sandy.write_line("", mat, 0, 0, 0)
                string += "\n"
            string += sandy.write_line("", 0, 0, 0, 0)
            string += "\n"
        string += sandy.write_line("", -1, 0, 0, 0)
        return string

    def to_file(self, filename, mode="w", **kwargs):
        text = self.write_string(**kwargs)
        with open(filename, mode) as f:
            f.write(text)

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

    def get_pendf(self, njoy, debug=False):
        """
        #>>> file = os.path.join(sandy.data.__path__[0], "h1.endf")
        #>>> endf6 = sandy.Endf6.from_file(file)
        #>>> pendf = endf6.get_pendf("/home/lfiorito/vc/NJOY2016/bin/njoy")
        #>>> pendf
        MAT  MF  MT 
        125  1   451     1.001000+3 9.991673-1          2          0  ...
             2   151     1.001000+3 9.991673-1          0          0  ...
             3   1       1.001000+3 9.991673-1          0         99  ...
                 2       1.001000+3 9.991673-1          0          0  ...
                 102     1.001000+3 9.991673-1          0          0  ...
        dtype: object

        #>>> pendf.kind
        'pendf'
        """
        with tempfile.TemporaryDirectory() as td:
            endf6file = os.path.join(td, "endf6_file")
            self.to_file(endf6file)
            outputs = sandy.njoy.process(
                    endf6file,
                    broadr=False,
                    thermr=False,
                    unresr=False,
                    heatr=False,
                    gaspr=False,
                    purr=False,
                    errorr=False,
                    acer=False,
                    wdir=".",
                    keep_pendf=True,
                    exe=njoy,
                    temperatures=[0],
                    suffixes=[0],
                    err=0.005
                    )[2]  # keep only pendf filename
            pendffile = outputs["tape30"]
            pendf = Endf6.from_file(pendffile)
        return pendf

    def merge_pendf(self, pendf):
        """
        Merge endf-6 file content with that of a pendf file.

        Parameters
        ----------
        pendf : `sandy.Endf6`
            `Endf6` object containing pendf tape.

        Returns
        -------
        `sandy.Endf6`
            `Endf6` object with MF3 and MF1MT451 from PENDF

        Notes
        -----
        .. note:: the `Endf6` object in output has attribute `kind='pendf'`.

        .. note:: the `Endf6` object in output contains all sections from the
                  original endf-6 file, but for all MF=3 and MF=1 MT=451.

        Examples
        --------
        >>> file = os.path.join(sandy.data.__path__[0], "h1.endf")
        >>> endf6 = sandy.Endf6.from_file(file)
        >>> file = os.path.join(sandy.data.__path__[0], "h1.pendf")
        >>> pendf = sandy.Endf6.from_file(file)
        >>> merged = endf6.merge_pendf(pendf)
        >>> merged
        MAT  MF  MT 
        125  1   451     1.001000+3 9.991673-1          2          0  ...
             2   151     1.001000+3 9.991673-1          0          0  ...
             3   1       1.001000+3 9.991673-1          0         99  ...
                 2       1.001000+3 9.991673-1          0          0  ...
                 102     1.001000+3 9.991673-1          0          0  ...
             4   2       1.001000+3 9.991673-1          0          1  ...
             6   102     1.001000+3 9.991673-1          0          2  ...
             33  1       1.001000+3 9.991673-1          0          0  ...
                 2       1.001000+3 9.991673-1          0          0  ...
                 102     1.001000+3 9.991673-1          0          0  ...
        dtype: object

        The cross section in the merged file come from the pendf.
        >>> assert merged.data[125, 3, 1] == pendf.data[125, 3, 1]
        >>> assert merged.data[125, 3, 1] != endf6.data[125, 3, 1]

        The new file is also a pendf.
        >>> merged.kind
        'pendf'
        """
        if pendf.kind != "pendf":
            raise sandy.Error("given file is not a pendf")
        section_3_endf6 = self.filter_by(listmf=[3]).data
        section_3_pendf = pendf.filter_by(listmf=[3]).data
        section_1451_pendf = pendf.filter_by(listmf=[1], listmt=[451]).data
        return self.delete_sections(section_3_endf6) \
                   .add_sections(section_3_pendf) \
                   .add_sections(section_1451_pendf)

    def get_errorr(self, njoy, pendf=None):
        """
        .. important:: this module uses `read_formatted_file`, which is no
                       longer supported.

        Examples
        --------
        #>>> file = os.path.join(sandy.data.__path__[0], "h1.endf")
        #>>> endf6 = sandy.Endf6.from_file(file)
        #>>> endf6.get_errorr("/home/lfiorito/vc/NJOY2016/bin/njoy")
                                                                 TEXT
        MAT MF MT 
        125 1  451   1.001000+3 9.991673-1          6          0  ...
            3  1     1.001000+3 0.000000+0          0          0  ...
               2     1.001000+3 0.000000+0          0          0  ...
               102   1.001000+3 0.000000+0          0          0  ...
            33 1     1.001000+3 9.991673-1          0          0  ...
               2     1.001000+3 9.991673-1          0          0  ...
               102   1.001000+3 9.991673-1          0          0  ...

        #>>> file = os.path.join(sandy.data.__path__[0], "h1.pendf")
        #>>> pendf = sandy.Endf6.from_file(file)
        #>>> pendf.get_errorr("/home/lfiorito/vc/NJOY2016/bin/njoy")

        #>>>endf6.get_errorr("/home/lfiorito/vc/NJOY2016/bin/njoy", pendf=pendf)
        """
        if self.kind != "endf6":
            raise sandy.Error("stop because `kind!='endf6'`")
        with tempfile.TemporaryDirectory() as td:
            endf6file = os.path.join(td, "endf6_file")
            self.to_file(endf6file)
            if pendf:
                pendftape = os.path.join(td, "pendf_file")
                pendf.to_file(pendftape)
            else:
                pendftape = None
            outputs = sandy.njoy.process(
                    endf6file,
                    pendftape=pendftape,
                    broadr=False,
                    thermr=False,
                    unresr=False,
                    heatr=False,
                    gaspr=False,
                    purr=False,
                    errorr=True,
                    acer=False,
                    wdir=".",
                    keep_pendf=False,
                    exe=njoy,
                    temperatures=[0],
                    suffixes=[0],
                    err=0.005
                    )[2]  # keep only pendf filename
            errorr = sandy.read_formatted_file(outputs["tape33"])
        return errorr
