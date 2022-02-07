"""
This module contains all classes and functions necessary to work with
the ALEPH-2 output file.
"""

__author__ = "ALEPH development team"
__all__ = [
    "OutputFile",
    "read_output",
    ]

import re
import io
import logging
import itertools
from collections import namedtuple

import numpy as np
import pandas as pd

from sandy.utils import grouper


summary_header = "\n".join((
    "\*" * 80,
    "\*                            ALEPH problem summary                             \*",
    "\*" * 80,
))

table_header = "\n".join((
    "\*" * 80,
    "\*\-\-\- Table (?P<table_number>[0-9\s]{2}) \-\-\- (?P<table_name>.*?)\*",
    "\*" * 80,
))
material_header = "\*\*\*\*\*   Material (?P<material_number>[0-9\s]{7})\n"
element_header = "\*\*\*\*\*   Per element.*\n"
contributor_header = "Main contributors.*\n\n\n"

table_footer = "^\s+Total"

PATTERN_TIMEKEFF = re.compile("\s+Global neutronics parameters\n\s+\-+\n\s+Time\s+\(days\)\s+(?P<data>.*?)\n")
PATTERN_BURNUP = re.compile("^\sFuel burnup \(MWd/kg HM\)\s+(?P<data>.*?)$", flags=re.MULTILINE)
PATTERN_KEFF = re.compile("^\s+Keff  eff. mult. factor\s+(?P<data>.*?)$", flags=re.MULTILINE)
PATTERN_DKEFF = re.compile("^\s+Relative std. deviation\s+(?P<data>.*?)$", flags=re.MULTILINE)

PATTERN_BEGIN = re.compile("\s+Irradiated materials\n\s+\-{20}")
PATTERN_END = re.compile("\s+\*{28}\n\s+\* Total over all materials \*\n\s+\*{28}")

PATTERN_MATERIAL = re.compile("^\*{5}\s+Material\s+(?P<mat>.*?)$", flags=re.MULTILINE)
PATTERN_CELLS = re.compile("\n\n\s+Cells\s+=\s+(?P<data>(?:[0-9, ]+\n)+)\n")
PATTERN_VOLUME = re.compile("^\s+Volumes \(cm3\)\s+=\s+(?P<data>.*?) $", flags=re.MULTILINE)

PATTERN_TIME = re.compile("^\s+Time\s+\(days\)\s+(?P<data>.*?)$", flags=re.MULTILINE)
PATTERN_WDENSITY = re.compile("^\s+Density\s+\(g/cm3\)\s+(?P<data>.*?)$", flags=re.MULTILINE)
PATTERN_ADENSITY = re.compile("^\s+Density\s+\(at/\(b\*cm\)\)\s+(?P<data>.*?)$", flags=re.MULTILINE)
PATTERN_TEMPERATURE = re.compile("^\s+Temperature\s+\(K\)\s+(?P<data>.*?)$", flags=re.MULTILINE)
PATTERN_SOURCE = re.compile("^\s+Source strength\s+\(part/s\)\s+(?P<data>.*?)$", flags=re.MULTILINE)
PATTERN_NFLUX = re.compile("^\s+Neutron flux\s+\(n/\(cm2\*s\)\)\s+(?P<data>.*?)$", flags=re.MULTILINE)
PATTERN_POWER = re.compile("^\s+Thermal power\s+\(MW\)\s+(?P<data>.*?)$", flags=re.MULTILINE)
PATTERN_BURNUP = re.compile("^\s+Fuel burnup\s+\(MWd/kg HM\)\s+(?P<data>.*?)$", flags=re.MULTILINE)
PATTERN_TOTBURNUP = re.compile("^\s+Total burnup \(MWd/kg HM\):\s+(?P<data>.*?)$", flags=re.MULTILINE)
PATTERN_FISSIONS = re.compile("^\s+Integral number of fissions\s+(?P<data>.*?)$", flags=re.MULTILINE)
PATTERN_TOTFISSIONS = re.compile("^\s+Total number of fissions:\s+(?P<data>.*?)$", flags=re.MULTILINE)

MATERIAL_KEYS = [
    "ID",
    "cells",
    "volume",
    "time",
    "weight_density",
    "atomic_density",
    "temperature",
    "source_strenght",
    "neutron_flux",
    "power",
]

FISS_MATERIAL_KEYS = MATERIAL_KEYS + [
    "burnup",
    "cumulative_burnup",
    "total_burnup",
    "fissions",
    "total_fissions",
]

Mat = namedtuple("Material", MATERIAL_KEYS)
FissMat = namedtuple("FissileMaterial", FISS_MATERIAL_KEYS)


class OutputFile():
    """
    Class dedicated to store the content of the `'output'` ile produced by
    in any succesfull ALEPH-2 run.

    Attributes
    ----------
    data : `dict`
        dictionary of ALEPH-2 output sections
    burnup : `pandas.Series`
        series with time-dependent burnup
    cumulative_burnup : `pandas.Sereis`
        series with time-dependent cumulative burnup
    keff : `pandas. DataFrame`
        dataframe with time-dependent keff and associated statistical
        error
    summary : `str`
        summary of the ALEPH-2 output (first output section before the tables)
    text : `str`
        full ASCII text of ALEPH-2 output

    Methods
    -------
    from_file
        initialize object reading ALEPH-2 output from file
    from_string
        initialize object reading ALEPH-2 output from string
    get_burnup
        get array of burnup values
    get_keff
        get array of keff values
    get_keff_staterr
        get array of keff statistical errors
    get_table
        for a specified table number, return a dictionary of tables
        indexed by material number
    get_time
        get array of irradiation/decay time steps (in days)
    parse_tables
    """

    def __init__(self, text):
        self.text = text

    def parse_output(self, table_index="ZAM"):
        self.data = parse_output_for_tables(self.text, index=table_index)

    @property
    def summary(self):
        """
        Returns
        -------
        `str`
            summary of the ALEPH-2 output (first output section before the
                                           tables).
        """
        return self.data[0]

    def get_burnup(self):
        """
        Returns
        -------
        `pandas.Series`
            series with time-dependent burnup.
        """
        mats = self.get_materials()
        bu = pd.DataFrame({
            m: mat.burnup for m, mat in mats.items()
            if hasattr(mat, "burnup")
            })
        days = self.get_time()
        bu.index = pd.Index(days, name="days")
        bu.columns.name = "material"
        return bu

    def get_cumulative_burnup(self):
        """
        Returns
        -------
        `pandas.Series`
            series with time-dependent cumulative burnup.
        """
        mats = self.get_materials()
        bu = pd.DataFrame({
            m: mat.cumulative_burnup for m, mat in mats.items()
            if hasattr(mat, "cumulative_burnup")
            })
        days = self.get_time()
        bu.index = pd.Index(days, name="days")
        bu.columns.name = "material"
        return bu

    @property
    def keff(self):
        """
        Returns
        -------
        `pandas.DataFrame`
        dataframe with time-dependent keff and associated statistical
        error
        """
        keff = self.get_keff()
        dkeff = self.get_keff_staterr()
        days = self.get_time()
        index = pd.Index(days, name="days")
        return pd.DataFrame({
            "keff": keff,
            "staterr": dkeff,
            }, index=index)

    @property
    def data(self):
        """
        Dictionary of ALEPH-2 output sections.
        Keys are the table numbers, values are dictionaries with
        `pandas.DataFrame` objects indexed by material number.

        Raises
        ------
        AttributeError
            if tables were not parsed you must first run `parse_output`.

        Notes
        -----
        `data[0]` is the ALEPH-2 summary also accessible with attribute
        `summary`.

        Returns
        -------
        `dict`
            dictionary of ALEPH-2 output sections.
        """
        if not hasattr(self, "_data"):
            msg = "tables were not parsed, run first 'parse_tables'"
            raise AttributeError(msg)
        return self._data

    @data.setter
    def data(self, data):
        self._data = data

    def get_materials(self):
        if not hasattr(self, "_data"):
            self.parse_output()
        text = self.summary
        text = PATTERN_BEGIN.split(text, maxsplit=1)[1]
        text = PATTERN_END.split(text, maxsplit=1)[0]
        return parse_materials_output(text)

    def get_time(self):
        """
        Get array of irradiation/decay time steps (in days).

        Returns
        -------
        `numpy.ndarray` of `float` values
            array of time values.
        """
        match = PATTERN_TIMEKEFF.search(self.summary)
        data = match.group("data")
        return np.array([float(x) for x in data.split()])

    def get_keff(self):
        """
        Get array of keff values.

        Returns
        -------
        `numpy.ndarray` of `float` values
            array of keff values.
        """
        match = PATTERN_KEFF.search(self.summary)
        data = match.group("data")
        return np.array([float(x) for x in data.split()])

    def get_keff_staterr(self):
        """
        Get array of keff statistical errors.

        Returns
        -------
        `numpy.ndarray` of `float` values
            array of keff statistical errors.
        """
        match = PATTERN_DKEFF.search(self.summary)
        data = match.group("data")
        return np.array([float(x) for x in data.split()])

    def parse_tables(self):
        data = {}
        inp, rest = re.split(summary_header, self.text, maxsplit=1)
        tab_items = re.split(table_header, rest)
        summary, tab_items = tab_items[0], tab_items[1:]
        tab_numbers = map(int, tab_items[0::3])
        tab_titles = map(lambda x: x.strip(), tab_items[1::3])
        tab_texts = tab_items[2::3]
        tab_zipped = zip(tab_numbers, tab_titles, tab_texts)
        for num, title, tab_txt in tab_zipped:
            if num > 8 and num != 11:
                continue
            data[num] = Table(num, title, tab_txt)
        self.data = data

    def get_table(self, table_number):
        """
        For a specified table number, return a dictionary of tables
        indexed by material number.

        Parameters
        ----------
        table_number : `ïnt`
            Table number found in the ALEPH-2 output file.

        Returns
        -------
        `dict`
            dictionary of tables.
        """
        return self.data[table_number]["materials"]

    @classmethod
    def from_string(cls, string):
        """
        Initialize object reading ALEPH-2 output from string.

        Parameters
        ----------
        string: `str`
            string containing the ALEPH-2 output

        Returns
        -------
        `OutputFile`
            Container for ALEPH-2 output file content
        """
        return cls(string)

    @classmethod
    def from_file(cls, file, encoding="utf8", errors='ignore'):
        """
        Initialize object reading ALEPH-2 output from file.

        Parameters
        ----------
        file : `str`
            filename
        encoding : `str`, optional, default is `'utf8'`
            encoding used to open the file and passed to python function `open`
        errors : `str`, optional, default is `'ignore'`
            option passed to python function `open` and that specifies how
            encoding and decoding errors are to be handled

        Returns
        -------
        `OutputFile`
            Container for ALEPH-2 output file content
        """
        with open(file, encoding=encoding, errors=errors) as f:
            string = f.read()
        return cls(string)


class Table():

    def __repr__(self):
        return f"ALEPH-2 table: {self.title}"

    def __init__(self, num, title, text):
        self.number = num
        self.title = title
        self.text = text
        self.materials = self._get_tables()

    def _get_tables(self):
        out = {}
        mat_items = re.split(material_header, self.text)[1:]
        mat_number = map(int, mat_items[0::2])
        mat_texts = mat_items[1::2]
        mat_zipped = zip(mat_number, mat_texts)
        for mat, mat_txt in mat_zipped:
            if self.number == 1 or self.number == 2:
                SubTable = namedtuple(
                    f"Table{self.number}",
                    "by_nuclide by_element",
                    )
                content, rest = re.split(
                    table_footer,
                    mat_txt,
                    maxsplit=1,
                    flags=re.MULTILINE,
                    )
                df1 = parse_table_nuclide(content, index="ZAM")
                content = re.split(
                    element_header,
                    rest,
                    maxsplit=1,
                    flags=re.MULTILINE,
                    )[1]
                df2 = parse_table_element(content)
                df = SubTable(df1, df2)
            elif self.number == 3 or self.number == 4:
                SubTable = namedtuple(
                    f"Table{self.number}",
                    "by_nuclide abs_contributors rel_contributors",
                    )
                content, rest = re.split(
                    table_footer,
                    mat_txt,
                    maxsplit=1,
                    flags=re.MULTILINE,
                    )
                df1 = parse_table_nuclide(content, index="ZAM")
                absol, rel = re.split(
                    contributor_header,
                    rest,
                    maxsplit=2,
                    flags=re.MULTILINE,
                    )[1:]
                df2 = parse_table_contributor(absol)
                df3 = parse_table_contributor(rel)
                df = SubTable(df1, df2, df3)
            elif self.number == 5 or self.number == 6 or self.number == 7:
                SubTable = namedtuple(
                    f"Table{self.number}",
                    "by_nuclide by_energy count",
                    )
                content, rest = re.split(
                    table_footer,
                    mat_txt,
                    maxsplit=1,
                    flags=re.MULTILINE,
                    )
                df1 = parse_table_nuclide(content, index="ZAM")
                rest = "\n".join(rest.splitlines()[5:])
                content, rest = re.split(
                    table_footer,
                    rest,
                    maxsplit=1,
                    flags=re.MULTILINE,
                    )
                df2 = parse_table_energy(content)
                rest = "\n".join(rest.splitlines()[6:])
                content, rest = re.split(
                    table_footer,
                    rest,
                    maxsplit=1,
                    flags=re.MULTILINE,
                    )
                df3 = parse_table_nuclide(content, index="ZAM")
                df = SubTable(df1, df2, df3)
            elif self.number == 8:
                SubTable = namedtuple(
                    f"Table{self.number}",
                    "by_energy by_nuclide",
                    )
                content, rest = re.split(
                    table_footer,
                    mat_txt,
                    maxsplit=1,
                    flags=re.MULTILINE,
                    )
                df1 = parse_table_energy(content)
                rest = "\n".join(rest.splitlines()[6:])
                content, rest = re.split(
                    table_footer,
                    rest,
                    maxsplit=1,
                    flags=re.MULTILINE,
                    )
                df2 = parse_table_nuclide(content, index="ZAM")
                df = SubTable(df1, df2)
            elif self.number == 11:
                SubTable = namedtuple(
                    f"Table{self.number}",
                    "by_energy by_nuclide",
                    )
                content, rest = re.split(
                    table_footer,
                    mat_txt,
                    maxsplit=1,
                    flags=re.MULTILINE,
                    )
                df1 = parse_table_energy(content)
                rest = "\n".join(rest.splitlines()[4:])
                content, rest = re.split(
                    table_footer,
                    rest,
                    maxsplit=1,
                    flags=re.MULTILINE,
                    )
                content = content[:40] + "BR NUBAR" + content[40:]
                df2 = parse_table_nuclide(content, index="ZAM", data_row=2).iloc[2:]
                df2.index = df2.index.astype(float)
                df = SubTable(df1, df2)
                # Type A and Type B are treated as different materials
                if mat in out:
                    mat *= 1000  # this is Type B
            out[mat] = df
        return out


def parse_output_for_tables(text, index="ZAM"):
    if index == "ZAM":
        index_col = 1
    else:
        index_col = 0

    items = re.split(summary_header, text, maxsplit=1)
    table_text = items[1]

    tab_items = re.split(table_header, table_text)
    summary = tab_items[0]
    dct = {
        0: summary,
        }
    for x, y, z in zip(tab_items[1::3], tab_items[2::3], tab_items[3::3]):
        table_number = int(x)
        if table_number == 7:
            logging.warning("skip Table 7")
            continue
        table_name = y.strip()
        mat_items = re.split(material_header, z)[1:]
        mat_dct = {}
        for i, j in zip(mat_items[::2], mat_items[1::2]):  # N_MAT + TEXT
            mat_number = int(i)
            content = re.split(table_footer, j, flags=re.MULTILINE)[0]
            lines = content.splitlines()
            string = io.StringIO(lines[0])
            columns = pd.read_csv(
                string,
                sep="\s+",
                header=None,
                index_col=0,
                ).iloc[:, 1:].squeeze()
            if table_number in [8, 11]:
                begin = 2
                string = io.StringIO("\n".join(lines[begin:]))
                df = pd.read_csv(
                    string,
                    sep="\s+",
                    header=None,
                    index_col=0,
                    )
                df.index.name = "Energy"
            else:
                begin = 3
                string = io.StringIO("\n".join(lines[begin:]))
                df = pd.read_csv(
                    string,
                    sep="\s+",
                    header=None,
                    index_col=index_col,
                    ).iloc[:, 1:]
                df.index.name = index
            df.columns = columns
            mat_dct[mat_number] = df
        dct[table_number] = {
            "name": table_name,
            "materials": mat_dct,
            }
    return dct


def parse_table_nuclide(text, index="ZAM", data_row=3, **kwargs):
    index_col = 1 if index.lower() == "zam" else 0
    begin = data_row
    lines = text.splitlines()
    string = io.StringIO(lines[0])
    columns = pd.read_csv(
        string,
        sep="\s+",
        header=None,
        index_col=0,
        ).iloc[:, 1:].squeeze()
    if len(lines[begin:]) == 0:
        df = pd.DataFrame(columns=columns)
    else:
        string = io.StringIO("\n".join(lines[begin:]))
        df = pd.read_csv(
            string,
            sep="\s+",
            header=None,
            index_col=index_col,
            ).iloc[:, 1:]
    df.index.name = index
    df.columns = columns
    return df.T


def parse_table_energy(text, index="Energy", **kwargs):
    index_col = 0
    begin = 2
    lines = text.splitlines()
    string = io.StringIO(lines[0])
    columns = pd.read_csv(
        string,
        sep="\s+",
        header=None,
        index_col=0,
        ).iloc[:, 1:].squeeze()
    if len(lines[begin:]) == 0:
        df = pd.DataFrame(columns=columns)
    else:
        string = io.StringIO("\n".join(lines[begin:]))
        df = pd.read_csv(
            string,
            sep="\s+",
            header=None,
            index_col=index_col,
            )
    df.index.name = index
    df.columns = columns
    return df.T


def parse_table_element(text, index="Element", **kwargs):
    index_col = 1
    begin = 2
    lines = text.splitlines()
    string = io.StringIO(lines[0])
    columns = pd.read_csv(
        string,
        sep="\s+",
        header=None,
        index_col=0,
        ).iloc[:, 1:].squeeze()
    if len(lines[begin:]) == 0:
        df = pd.DataFrame(columns=columns)
    else:
        string = io.StringIO("\n".join(lines[begin:]))
        df = pd.read_csv(
            string,
            sep="\s+",
            header=None,
            index_col=index_col,
            ).iloc[:, 1:]
    df.index.name = index
    df.columns = columns
    return df.T


def parse_table_contributor(text, index="Main contributor", **kwargs):
    index_col = 0
    begin = 1
    lines = text.splitlines()
    string = io.StringIO(lines[0])
    columns = pd.read_csv(
        string,
        sep="\s+",
        header=None,
        index_col=0,
        ).iloc[:, 1:].squeeze()
    if len(lines[begin:]) == 0:
        df = pd.DataFrame(columns=columns)
    else:
        string = io.StringIO("\n".join(lines[begin:]))
        df = pd.read_csv(
            string,
            sep="\s+",
            header=None,
            index_col=index_col,
            )
    df.columns = columns
    df.index.name = df.columns.name
    df.columns.name = index
    return df


def parse_table(text, index="ZAM"):
    """
    Given the text of a ALEPH-2 output table, either expressing
    time-energy dependence or time-nuclide dependence, convert the tabulated
    values into a dataframe.

    Parameters
    ----------
    text : `str`
        string containing the text of a ALEPH-2 output table.
    index : `str`, optional, default is `"ZAM"`
        index of the dataframe. If the output values are tabulated as a
        function of energy, then `index` **must** be `"Energy"`. Else the
        `index` can be either `"ZAM"` or others. In the first case the ZAM
        number will be used as row index. In the second case The Z-SYM-Am
        identifier will be used.

    Returns
    -------
    df : `pandas.DataFrame`
        dataframe with time-energy or time-nuclide dependent outputs.

    """
    begin = 3
    index_col = 0
    if index == "ZAM":
        index_col = 1
    elif index == "Energy":
        begin = 2
    lines = text.splitlines()
    string = io.StringIO(lines[0])
    columns = pd.read_csv(
        string,
        sep="\s+",
        header=None,
        index_col=0,
        ).iloc[:, 1:].squeeze()
    string = io.StringIO("\n".join(lines[begin:]))
    df = pd.read_csv(
        string,
        sep="\s+",
        header=None,
        index_col=index_col,
        ).iloc[:, 1:]
    df.index.name = index
    df.columns = columns
    return df


def _search_pattern_cells(text):
    cells = PATTERN_CELLS.search(text).group("data") \
        .replace("\n", "") \
        .replace(" ", "") \
        .split(",")
    return list(map(int, cells))


def _search_pattern_float(text, pattern):
    data = pattern.search(text).group("data")
    return float(data)


def _search_pattern_flist(text, pattern):
    data = pattern.search(text).group("data").split()
    return list(map(float, data))


def parse_materials_output(text):
    materials = {}
    for k, v in grouper(PATTERN_MATERIAL.split(text)[1:], 2):
        mat_number = int(k)
        cells = _search_pattern_cells(v)
        volume = _search_pattern_float(v, PATTERN_VOLUME)
        time = _search_pattern_flist(v, PATTERN_TIME)
        wdensity = _search_pattern_flist(v, PATTERN_WDENSITY)
        adensity = _search_pattern_flist(v, PATTERN_ADENSITY)
        temperature = _search_pattern_flist(v, PATTERN_TEMPERATURE)
        src = _search_pattern_flist(v, PATTERN_SOURCE)
        nflux = _search_pattern_flist(v, PATTERN_NFLUX)
        power = _search_pattern_flist(v, PATTERN_POWER)
        if PATTERN_FISSIONS.search(v):
            burnup = _search_pattern_flist(v, PATTERN_BURNUP)
            cumburnup = itertools.accumulate(burnup)
            totburnup = _search_pattern_float(v, PATTERN_TOTBURNUP)
            fissions = _search_pattern_flist(v, PATTERN_FISSIONS)
            totfissions = _search_pattern_float(v, PATTERN_TOTFISSIONS)
            material = FissMat(
                ID=mat_number,
                cells=cells,
                volume=volume,
                time=pd.Series(
                    time,
                    name="time [days]",
                    ),
                weight_density=pd.Series(
                    wdensity,
                    index=time,
                    name="weight density [g/cm3]",
                    ),
                atomic_density=pd.Series(
                    adensity,
                    index=time,
                    name="atomic density [g/cm/b]",
                    ),
                temperature=pd.Series(
                    temperature,
                    index=time,
                    name="temperature [K]",
                    ),
                source_strenght=pd.Series(
                    src,
                    index=time,
                    name="source strength [part/s]",
                    ),
                neutron_flux=pd.Series(
                    nflux,
                    index=time,
                    name="neutron flux [n/cm2/s]",
                    ),
                power=pd.Series(
                    power,
                    index=time,
                    name="thermal power [MW]",
                    ),
                burnup=pd.Series(
                    burnup,
                    index=time,
                    name="burnup [GWd/tHM]",
                    ),
                cumulative_burnup=pd.Series(
                    cumburnup,
                    index=time,
                    name="cumulative burnup [GWd/tHM]",
                    ),
                total_burnup=totburnup,
                fissions=pd.Series(
                    fissions,
                    index=time,
                    name="fissions [#]",
                    ),
                total_fissions=totfissions,
            )
        else:
            material = Mat(
                ID=mat_number,
                cells=cells,
                volume=volume,
                time=pd.Series(
                    time,
                    name="time [days]",
                    ),
                weight_density=pd.Series(
                    wdensity,
                    index=time,
                    name="weight density [g/cm3]",
                    ),
                atomic_density=pd.Series(
                    adensity,
                    index=time,
                    name="atomic density [g/cm/b]",
                    ),
                temperature=pd.Series(
                    temperature,
                    index=time,
                    name="temperature [K]",
                    ),
                source_streght=pd.Series(
                    src,
                    index=time,
                    name="source strength [part/s]",
                    ),
                neutron_flux=pd.Series(
                    nflux,
                    index=time,
                    name="neutron flux [n/cm2/s]",
                    ),
                power=pd.Series(
                    power,
                    index=time,
                    name="thermal power [MW]",
                    ),
                )
        materials[mat_number] = material
    return materials


def read_output(file):
    return OutputFile.from_file(file)
