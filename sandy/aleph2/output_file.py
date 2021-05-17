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
            msg = "tables were not parsed, run first 'parse_output'"
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
        for i, j in zip(mat_items[::2], mat_items[1::2]):
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
