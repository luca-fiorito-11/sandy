"""
This module contains all classes and functions to provide a API for a MCNP
*meshtal* file.
"""

import re
import io
import logging

import pandas as pd
import numpy as np


__author__ = "Luca Fiorito"

__all__ = [
        "MshtTally",
        ]


class MshtTally():
    """
    Container for a Meshtal tally.

    Attributes
    ----------
    data : `pandas.DataFrame`
        tabulated mesthal tally with column names according to mesthal output
        file.

    Methods
    -------
    filter_by
        filter meshtally data according to given query
    from_file
        parse a meshtal output file and return tallies
    get_mesh
        pivot meshtally data according to two axes (e.g. `'X'` and `'Y'`)
        and return the coordinates and values on a meshgrid
    pivot_table
        pivot meshtally data according to two axes (e.g. `'X'` and `'Y'`)
    plot
        plot meshtally data using `matplotlib`
    """

    def __init__(self, data, *args, **kwargs):
        self.data = data

    def pivot_table(self, x="X", y="Y", values="Result", aggfunc="sum"):
        """
        Pivot meshtally data according to two axes (e.g. `'X'` and `'Y'`)

        Parameters
        ----------
        x : `str`, optional, default is `'X'`
            data on the plot x-axis.
            Other options are `'X'`, '`Y`', `'Z'`, `'Energy'`.
        y : `str`, optional, default is `'Y'`
            data on the plot y-axis.
            Other options are `'X'`, '`Y`', `'Z'`, `'Energy'`.
        values : `str`, optional, default is `Result`
            data to plot.
            Other options are:
                - `'RelError'` : relative error (fraction)
                - `'RsltVol'` : `'Result'` mutliplied by `''Volume'`, only
                                if *out=cf* was used in *fmesh*
        aggfunc : function, optional, default is `'sum'`
            function to aggregate values (as in `pandas.pivot_table`).
            Other options are `'mean'`, `'min'`, `'max'`, '`std`'.

        Returns
        -------
        `pandas.DataFrame`
            An Excel-style pivot table
        """
        if x == y:
            raise ValueError(f"x='{x}' and y='{y}' must be different")
        return self.data.pivot_table(
                values=values,
                index=x,
                columns=y,
                aggfunc=aggfunc,
                )

    def filter_by(self, query):
        """
        Filter meshtally data according to given query.

        Parameters
        ----------
        query : `str`
            The query string used to filter `data` columns.
            Examples are `'X > -10 & X < 10'` or `'Z == 0'`.

        Returns
        -------
        `MshTally`
            filtered meshtally object
        """
        out = self.data.query(query)
        if out.empty:
            raise ValueError("no cell matches applied query")
        return self.__class__(out)

    def get_mesh(self, x="X", y="Y", values="Result", aggfunc="sum"):
        """
        Pivot meshtally data according to two axes (e.g. `'X'` and `'Y'`)
        and return the coordinates and values on a meshgrid.

        Parameters
        ----------
        x : `str`, optional, default is `'X'`
            data on the plot x-axis.
            Other options are `'X'`, '`Y`', `'Z'`, `'Energy'`.
        y : `str`, optional, default is `'Y'`
            data on the plot y-axis.
            Other options are `'X'`, '`Y`', `'Z'`, `'Energy'`.
        values : `str`, optional, default is `Result`
            data to plot.
            Other options are:
                - `'RelError'` : relative error (fraction)
                - `'RsltVol'` : `'Result'` mutliplied by `''Volume'`, only
                                if *out=cf* was used in *fmesh*
        aggfunc : function, optional, default is `'sum'`
            function to aggregate values (as in `pandas.pivot_table`).
            Other options are `'mean'`, `'min'`, `'max'`, '`std`'.

        Returns
        -------
        Three *NxM* arrays
            2D arrays for `x`, `y` and `values` that can be used for a
            meshplot
        """
        tab = self.pivot_table(x=x, y=y, values=values, aggfunc=aggfunc)
        X, Y = np.meshgrid(tab.index.values, tab.columns.values)
        return X.T, Y.T, tab.values

    def plot(self, ax, x="X", y="Y", values="Result",
             aggfunc="sum", query=None,
             **kwargs):
        """
        Plot meshtally data using `matplotlib`.

        Parameters
        ----------
        ax : `matplotlib.axes.Axes`
            `Axes` instance used for plotting
        x : `str`, optional, default is `'X'`
            data on the plot x-axis.
            Other options are `'X'`, '`Y`', `'Z'`, `'Energy'`.
        y : `str`, optional, default is `'Y'`
            data on the plot y-axis.
            Other options are `'X'`, '`Y`', `'Z'`, `'Energy'`.
        values : `str`, optional, default is `Result`
            data to plot.
            Other options are:
                - `'RelError'` : relative error (fraction)
                - `'RsltVol'` : `'Result'` mutliplied by `''Volume'`, only
                                if *out=cf* was used in *fmesh*
        aggfunc : function, optional, default is `'sum'`
            function to aggregate values (as in `pandas.pivot_table`).
            Other options are `'mean'`, `'min'`, `'max'`, '`std`'.
        query : `str`, optional, default is `None`
            The query string used to filter `data` columns.
            Examples are `'X > -10 & X < 10'` or `'Z == 0'`.
        **kwargs : keyword arguments
            Additionally, the following arguments are allowed.
            They are passed along to the `matplotlib.axes.Axes.pcolormesh`
            constructor.

        Returns
        -------
        `matplotlib.collections.QuadMesh`
            `matplotlib.axes.Axes.pcolormesh` output

        Notes
        -----
        Input parameter `ax` is affected by this method.
        """
        obj = self.filter_by(query) if query is not None else self
        X, Y, Z = obj.get_mesh(x=x, y=y, values=values, aggfunc=aggfunc)
        pcm = ax.pcolormesh(
                    X,
                    Y,
                    Z,
                    **kwargs
                    )
        ax.set_xlabel(f"{x}")
        ax.set_ylabel(f"{y}")
        return pcm

    @classmethod
    def from_file(cls, file="meshtal", encoding="utf8", errors='ignore'):
        """
        Parse a meshtal output file and return tallies.

        Parameters
        ----------
        filename : `str`, optional, default is `'meshtal'`
            name of the meshtal file
        encoding : `str`, optional, default is `'utf8'`
            option passed to python built-in function `open`
        errors : `str`, optional, default is `'ignore'`
            option passed to python built-in function `open`

        Returns
        -------
        `dict`
            dictionary of (tally numbers, `MshtTally` objects)
        """
        with open(file, encoding=encoding, errors=errors) as f:
            string = f.read()
        blocks = string.split("Mesh Tally Number")
        out = {}
        for block in blocks[1:]:
            sections = block.split("\n\n")
            idx = int(sections[0].splitlines()[0])
            data = re.sub(
                    "Rslt \* Vol",
                    "RsltVol",
                    re.sub("Rel Error", "RelError", sections[2]),
                    )
            df = pd.read_csv(io.StringIO(data), sep="\s+")
            out[idx] = cls(df)
        return out


def from_file(*args, **kwargs):
    """
    Same as `MshtTally.from_file`.
    """
    logging.warning("DEPRECATED! Use 'MshtTally.from_file'")
    return MshtTally.from_file(*args, **kwargs)
