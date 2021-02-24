import pdb
import logging

import numpy as np
import pandas as pd

import sandy

__author__ = "Luca Fiorito"
__all__ = [
        "Pert",
        ]


class Pert():
    """
    Container for tabulated multigroup perturbation coefficients.

    Attributes
    ----------
    data : `pandas.Series`
        series of groupwise perturbation coefficients
    left : `pandas.Series`
        perturbation coefficients with left-bounds of the energy intervals as
        index
    mid : `pandas.Series`
        perturbation coefficients with mid-values of the energy intervals as
        index
    right : `pandas.Series`
        perturbation coefficients with right-bounds of the energy intervals as
        index

    Methods
    -------
    reshape
        interpolate perturbation coefficients over new energy grid structure
    """
    _indexname = "ENERGY"

    def __repr__(self):
        return self.data.__repr__()

    def __init__(self, series, **kwargs):
        self.data = pd.Series(series, dtype=float, **kwargs)

    @property
    def data(self):
        """
        Series of groupwise perturbation coefficients.

        Attributes
        ----------
        index : `pandas.IntervalIndex`
            energy bins (in eV) over which the perturbations are defined
        values : `numpy.array`
            perturation coeffcients as ratio values

        Returns
        -------
        `pandas.Series`
            multigroup perturbation coefficients

        Raises
        ------
        `sandy.Error`
            if `data` is not a `pandas.Series`
        `sandy.Error`
            if energy grid is not monotonically increasing
        """
        return self._data

    @data.setter
    def data(self, data):
        self._data = data
        if isinstance(self._data.index, pd.IntervalIndex):
            self._data.index = self._data.index.right
        index = self._data.index.astype(float)
        if not index.is_monotonic_increasing:
            raise sandy.Error("energy grid is not monotonically increasing")
        self._data.index = pd.IntervalIndex.from_breaks(index.insert(0, 0))
        self._data.index.name = self._indexname

    @property
    def right(self):
        """
        Pertrubations with right-bounds of the energy intervals as index.
        .. note:: the index of the series is a `pandas.Index` instance, not a
                  `pandas.IntervalIndex` instance.

        Returns
        -------
        `pandas.Series`
            perturbations
        """
        return pd.Series(self.data.values, index=self.data.index.right)

    @property
    def left(self):
        """
        Pertrubations with left-bounds of the energy intervals as index.
        .. note:: the index of the series is a `pandas.Index` instance, not a
                  `pandas.IntervalIndex` instance.

        Returns
        -------
        `pandas.Series`
            perturbations
        """
        return pd.Series(self.data.values, index=self.data.index.left)

    @property
    def mid(self):
        """
        Pertrubations with mid-values of the energy intervals as index.
        .. note:: the index of the series is a `pandas.Index` instance, not a
                  `pandas.IntervalIndex` instance.

        Returns
        -------
        `pandas.Series`
            perturbations
        """
        return pd.Series(self.data.values, index=self.data.index.mid)

    def reshape(self, eg, inplace=False):
        """
        Interpolate perturbation over new energy grid structure using `bfill`
        method.

        .. note:: value `1` is given to extrapolated energies

        Parameters
        ----------
        eg : array-like object
            monotonic energy grid with non-negative values
        inplace : `bool`, optional, default is `False`
            apply changes **inplace** and return `None`

        Returns
        -------
        `sandy.Pert`
            perturbation coefficients reshaped over a union of the original and
            new energy grids

        Raises
        ------
        `aleph.Error`
            if the given energy grid is not monotonic
        `value.Error`
            if negative values are found in the given energy grid
        """
        index = pd.Index(eg)
        if not index.is_monotonic_increasing:
            raise sandy.Error("energy grid is not monotonic increasing")
        if (index < 0).any():
            raise ValueError("found negative values in the energy grid")
        enew = self.right.index.union(index).unique().astype(float).values
        # remove zero if any, it will be automatically added by `Pert`
        enew = enew[enew != 0]
        # this part prevents errors in "scipy.interp1d" when x.size == 1
        x = self.right.index.values
        y = self.right.values
        if y.size == 1:
            # this must be done after that enew is created
            x = np.insert(x, 0, 0)
            y = np.insert(y, 0, 0)
        pertnew = sandy.shared.reshape_bfill(
                          x,
                          y,
                          enew,
                          left_values=1,
                          right_values=1,
                          )
        series = pd.Series(pertnew, index=enew)
        if not inplace:
            return Pert(series)
        self.data = series

    def truncate(self, low=0.0, high=2.0):
        """
        Truncate perturbation values when they exceed a defined lower or
        upper value.
        Every value outside the boundary is replaced with the boundary value
        that was exceeded.

        Parameters
        ----------
        low : `float`, optional
            lower limit for truncation. The default is 0.0.
        high : `float`, optional
            upper limit for truncation. The default is 2.0.

        Returns
        -------
        `sandy.Pert`
            perturbation instance with truncated values.
        """
        data = self.data.copy()
        data = data.where(data <= high, high).where(data >= low, low)
        return self.__class__(data)

    def recenter(self, low=0.0, high=2.0, value=1.):
        """
        Truncate perturbation values when they exceed a defined lower or
        upper value.
        Every value outside the boundary is replaced with the cental
        perturbation vaule, i.e. 1, or with a value of choice.

        Parameters
        ----------
        low : `float`, optional
            lower limit for truncation. The default is 0.0.
        high : `float`, optional
            upper limit for truncation. The default is 2.0.
        value : `float`, optional
            value used to replace any perturbation exceeding the
            limiting domain. The default is 1.0.

        Returns
        -------
        `sandy.Pert`
            perturbation instance with truncated values.
        """
        data = self.data.copy()
        data = data.where(data <= high, value).where(data >= low, value)
        return self.__class__(data)

    @classmethod
    def from_file(cls, file, sep=None, **kwargs):
        """
        Initialize `Pert` object reading perturbations from file.

        Parameters
        ----------
        file : `str`
            file name (absolute or relative path)
        sep : `str`, optional, default `"\s+"`
            column separator. By default it takes blankspaces as separators.
            .. note:: for `csv` files use `","`
        **kwargs : `pandas.read_csv` properties, optional

        Returns
        -------
        `Pert`
            Container for binned perturbations
        """
        data = np.genfromtxt(file, dtype=float, delimiter=sep, **kwargs)
        if data.ndim < 2:
            raise sandy.Error("at least 2 columns should be given in the file")
        series = pd.Series(data[:, 1], index=data[:, 0])
        return Pert(series)
