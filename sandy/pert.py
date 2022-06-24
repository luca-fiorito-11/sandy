import pdb
import logging

import numpy as np
import pandas as pd

import sandy

__author__ = "Luca Fiorito"
__all__ = [
        "Pert",
        ]


col = pd.MultiIndex.from_arrays([[2631, 2631], [5, 2]], names=('MAT', 'MT'))
pert = pd.DataFrame([[1, 1.05], [1.05, 1]],
                    index=pd.IntervalIndex.from_breaks(pd.Index([10, 100])
                                                       .insert(0, 0)),
                    columns=col)

pert_ind = pd.Series([1, 1.05], index=pd.IntervalIndex
                                        .from_breaks(pd.Index([10, 100])
                                                       .insert(0, 0)))


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
    from_bin
        generate a `Pert` object from a pair of energy bins
    reshape
        interpolate perturbation coefficients over new energy grid structure
    """
    _indexname = "ENERGY"

    def __repr__(self):
        return self.data.__repr__()

    def __init__(self, df, **kwargs):
        self.data = pd.DataFrame(df, dtype=float, **kwargs)

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
        `pandas.DataFrame`
            multigroup perturbation coefficients

        Raises
        ------
        `sandy.Error`
            if `data` is not a `pandas.Series`
        `sandy.Error`
            if energy grid is not monotonically increasing

        Examples
        --------
        >>> sandy.Pert(pert_ind).right
        	                  0
        1.00000e+01	1.00000e+00
        1.00000e+02	1.05000e+00
        >>> sandy.Pert(pert).right
        MAT	        2631
        MT	        5	        2
        1.00000e+01	1.00000e+00	1.05000e+00
        1.00000e+02	1.05000e+00	1.00000e+00
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

        Examples
        --------
        >>> sandy.Pert(pert_ind).right
        	                  0
        1.00000e+01	1.00000e+00
        1.00000e+02	1.05000e+00

        >>> sandy.Pert(pert).right
        MAT	        2631
        MT	        5	        2
        1.00000e+01	1.00000e+00	1.05000e+00
        1.00000e+02	1.05000e+00	1.00000e+00
        """
        return pd.DataFrame(self.data.values,
                            index=self.data.index.right,
                            columns=self.data.columns)

    @property
    def left(self):
        """
        Pertrubations with left-bounds of the energy intervals as index.
        .. note:: the index of the series is a `pandas.Index` instance, not a
                  `pandas.IntervalIndex` instance.

        Returns
        -------
        `pandas.DataFrame`
            perturbations

        Examples
        --------
        >>> sandy.Pert(pert_ind).left
                              0
        0.00000e+00	1.00000e+00
        1.00000e+01	1.05000e+00

        >>> sandy.Pert(pert).left
        MAT	        2631
        MT	        5	        2
        0.00000e+00	1.00000e+00	1.05000e+00
        1.00000e+01	1.05000e+00	1.00000e+00
        """
        return pd.DataFrame(self.data.values,
                            index=self.data.index.left,
                            columns=self.data.columns)

    @property
    def mid(self):
        """
        Pertrubations with mid-values of the energy intervals as index.
        .. note:: the index of the series is a `pandas.Index` instance, not a
                  `pandas.IntervalIndex` instance.

        Returns
        -------
        `pandas.DataFrame`
            perturbations

        Examples
        --------
        >>> sandy.Pert(pert_ind).mid
                              0
        5.00000e+00	1.00000e+00
        5.50000e+01	1.05000e+00

        >>> sandy.Pert(pert).mid
        MAT	        2631
        MT	        5	        2
        5.00000e+00	1.00000e+00	1.05000e+00
        5.50000e+01	1.05000e+00	1.00000e+00
        """
        return pd.DataFrame(self.data.values,
                            index=self.data.index.mid,
                            columns=self.data.columns)

    def reorder(self, columns):
        """
        Sort the levels in the columns of `Sandy.Pert`

        Parameters
        ----------
        columns : `pd.Multiindex`
            New columns order.

        Returns
        -------
        `sandy.Pert`
            Pert order with sort levels.

        Examples
        --------
        >>> col = pd.MultiIndex.from_arrays([[5, 2], [2631, 2631]], names=('MT', 'MAT'))
        >>> pert = pd.DataFrame([[1, 1.05], [1.05, 1]], index =pd.IntervalIndex.from_breaks(pd.Index([1.94000e+08, 1.96000e+08+1]).insert(0, 0)), columns=col)
        >>> col_new = pd.MultiIndex.from_arrays([[2631, 2631], [2, 5]], names=('MAT', 'MT'))
        >>> sandy.Pert(pert).reorder(col_new)
        MAT                               2631
        MT                                   2           5
        ENERGY
        (0.0, 194000000.0]         1.05000e+00 1.00000e+00
        (194000000.0, 196000001.0] 1.00000e+00 1.05000e+00
        """
        data = self.data
        for name in columns.names:
            if name not in self.data.columns.names:
                print("Levels do not match")
                return
        data_reorder = data.reorder_levels(columns.names, axis=1)
        col = columns.intersection(data_reorder.columns)
        return sandy.Pert(data_reorder[col])

    def reshape(self, eg, inplace=False, right_values=1):
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

        Examples
        --------
        >>> sandy.Pert(pert).reshape([0, 1, 5,  25, 50, 100, 1.95000e+08]).data
        MAT	                    2631
        MT	                    5	        2
                      ENERGY
                  (0.0, 1.0]	1.00000e+00	1.05000e+00
                  (1.0, 5.0]	1.00000e+00	1.05000e+00
                 (5.0, 10.0]	1.00000e+00	1.05000e+00
                (10.0, 25.0]	1.05000e+00	1.00000e+00
                (25.0, 50.0]	1.05000e+00	1.00000e+00
               (50.0, 100.0]	1.05000e+00	1.00000e+00
        (100.0, 195000000.0]	1.00000e+00	1.00000e+00

        >>> sandy.Pert(pert_ind).reshape([0, 1, 5,  25, 50, 100, 1.95000e+08]).data
                                       0
        ENERGY
        (0.0, 1.0]           1.00000e+00
        (1.0, 5.0]           1.00000e+00
        (5.0, 10.0]          1.00000e+00
        (10.0, 25.0]         1.05000e+00
        (25.0, 50.0]         1.05000e+00
        (50.0, 100.0]        1.05000e+00
        (100.0, 195000000.0] 1.00000e+00
        """
        eg_ = pd.Index(eg)
        df = self.right
        if not eg_.is_monotonic_increasing:
            raise sandy.Error("energy grid is not monotonic increasing")
        if (eg_ < 0).any():
            raise ValueError("found negative values in the energy grid")
        enew = df.index.union(eg_).unique().astype(float).values
        # remove zero if any, it will be automatically added by `Pert`
        enew = enew[enew != 0]
        # this part prevents errors in "scipy.interp1d" when x.size == 1
        name = df.index.name
        if df.shape[1] == 1:
            index = df.index.values
            values = df.values
            if values.size == 1:
                # this must be done after that enew is created
                index = np.insert(index, 0, 0)
                values = np.insert(values, 0, 0)
            pertnew = sandy.shared.reshape_bfill(
                          index,
                          values,
                          enew,
                          left_values="first",
                          right_values=right_values,
                          )
            data = pd.DataFrame(pertnew, index=enew, columns=df.columns)
        else:
            data = df.apply(lambda x: sandy.shared.reshape_bfill(
                                    x.index.values,
                                    x.values,
                                    enew,
                                    left_values="first",
                                    right_values=right_values,
                                    )).set_index(enew).rename_axis(name)
        data = data.loc[(data.index >= eg_.min()) &
                        (data.index <= eg_.max()), :]
        if not inplace:
            return self.__class__(data)
        self.data = data

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

        Examples
        --------
        >>> pert = pd.DataFrame([[-1, 2.55], [2.55, -1]], index =pd.IntervalIndex.from_breaks(pd.Index([10, 100]).insert(0, 0)), columns=col)
        >>> sandy.Pert(pert).truncate().data
        MAT	            2631
        MT	            5	        2
               ENERGY
          (0.0, 10.0]	0.00000e+00	2.00000e+00
        (10.0, 100.0]	2.00000e+00	0.00000e+00

        >>> pert_ind = pd.Series([-1, 2.05], index=pd.IntervalIndex.from_breaks(pd.Index([10, 100]).insert(0, 0)))
        >>> sandy.Pert(pert_ind).truncate().data
        	                      0
               ENERGY
          (0.0, 10.0]	0.00000e+00
        (10.0, 100.0]	2.00000e+00
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

        Examples
        --------
        >>> pert = pd.DataFrame([[-1, 2.55], [2.55, -1]], index =pd.IntervalIndex.from_breaks(pd.Index([10, 100]).insert(0, 0)), columns=col)
        >>> sandy.Pert(pert).recenter().data
        MAT	            2631
        MT	            5	        2
               ENERGY
          (0.0, 10.0]	1.00000e+00	1.00000e+00
        (10.0, 100.0]	1.00000e+00	1.00000e+00

        >>> pert_ind = pd.Series([-1, 2.05], index=pd.IntervalIndex.from_breaks(pd.Index([10, 100]).insert(0, 0)))
        >>> sandy.Pert(pert_ind).recenter().data
                                  0
               ENERGY
          (0.0, 10.0]	1.00000e+00
        (10.0, 100.0]	1.00000e+00
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

    @classmethod
    def from_bin(cls, elow, ehigh, coeff):
        """
        Generate a `Pert` object from a pair of energy bins `(elow, ehigh]`.

        Parameters
        ----------
        elow : `float`
            lower energy boundary in eV.
        ehigh : `float`
            upper energy boundary in eV.
        coeff : `float`
            perturbation coefficient.

        Returns
        -------
        pert : `sandy.Pert`
            perturbation object.

        Examples
        --------
        >>> Pert.from_bin(1e-5, 1e-4, 0.05)
                                     0
                 ENERGY
           (0.0, 1e-05]    1.00000e+00
        (1e-05, 0.0001]    1.05000e+00
        """
        pert = cls([1, 1 + coeff], index=[elow, ehigh])
        return pert
