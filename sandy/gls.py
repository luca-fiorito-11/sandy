# -*- coding: utf-8 -*-
"""
This module contains all classes and functions specific for processing GLS.
"""
import numpy as np
import pandas as pd
import sandy

__author__ = "Aitor Bengoechea"
__all__ = [
        "GLS",
        ]


class GLS():
    """
    Container for processing GLS.

    Attributes
    ----------
    data : `pandas.Series`
        Vector where GLS is going to be applied.

    Methods
    -------
    gls_update:
        Perform GlS update.
    """

    def __repr__(self):
        return self.data.__repr__()

    def __init__(self, *args, **kwargs):
        self.data = pd.Series(*args, **kwargs)

    @property
    def data(self):
        """
        Vector where GLS is going to be applied as pandas.Series

        Attributes
        ----------
        index : `pandas.Index` or `pandas.MultiIndex`
            indices

        Returns
        -------
        `pandas.Series`
            Vector where GLS is going to be applied.

        """
        return self._data

    @data.setter
    def data(self, data):
        self._data = pd.Series(data, dtype=float)

    def _gls_sensitivity(self, S, x_p, threshold=None):
        """
        Method to obtain sensitivity according to GLS

        Parameters
        ----------
        S : 2D iterable
            Sensitivity square matrix (MXM).
        x_p : 1D iterable
            Extra vector (MX1).
        threshold : `int`, optional
            threshold to avoid numerical fluctuations. The default is None.

        Returns
        -------
        `GLS`
            GLS sensitivity for a given Vy and S.

        Example
        -------
        >>> S = np.array([[1, 2], [3, 4]])
        >>> sen = sandy.GLS([1, 1])
        >>> sen._gls_sensitivity(S, np.array([1, 1]))
                      0	          1
        0	-2.00000e-01	2.28571e-01
        1	2.00000e-01	5.71429e-02

        >>> S = pd.DataFrame([[1, 2], [3, 4]], columns =[1, 2],index=[3, 4])
        >>> sen = sandy.GLS([1, 1])
        >>> sen._gls_sensitivity(S,pd.Series([1, 1], index=[1, 2]))
                      1	          2
        1	-2.00000e-01	2.28571e-01
        2	2.00000e-01	5.71429e-02
        """
        columns = pd.DataFrame(S).columns
        S_ = pd.DataFrame(S).values
        x_p_ = pd.DataFrame(np.diag(pd.Series(x_p))).values
        # GLS_sensitivity:
        x = np.diag(self.data)
        M = S_.T.dot(x).dot(S_) + x_p_
        M_inv = sandy.CategoryCov(M).invert()
        sensitivity = x.dot(S_).dot(M_inv.data.values)
        if threshold is not None:
            sensitivity[sensitivity < threshold] = 0
        return pd.DataFrame(sensitivity, index=columns, columns=columns)

    def gls_update(self, S, Vx, Vy, x_p, x, threshold=None):
        """
        Perform GlS update for a given variance and sensitivity.

        Parameters
        ----------
        Vx : 1D iterable
            Covariance vector (MX1).
        Vy : 1D iterable
            Extra Covariance vector (MX1).
        S : 2D iterable
            Sensitivity square matrix (MXM).
        x : 1D iterable
            x vector.
        x_p : 1D iterable
            x vector perturbed.
        threshold : `int`, optional
            Thereshold to avoid numerical fluctuations. The default is None.

        Returns
        -------
        `GLS`
            GLS method apply to a GLS object for a given Vy, x, x_p and S.

        Example
        -------
        >>> S = np.array([[1, 2], [3, 4]])
        >>> y = GLS([1, 1])
        >>> Vx = pd.Series([1, 1])
        >>> Vy = pd.Series([1, 1], index=[1, 2])
        >>> x = np.array([1, 1])
        >>> x_p = np.array([2, 2])
        >>> y.gls_update(S, Vx, Vy, x_p, x)
        0   1.02857e+00
        1   1.25714e+00
        dtype: float64
        """
        index = self.data.index
        x_ = pd.Series(x).values
        x_p_ = pd.Series(x_p).values
        delta = x_p_ - x_
        Vx_ = sandy.GLS(Vx)
        y = self.data.values
        A = Vx_._gls_sensitivity(S, Vy, threshold).values
        y_new = y + A.dot(delta)
        return self.__class__(y_new, index=index)
