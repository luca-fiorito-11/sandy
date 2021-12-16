# -*- coding: utf-8 -*-
"""
This module contains all classes and functions specific for processing
covariance matrix with scipy.sparse
"""
import pandas as pd
import numpy as np
import scipy
import scipy.sparse as sps
import sandy

__author__ = "Aitor Bengoechea"


class cov:
    def __init__(self, datos):
        self.data = datos

    @property
    def data(self):
        """
        Introduce a DataFrame and transfor into scipy sparse.

        Returns
        -------
        `dict`
            hierarchical RDD content
        """
        return self._data

    @data.setter
    def data(self, data):
        """[summary]

        Args:
            data (pd.DataFrame): [description]
        """
        data_ = np.array(data)
        self._data = sps.csr_matrix(data_)

    @property
    def size(self):
        return self.data.shape[0]

    def eig(self, sort=True):
        return get_eig(self.data, sort=sort)

    @property
    def std(self):
        cova = self.data
        std = np.sqrt(cova.diagonal())
        return sps.csr_matrix((std))

    @property
    def corr(self):
        cov = self.data
        coeff = self.std
        corr = cov/coeff.multiply(coeff.T)
        return corr

    def invert(self):
        cov = self.data
        cov_inv = sps.linalg.inv(cov).toarray()
        return self.__class__(cov_inv)

    def decompose(self):
        E, V = self.eig(sort=False)
        E[E < 0] = 0
        M = sps.csr_matrix(V).dot(sps.diags(np.sqrt(E))).toarray()
        Q, R = scipy.linalg.qr(M.T)
        return R.T

    @classmethod
    def from_var(cls, var):
        var_ = pd.Series(var).values
        cov = sps.diags(var_)
        return cls(cov.toarray())

    @classmethod
    def from_stdev(cls, std):
        std_ = pd.Series(std).values
        var = std_ * std_
        return cov.from_var(var)

    def _gls_Vy_calc(self, S):
        S_ = pd.DataFrame(S).values
        S_ = sps.csr_matrix(S_)
        Vx_prior = self.data
        Vy_calc = S_.dot(Vx_prior).dot(S_.T)
        return Vy_calc

    def _gls_G(self, S, Vy_extra):
        Vy_extra_ = pd.DataFrame(Vy_extra).values
        Vy_extra_ = cov(Vy_extra_).data
        # GLS_sensitivity:
        Vy_calc = self._gls_Vy_calc(S)
        G = Vy_calc + Vy_extra_
        return G

    def _gls_G_inv(self, S, Vy_extra):
        G = self._gls_G(S, Vy_extra).toarray()
        return cov(G).invert()

    def _gls_general_sensitivity(self, S, Vy_extra, threshold=None):
        S_ = sps.csr_matrix(pd.DataFrame(S).values)
        Vx_prior = self.data
        # GLS_sensitivity:
        G_inv = self._gls_G_inv(S, Vy_extra).data
        sensitivity = Vx_prior.dot(S_.T).dot(G_inv)
        if threshold is not None:
            sensitivity[sensitivity < threshold] = 0
        return sensitivity


def get_eig(cov, sort=True):
    try:
        E, V = sps.linalg.eig(cov)
    except:  # Because sometimes sps gives a error to send you to this method
        E, V = scipy.linalg.eig(cov.toarray())
    E = pd.Series(E.real)
    V = pd.DataFrame(V.real)
    if sort:
        idx = E.sort_values(ascending=False).index
        E = E.iloc[idx].reset_index(drop=True)
        V = V.iloc[idx].reset_index(drop=True)
    return E.values, V.values
