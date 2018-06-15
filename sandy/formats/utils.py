# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 09:19:24 2018

@author: fiorito_l
"""

import pandas as pd
import numpy as np



class Section(dict):
    pass



class Xs(pd.DataFrame):

    redundant_xs = {107 : range(800,850),
                    106 : range(750,800),
                    105 : range(700,750),
                    104 : range(650,700),
                    103 : range(600,650),
                    101 : range(102,118),
                    18 : (19,20,21,38),
                    27 : (18,101),
                    4 : range(50,92),
                    3 : (4,5,11,16,17,*range(22,38),41,42,44,45),
                    1 : (2,3),
                    452 : (455,456)}

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.name = "E"
        self.columns.names = ["MAT", "MT"]

    def reconstruct_sums(self, drop=True):
        """
        Reconstruct redundant xs.
        """
        frame = self.copy()
        for mat in frame.columns.get_level_values("MAT").unique():
            for parent, daughters in sorted(Xs.redundant_xs.items(), reverse=True):
                daughters = [ x for x in daughters if x in frame[mat].columns]
                if daughters:
                    frame[mat,parent] = frame[mat][daughters].sum(axis=1)
            # keep only mts present in the original file
            if drop:
                todrop = [ x for x in frame[mat].columns if x not in self.columns.get_level_values("MT") ]
                frame.drop(pd.MultiIndex.from_product([[mat], todrop]), axis=1, inplace=True)
        return Xs(frame)

    def update_tape(self, tapein):
        from copy import deepcopy
        tape = pd.DataFrame(index=tapein.index.copy(), columns=tapein.columns.copy())
        for k,row in tapein.iterrows():
            tape.loc[k].DATA = deepcopy(row.DATA)
            tape.loc[k].TEXT = deepcopy(row.TEXT)
        for mat, mt in self:
            mf = 1 if mt in (452,455,456) else 3
            name = 'NUBAR' if mt in (452,455,456) else 'XS'
            if (mat, mf, mt) not in tape.index:
                continue
            # Cut threshold xs
            iNotZero = next((i for i, x in enumerate(self[mat,mt]) if x), None)
            if iNotZero > 0:
                SeriesXs = self[mat,mt].iloc[iNotZero-1:]
            else:
                SeriesXs = self[mat,mt]
            # Assume all xs have only 1 interpolation region and it is linear
            tape.DATA.loc[mat,mf,mt][name] = SeriesXs
            tape.DATA.loc[mat,mf,mt]["NBT"] = [len(SeriesXs)]
            tape.DATA.loc[mat,mf,mt]["INT"] = [2]
        return Endf6(tape)

    def perturb(self, pert, **kwargs):
        frame = self.copy()
#        indexName = Xs.index.name
        # Add extra energy points
#        if "energy_point" in kwargs:
#            Xs = Xs.reindex(Xs.index.union(kwargs["energy_point"])).interpolate(method="slinear").fillna(0)
#        Xs.index.name = indexName
        for mat, mt in frame:
            if mat not in pert.index.get_level_values("MAT").unique():
                continue
            lmtp = pert.loc[mat].index.get_level_values("MT").unique()
            mtPert = None
            if mt in lmtp:
                mtPert = mt
            else:
                for parent, daughters in sorted(self.__class__.redundant_xs.items(), reverse=True):
                    if mt in daughters and not list(filter(lambda x: x in lmtp, daughters)) and parent in lmtp:
                        mtPert = parent
                        break
            if not mtPert:
                continue
            P = pert.loc[mat,mtPert]
            P = P.reindex(P.index.union(frame[mat,mt].index)).ffill().fillna(1).reindex(frame[mat,mt].index)
            frame[mat,mt] = frame[mat,mt].multiply(P, axis="index")
            # Negative values are set to zero
            frame[mat,mt][frame[mat,mt] <= 0] = 0
        return Xs(frame).reconstruct_sums()



class XsCov(pd.DataFrame):
    """
    columns =  (MATi,MTj) ... (MATm,MTn)
    index = E1, E2, ..., El
    """

    pass


def triu_matrix(arr, size):
    """
    Given the upper triangular values of a **square symmetric** matrix in
    an array, return the full matrix.

    Inputs:
        - arr :
            (1d array) array with the upper triangular values of the matrix
        - size :
            (int) dimension of the matrix

    Outputs:
        - matrix :
            (2d array) reconstructed 2d-array with symmetric matrix
    """
    matrix = np.zeros([size, size])
    indices = np.triu_indices(size)
    matrix[indices] = arr
    matrix += np.triu(matrix, 1).T
    return matrix



def up2down(C):
    """
    Given a covariance matrix in input, copy the upper triangular part to the
    lower triangular part.

    Inputs:
        - C :
            (2d-array) input covariance matrix

    Outputs:
        - C1 :
            (2d-array) output covariance matrix
    """
    U = np.triu(C)
    L = np.triu(C, 1).T
    C1 = U + L
    return C1



def corr2cov(corr, s):
    dim = corr.shape[0]
    S = np.repeat(s, dim).reshape(dim, dim)
    cov = S.T * (corr * S)
    cov = up2down(cov)
    return cov