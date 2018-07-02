# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 09:19:24 2018

@author: fiorito_l
"""

import pandas as pd
import numpy as np
import logging, pdb


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

    def perturb(self, pert, **kwargs):
        frame = self.copy()
#        indexName = Xs.index.name
        # Add extra energy points
#        if "energy_point" in kwargs:
#            Xs = Xs.reindex(Xs.index.union(kwargs["energy_point"])).interpolate(method="slinear").fillna(0)
#        Xs.index.name = indexName
        for mat, mt in frame:
            if mat not in pert.index.get_level_values("MAT").unique(): continue
            lmtp = pert.loc[mat].index.get_level_values("MT").unique()
            mtPert = None
            if mt in lmtp:
                mtPert = mt
            else:
                for parent, daughters in sorted(self.__class__.redundant_xs.items(), reverse=True):
                    if mt in daughters and not list(filter(lambda x: x in lmtp, daughters)) and parent in lmtp:
                        mtPert = parent
                        break
            if not mtPert: continue
            P = pert.loc[mat,mtPert]
            P = P.reindex(P.index.union(frame[mat,mt].index)).ffill().fillna(1).reindex(frame[mat,mt].index)
            frame[mat,mt] = frame[mat,mt].multiply(P, axis="index")
            # Negative values are set to zero
            frame[mat,mt][frame[mat,mt] <= 0] = 0
        return Xs(frame).reconstruct_sums()

    def macs(self, E0=0.0253, Elo=1E-5, Ehi=1E1):
        """
        Calculate Maxwellian averaged cross sections.
        """
        from math import sqrt, pi
        from ..integrals.macs import maxw_int, maxw_xs_int
        # add points to the index
        index = set(self.index.values)
        index.update([Elo, Ehi])
        index = np.array(sorted(index))
        index = index[(index >= Elo) & (index <= Ehi)]
        xs = self.reindex(index).interpolate(method='slinear', axis=0).fillna(0)
        data = [[E0,
                 xs.index[i],
                 xs.index[i+1],
                 maxw_int(E0, xs.index[i], xs.index[i+1])
                 ] for i in range(len(xs)-1)]
        dframe = pd.DataFrame(data, columns=["E0", "E1", "E2", "INT"])
        cond = dframe.E1/E0 >= 1e-5
        records = []
        for (mat,mt),x in xs.items():
            data = [[E0,
                     x.index[i],
                     x.iloc[i],
                     x.index[i+1],
                     x.iloc[i+1],
                     maxw_xs_int(E0, x.index[i], x.iloc[i], x.index[i+1], x.iloc[i+1])
                     ] for i in range(len(x)-1)]
            nframe = pd.DataFrame(data, columns=["E0", "E1", "S1", "E2", "S2", "INT"])
            N = nframe[cond].INT.sum(); D = dframe[cond].INT.sum()
            I = N / D * (2/sqrt(pi))
            skipped = "{}/{}".format(sum(cond==False), len(dframe))
            records.append([mat, mt, I, D, Elo, Ehi, E0, skipped])
        return pd.DataFrame(records, columns=["MAT", "MT", "MACS", "FLUX", "Elo", "Ehi", "E0","SKIPPED"])


class Edistr(pd.DataFrame):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.names = ["MAT", "MT", "K", "EIN"]
        self.columns.name = "EOUT"

    def add_points(self, extra_points):
        """
        Add additional incoming energy points.
        """
        frame = self.copy()
        List = []
        for (mat,mt,k),df in frame.groupby(["MAT","MT","K"]):
            grid = sorted((set(df.loc[mat, mt, k].index) | set(extra_points)))
            df = df.reset_index().set_index("EIN").reindex(grid).interpolate(method='slinear').fillna(0).reset_index()
            df["MAT"] = np.round(df.MAT.values).astype(int)
            df["MT"] = np.round(df.MT.values).astype(int)
            df["K"] = np.round(df.K.values).astype(int)
            df = df.set_index(["MAT","MT","K","EIN"])
            List.append(df)
        return Edistr(pd.concat(List, axis=0))

    def normalize(self):
        """
        Normalize each outgoing energy distribution to 1.
        """
        List = []#pd.DataFrame([v/v.sum() for i,v in self.iterrows()])
        for i,v in self.iterrows():
            dx = v.index.values[1:] - v.index.values[:-1]
            y = (v.values[1:]+v.values[:-1])/2
            List.append(v/y.dot(dx))
        frame = pd.DataFrame(List)
        frame.index = pd.MultiIndex.from_tuples(frame.index)
        return Edistr(frame)

    def perturb(self, pert, **kwargs):
        frame = self.copy()
        for (mat,mt,k),S in self.groupby(["MAT", "MT", "K"]):
            if (mat,mt) not in pert.index: continue
            for ein,edistr in S.loc[mat,mt,k].iterrows():
                for (elo,ehi),P in pert.loc[mat,mt].groupby(["ELO","EHI"]):
                    if ein >= elo and ein <= ehi:
                        P = P[elo,ehi]
                        eg = sorted(set(edistr.index) | set(P.index))
                        if len(eg) != len(P):pdb.set_trace()
                        P = P.reindex(eg).ffill().fillna(0).reindex(edistr.index)
                        pedistr = edistr + P
                        frame.loc[mat,mt,k,ein] = pd.Series(np.where(pedistr>0, pedistr, edistr), index=pedistr.index)
        return Edistr(frame).normalize()


class XsCov(pd.DataFrame):
    """
    columns =  (MATi,MTj) ... (MATm,MTn)
    index = E1, E2, ..., El
    """

    def to_matrix(self):
        return self.index, Cov(self.values)

    def get_samples(self, nsmp, **kwargs):
        from ..functions import div0
        index, cov = self.to_matrix()
        frame = pd.DataFrame(cov.sampling(nsmp) + 1, index=index, columns=range(1,nsmp+1))
        frame.columns.name = 'SMP'
        if "eig" in kwargs:
            if kwargs["eig"] > 0:
                eigs = cov.eig()[0]
                idxs = np.abs(eigs).argsort()[::-1]
                dim = min(len(eigs), kwargs["eig"])
                eigs_smp = Cov(np.cov(frame.values)).eig()[0]
                idxs_smp = np.abs(eigs_smp).argsort()[::-1]
                print("MF[31,33] eigenvalues:\n{:^10}{:^10}{:^10}".format("EVAL", "SAMPLES","DIFF %"))
                diff = div0(eigs[idxs]-eigs_smp[idxs_smp], eigs[idxs], value=np.NaN)*100.
                E = ["{:^10.2E}{:^10.2E}{:^10.1F}".format(a,b,c) for a,b,c in zip(eigs[idxs][:dim], eigs_smp[idxs_smp][:dim], diff[:dim])]
                print("\n".join(E))
        return frame

    def macs(self, E0=0.0253, Elo=1E-5, Ehi=1E1):
        from ..integrals.macs import maxw_int
        records = []
        for (mat,mt),sec in self.groupby(["MAT","MT"]):
            C = sec[mat,mt].loc[mat,mt]
            E = set(C.index.values)
            E.update([Elo, Ehi])
            E = np.array(sorted(E))
            E = E[(E >= Elo) & (E <= Ehi)]
            C = C.reindex(E).ffill().fillna(0).T.reindex(E).ffill().fillna(0)
            data = [[E0,
                     E[i],
                     E[i+1],
                     maxw_int(E0, E[i], E[i+1])
                     ] for i in range(len(E)-1)]
            dframe = pd.DataFrame(data, columns=["E0", "E1", "E2", "INT"])
            cond = dframe.E1/E0 >= 1e-5
            D = dframe[cond].INT.sum()
            S = dframe[cond].INT / D
            rvar = S.dot(C.values[:-1,:-1][cond][:,cond].dot(S))
            rstd = np.sqrt(rvar)
            skipped = "{}/{}".format(sum(cond==False), len(dframe))
            records.append([mat, mt, rvar, rstd, D, Elo, Ehi, E0, skipped])
        return pd.DataFrame.from_records(records, columns=["MAT", "MT", "VAR", "STD", "FLUX", "Elo", "Ehi", "E0","SKIPPED"])


class EdistrCov(pd.DataFrame):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.names = ["MAT", "MT", "ELO", "EHI", "EOUT"]
        self.columns.names = ["MAT", "MT", "ELO", "EHI", "EOUT"]

    def to_matrix(self):
        return self.index, Cov(self.values)

    def get_samples(self, nsmp, **kwargs):
        from ..functions import div0
        index, cov = self.to_matrix()
        frame = pd.DataFrame(cov.sampling(nsmp), index=index, columns=range(1,nsmp+1))
        frame.columns.name = 'SMP'
        if "eig" in kwargs:
            if kwargs["eig"] > 0:
                eigs = cov.eig()[0]
                idxs = np.abs(eigs).argsort()[::-1]
                dim = min(len(eigs), kwargs["eig"])
                eigs_smp = Cov(np.cov(frame.values)).eig()[0]
                idxs_smp = np.abs(eigs_smp).argsort()[::-1]
                print("MF35 eigenvalues:\n{:^10}{:^10}{:^10}".format("EVAL", "SAMPLES","DIFF %"))
                diff = div0(eigs[idxs]-eigs_smp[idxs_smp], eigs[idxs], value=np.NaN)*100.
                E = ["{:^10.2E}{:^10.2E}{:^10.1F}".format(a,b,c) for a,b,c in zip(eigs[idxs][:dim], eigs_smp[idxs_smp][:dim], diff[:dim])]
                print("\n".join(E))
        return frame



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


class Cov(np.ndarray):

    def __new__(cls, arr):
        obj = np.ndarray.__new__(cls, arr.shape, arr.dtype)
        obj[:] = arr[:]
        return obj

    @property
    def prefix(self):
        return "COV : "

    @property
    def dim(self):
        """
        Dimension of the covariance.
        """
        length = self.shape[0]
        return length

    @property
    def var(self):
        r"""
        Variance array.
        """
        var = np.diag(np.array(self))
        return var

    @property
    def std(self):
        r"""
        Standard deviation array.
        """
        var = self.var
        if (var < 0).any():
            raise ValueError("Variances must be non-negative")
        std = np.sqrt(var)
        return std

    @property
    def nnegvar(self):
        r"""
        Number of negative variances.
        """
        return np.flatnonzero(self.var < 0).size

    @property
    def nzerovar(self):
        r"""
        Number of zero variances.
        """
        return np.flatnonzero(self.var == 0).size

    def empty_off_diagonals(self):
        r"""
        Remove off-diagonal elements.

        Outputs:
            - :``C``: :
                (``cov.Cov instance``) covariance with empty off-diagonals
        """
        logging.info(self.prefix + "'no_correlations' option is requested, delete off-diagonal terms")
        C = Cov(np.diag(np.diag(self)))
        return C

    def is_symmetric(self):
        r"""
        Check if covariance is symmetric.

        If it is nearly symmetric (rtol=1e-5), then we copy the upper
        triangular part to the lower triangular part and we make it
        symmetric.

        Outputs:
            - :``check``: :
                (boolean) ``True`` if matrix is symmetric, else ``False``
        """
        check = True
        if not (self.T == self).all():
            check = False
            if np.isclose(self.T, self).all():
                check = True
                self[:] = up2down(self)
        return check

    def reduce_size(self):
        """
        Reduce matrix dimensions when zeros are found on the diagonal.

        Outputs:
            * :``nonzero_idxs``: :
                (1d array) positions of the original diagonal matrix where the
                coefficients were not zero
            * :``cov_reduced``: :
                (``cov.Cov`` instance) reduced covariance matrix

        """
        nonzero_idxs =  np.flatnonzero(np.diag(self))
        cov_reduced = self[nonzero_idxs][:,nonzero_idxs]
        return nonzero_idxs, cov_reduced

    def restore_size(self, nonzero_idxs, cov_reduced):
        """
        Restore original matrix dimensions from a reduced matrix and an array
        of positions to convert from reduced to original size.

        Inputs:
            * :``nonzero_idxs``: :
                (1d array) positions of the original diagonal matrix where the
                coefficients were not zero
            * :``cov_reduced``: :
                (``cov.Cov`` instance) reduced covariance matrix

        Outputs:
            * :``cov``: :
                (``cov.Cov`` instance) reduced covariance matrix increased
                to given size according to the indexes given in input
        """
        cov = Cov(np.zeros_like(self))
        for i,ni in enumerate(nonzero_idxs):
            cov[ni,nonzero_idxs] = cov_reduced[i]
        return cov

    def sampling(self, nsmp, pdf='normal'):
        r"""
        Extract random samples from the covariance matrix, either using
        the cholesky or the eigenvalue decomposition.

        Inputs:
            - :``nsmp``: :
                (integer) number of samples

        Outputs:
            - :``samples``: :
                (array) random samples
        """
        logging.debug(self.prefix + "Covariance matrix dimension is {} X {}".format(*self.shape))
        y = np.random.randn(self.dim, int(nsmp))
        nonzero_idxs, cov_reduced = self.reduce_size()
        nzeros = self.shape[0] - len(nonzero_idxs)
        if nzeros > 0:
            logging.debug(self.prefix + "Found {} zeros on the diagonal, reduce matrix dimension to {} X {}".format(nzeros, *cov_reduced.shape))
        try:
            L_reduced = cov_reduced.cholesky()
            logging.debug(self.prefix + "Cholesky decomposition was successful")
        except np.linalg.linalg.LinAlgError as exc:
            logging.debug(self.prefix + "Cholesky decomposition was not successful, proceed with eigenvalue decomposition")
            L_reduced = cov_reduced.eigendecomp()
        L = self.restore_size(nonzero_idxs, L_reduced)
        samples = np.array(L.dot(y), dtype=float)
        return samples

    @property
    def corr(self):
        r"""
        Correlation matrix.
        """
        from sandy.functions import div0
        if not self.is_symmetric():
            raise ValueError("Covariance matrix must be square and symmetric")
        coeff = div0(1, self.std)
        corr = np.multiply(np.multiply(self.T, coeff).T, coeff)
        return corr

    def cholesky(self):
        r"""
        Perform a Cholesky decomposition of the covariance matrix.

        Outputs:
            - :``L``: :
                (2d array) lower triangular matrix
        """
        from scipy.linalg import cholesky
        L = cholesky(self, lower=True, overwrite_a=False, check_finite=False)
        return L

    def eig(self):
        r"""
        Extract eigenvalues and eigenvectors of the covariance matrix.

        Outputs:
            - :``E``: :
                (1d-array) eigenvalues
            - :``V``: :
                (2d-array) eigenvectors
        """
        from scipy.linalg import eig
        E, V = eig(self)
        E, V = E.real, V.real
        return E, V

    def eigendecomp(self):
        r"""
        Perform an eigenvalue decomposition of the covariance matrix.

        Outputs:
            - :``L``: :
                (2d-array) lower triangular matrix
        """
        from scipy.linalg import qr
        E, V = self.eig()
        negative_eig = np.extract(E < 0, E)    # extract negative eigenvalues
        if len(negative_eig) != 0:
            largest_negative = max(abs(negative_eig))
            logging.debug(self.prefix + '{} negative eigenvalues were found and replaced with zero'.format(negative_eig.size))
            pos = sorted(abs(E),reverse=True).index(largest_negative) + 1
            logging.debug(self.prefix + 'Largest negative eigenvalue ranks {}/{}'.format(pos, E.size))
            logging.debug(self.prefix + 'eig(-)/eig_max = {}%'.format(largest_negative/max(abs(E))*100.))
        E[E<=0] = 0
        Esqrt = np.diag(np.sqrt(E))
        M = V.dot(Esqrt)
        Q,R = qr(M.T)
        L = R.T
        logging.debug(self.prefix + "Eigenvalue decomposition was successful")
        return L

    def plot(self):
        r"""
        Plot covariance matrix as a pseudocolor plot of a 2-D array.
        The colorbar is also added to the figure.
        """
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        pcm = ax.matshow(self.corr, vmin=-1, vmax=1, cmap='bwr', aspect='auto')
        # Resize the plot to make space for the colorbar
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, 0.7, box.height])
        # set labels
        ax.set_title('evaluated correlation matrix')
        ax.set_xlabel('energy (eV)')
        ax.set_ylabel('energy (eV)')
        # Plot the colorbar in desired position
        cbaxes = fig.add_axes([0.85, 0.1, 0.03, 0.8])
        plt.colorbar(pcm, cax=cbaxes)
        plt.show()
#        fig.show()

    def dump(self, fname):
        np.savetxt(fname, self, fmt='%.5e')