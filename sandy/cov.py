# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 13:35:07 2017

@author: lfiorito
"""
import numpy as np
import logging

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
    mytype = type(C)
    U = np.triu(C)
    L = np.triu(C, 1).T
    C1 = mytype(U + L)
    return C1



def corr2cov(corr, s):
    mytype = type(corr)
    dim = corr.shape[0]
    S = np.repeat(s, dim).reshape(dim, dim)
    cov = mytype(S.T * (corr * S))
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
        from sandy.sandy_input import options
        logging.debug(self.prefix + "Covariance matrix dimension is {} X {}".format(*self.shape))
        y = np.random.randn(self.dim, int(nsmp))
        nonzero_idxs, cov_reduced = self.reduce_size()
        nzeros = self.shape[0] - len(nonzero_idxs)
        if nzeros > 0:
            logging.debug(self.prefix + "Found {} zeros on the diagonal, reduce matrix dimension to {} X {}".format(nzeros, *cov_reduced.shape))
        if "no_correlations" in options:
            if options['no_correlations']:
                cov_reduced = cov_reduced.empty_off_diagonals()
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

class Cov33(Cov):
    """
    The same as the inherited object *Cov*, but with some extra methods.

    The class methods permit to initialize a *Cov33* instance directly when
    reading the ENDF-6 file or by merging already existing instances.

    Method *reshape* exploits the zero-interpolation law inherent in
    the covariance matrix.
    """

    @classmethod
    def merge(cls, *covs):
        """
        Create a new *Cov33* instance by merging multiple *Cov33* instances.

        The new instance has got attributes:
             - _xx :
                 union grid of the *_xx* attributes of the other instances
             - _yy :
                 union grid of the *_yy* attributes of the other instances
             - _cov :
                 sum of the covariance matrices reshaped according to *_xx*
                 and *_yy*

        Inputs:
            - covs :
                positional arguments containing *Cov33* instances
        """
        from sandy.functions import union_grid
        if len(covs) == 0:
            raise ValueError("At least one or more input arguments must be provided")
        for inst in covs:
            if not isinstance(inst, Cov33):
                raise TypeError("All arguments must be of 'Cov33' type")
        uxx = union_grid(*[inst.xx for inst in covs])
        covs = [ cov.reinit(uxx, uxx) for cov in covs ]
        ucov = np.sum([cov for cov in covs], 0)
        return cls(uxx, uxx, ucov)

    @classmethod
    def read_endf(cls, endf):
        """
        Extract data from ``ENDF-6`` file.
        """
        lb = endf.read_cont([0,1,2,4,5])
        mt = endf.read_control()[2]
        logging.debug(" - lb={}".format(lb))
        endf.back()
        if lb == 0 or lb == 1:
            return cls.lb1(endf)
        elif lb == 2:
            return cls.lb2(endf)
        elif lb == 4:
            return cls.lb4(endf)
        elif lb == 5 or lb == 7:
            return cls.lb5(endf)
        elif lb == 6:
            return cls.lb6(endf)
        elif lb == 8:
            return cls.lb8(endf)
        else:
            logging.error("MF33 MT{} : Covariances with 'lb={}' cannot yet be read by SANDY".format(mt, lb))
            sys.exit()

    @classmethod
    def lb1(cls, endf):
        r"""
        Classmethod to initialize the covariance from the text read in `LB=1`
        format of `ENDF-6`.

        ..ENDF-6 description::
            Fractional components correlated only within each :math:`E_k`
            interval.

        ..math::
            rcov(x,y) = f(E_x) \delta(E_x,E_y)

        Inputs:
            - endf :
                *endf.ENDF* instance

        Outputs:
            - inst :
                *Cov33* instance

        Section parameters:
            - ne :
                number of entries in the array {:math:`E_k`} defining
                (NE-1) energy intervals.
        """
        from sandy.cov import Cov
        _, _, lt, lb, nt, ne, data = endf.read_list()
        xx = yy = data[::2]
        f = data[1::2]
        cov = Cov(np.diag(f))
        inst = cls(xx, yy, cov, lb)
        return inst

    @classmethod
    def lb2(cls, endf):
        r"""
        Classmethod to initialize the covariance from the text read in `LB=2`
        format of `ENDF-6`.

        ..ENDF-6 description::
            Fractional components correlated over all :math:E k intervals.

        ..math::
            rcov(x,y) = f(E_x) f(E_y)

        Inputs:
            - endf :
                *endf.ENDF* instance

        Outputs:
            - inst :
                *Cov33* instance

        Section parameters:
            - ne :
                number of entries in the array {:math:`E_k`} defining
                (NE-1) energy intervals.
        """
        from sandy.cov import Cov
        _, _, lt, lb, nt, ne, data = endf.read_list()
        xx = yy = data[::2]
        f = np.array(data[1::2])
        cov = Cov(f*f.reshape(-1,1))
        return cls(xx, yy, cov, lb)

    @classmethod
    def lb4(cls, endf):
        r"""
        Classmethod to initialize the covariance from the text read in `LB=4`
        format of `ENDF-6`.

        Inputs

        ENDF-6 description
        ==================
        Fractional components correlated over all :math:`E_l` intervals
        within each :math:`E_k` interval.
        :math:`rcov(x,y) = f(Ex)*\delta(E_x,E_y)*f'(E_x)*f'(E_y):
        """
        from sandy.functions import union_grid
        _, _, lt, lb, nt, ne, data = endf.read_list()
        xx = data[::2][:-lt]  # first energy grid is until LT points to the end
        yy = data[::2][-lt:]  # second energy grid is the last LT points
        fx = np.array(data[1::2][:-lt])
        fy = np.array(data[1::2][-lt:])
        eg = union_grid(xx, yy)
        cov = fy.reshape(-1,1) * fy
        cov = cls(xx, yy, cov).reinit(eg)
        for i,e1 in enumerate(eg):
            for j,e2 in enumerate(eg):
                if j < i:
                    continue
                f = 0
                for k,e in enumerate(xx[:-1]):
                    if e1 >= e and e1 < xx[k+1] and e2 >= e and e2 < xx[k+1]:
                        f = fx[k]
                cov[i,j] = cov[j,i] = cov[i,j] * f
        return cls(xx, yy, cov, lb)

    @classmethod
    def lb5(cls, endf):
        r"""
        Classmethod to initialize the covariance from the text read in `LB=5`
        format of `ENDF-6`.
        ..`ENDF-6` description::
            Relative covariance matrix of some cross sections averaged over
            some energy intervals.
            LS=0 Asymmetric matrix, LS=1 Symmetric matrix.

        ..math:
            rcov(x,y) = g(Ex,Ey)

        Inputs:
            - endf :
                *endf.ENDF* instance

        Outputs:
            - inst :
                *Cov33* instance

        Section parameters:
            - ne :
                number of entries in the array {:math:`E_k`} defining
                (NE-1) energy intervals.
            - ls :
                flag stating if the matrix is asymmetric (`ls=0`) or
                symmetric (`ls=1`).
        """
        from sandy.cov import triu_matrix
        _, _, ls, lb, nt, ne, data = endf.read_list()
        xx = yy = data[:ne]
        f = np.array(data[ne:])
        if ls == 0: # to be tested
            cov = f.reshape(ne-1, ne-1)
        else:
            cov = triu_matrix(f, ne-1)
        # add zero row and column at the end of the matrix
        cov = np.insert(cov, cov.shape[0], [0]*cov.shape[1], axis=0)
        cov = np.insert(cov, cov.shape[1], [0]*cov.shape[0], axis=1)
        return cls(xx, yy, cov, lb)

    @classmethod
    def lb6(cls, endf):
        r"""
        Classmethod to initialize the covariance from the text read in ``LB=6``
        format of ``ENDF-6``.
        A relative covariance matrix interrelating the cross sections for
        two different reaction types or materials with (generally)
        different energy grids for its rows and columns.
        :math:`rcov(x,y) = g'(Ex,Ey)`

        Inputs

        Outputs

        ENDF-6 description
        ==================
        A relative covariance matrix interrelating the cross sections for
        two different reaction types or materials with (generally) different
        energy grids for its rows and columns.
        :math:`rcov(x,y) = g'(Ex,Ey)`
        """
        from sandy.functions import union_grid
        _, _, lt, lb, nt, ner, data = endf.read_list()
        nt = len(data)
        nec = (nt - 1)//ner
        er = data[:ner]
        ec = data[ner:ner+nec]
        f = np.array(data[ner+nec:])
        cov = f.reshape(ner-1, nec-1)
        # add zero row and column at the end of the matrix
        cov = np.insert(cov, cov.shape[0], [0]*cov.shape[1], axis=0)
        cov = np.insert(cov, cov.shape[1], [0]*cov.shape[0], axis=1)
        ug = union_grid(er, ec)
        inst = cls(er, ec, cov, lb)
        inst = inst.reinit(ug)
        return inst

    @classmethod
    def lb8(cls, endf):
        r"""
        Classmethod to initialize the covariance from the text read in `LB=8`
        format of `ENDF-6`.
        ..`ENDF-6` description::
            In general, each :math:`f_k` characterizes an uncorrelated
            contribution to the absolute variance of the indicated cross
            section averaged over any energy interval (subgroup) that includes
            a portion of the energy interval :math:`\Delta E_{k}`.
            The variance contribution :math:`Var(X_{jj})` from an `LB=8`
            sub-subsection to the processed group variance for the energy
            group :math:`(E_j , E_{j+1} ) is inversely proportional to its
            width :math:`\Delta E_j` when :math`(E_j , E_{ j+1} )` lies
            within :math:`(E_k , E_{k+1} )` and is obtained from the relation:
        ..math::
            Var(X_{jj} ) = f_k \Delta E_k/\Delta E_j
        """
        inst = cls.lb1(endf)
        inst.lb = 8
        return inst

    def __new__(cls, xx, yy, cov, *args, **kwargs):
        obj = Cov.__new__(cls, cov)
        return obj

    def __init__(self, xx, yy, matrix, lb=100):
        self.xx = xx
        self.yy = yy
        self.lb = lb

    @property
    def lb(self):
        r"""
        Flag whose numerical value determines how the covariance matrix is given.
        """
        return self._lb

    @lb.setter
    def lb(self, lb):
        if lb not in [0, 1, 2, 3, 4, 5, 6, 7, 8, 100]:
            raise NotImplementedError("Invalid 'lb' value for 'Cov33' object")
        self._lb = lb

    @property
    def xx(self):
        r"""
        Array of tabulated x-values.
        """
        return self._xx

    @xx.setter
    def xx(self, xx):
        from sandy.records import Grid
        self._xx = Grid(xx, size=self.shape[0])

    @property
    def yy(self):
        r"""
        Array of tabulated y-values.
        """
        return self._yy

    @yy.setter
    def yy(self, yy):
        from sandy.records import Grid
        self._yy = Grid(yy, size=self.shape[1])

    def corr2cov(self, s):
        dim = self.shape[0]
        S = np.repeat(s, dim).reshape(dim, dim)
        cov = self.__class__(self.xx, self.yy, S.T * (self * S), lb=self.lb)
        cov = cov.up2down()
        return cov

    def up2down(self):
        r"""
        Copy the upper triangular part to the lower triangular part.

        Outputs:
            - :``cov``: :
                (2d-array) output covariance matrix
        """
        U = np.triu(self)
        L = np.triu(self, 1).T
        cov = self.__class__(self.xx, self.yy, U + L, lb=self.lb)
        return cov

    def reinit(self, xx, yy=None):
        r"""
        Reshape covariance using zero interpolation.

        Inputs:
            - xx :
                new array of tabulated energy for 1st covariance dimension.
            - yy :
                new array of tabulated energy for 2nd covariance dimension.
                If not given it is set equal to *xx*.

        Outputs :
            - inst :
                new *Cov33* instance with arguments *xx*, *yy* and *cov*,
                where *cov* is the reshaped covariance matrix.
        """
        from sandy.functions import zero_interp
        xx = np.atleast_1d(xx)
        if yy is None:
            yy = xx
        yy = np.atleast_1d(yy)
        cov = zero_interp(self.xx, self, xx)
        cov = zero_interp(self.yy, cov.T, yy)
        cov = cov.T
        inst = self.__class__(xx, yy, cov, self.lb)
        return inst

    @property
    def Var(self):
        """
        Extract the variance array.

        Outputs:
            - var :
                (1d array) variance array
        """
        from sandy.records import Tab1
        var = np.diag(np.array(self))
        if (var < 0).any():
            raise ValueError("Variances must be non-negative")
        return Tab1(self.xx, var)

    @property
    def Std(self):
        """
        Extract the standard deviation array.

        Outputs:
            - std :
                (1d array) standard deviation array
        """
        from sandy.records import Tab1
        std = np.sqrt(np.array(self.Var))
        return Tab1(self.xx, std)

    def write_to_csv(self, file, c1='Energy (eV)', c2='Stdev (%)'):
        """
        Write standard deviation and correlation matrix to `csv` file.

        Inputs:
            - file :
                (string) csv file name
        """
        import csv
        todump = [[c1, c2] + self.xx.tolist()]
        for y,s,row in zip(self.yy, self.std*100., self.corr*1000.):
            todump.append([y] + [s] + row.tolist())
        with open(file, 'w') as fout:
            writer = csv.writer(fout, quoting=csv.QUOTE_ALL)
            writer.writerows(todump)

    def plot(self):
        r"""
        Plot covariance matrix as a pseudocolor plot of a 2-D array.
        The colorbar is also added to the figure.
        """
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        X, Y = np.meshgrid(self.xx, self.yy)
        pcm = ax.pcolormesh(X, Y, self.corr, vmin=-1, vmax=1, cmap='bwr')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_aspect(1) # Height is 0.5 times the width
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
        fig.show()

class Ucov(dict):
    r"""

    *data* must be structured as a dictionary::

        data = { (key1,key2) : Cov33, ... }

    Hence, each dictionary keyword must be a tuple with two elements.
    Furthermore, each key can also be a tuple or else (not ndarray),
    as long as *key1* and *key2* remain symmetric.

    Extract pointers for the union covariance matrix.

    Example,

    |              l=0   l=0   l=0   l=1   l=1
    |              e=1   e=2   e=3   e=2   e=3
    | l=0 e=1       x     x     x     x     x
    | l=0 e=2       x     x     x     x     x
    | l=0 e=3       x     x     x     x     x
    | l=1 e=2       x     x     x     x     x
    | l=1 e=3       x     x     x     x     x

    Then the pointers are `[0,0,0,1,1]`.

    Inputs:
        - mt1 :
            (int) "other" reaction type

    Ouputs:
        - pointers :
            1d array of pointers
    """
    def __setitem__(self, key, item):
        """
        Filter the items to set in a ``cov.Ucov`` dictionary according to
        the following rules:
            * raise error if key is not a tuple of two elements
            * raise error if item is not of type ``mf33.Cov33``
        """
        from sandy.mf33 import Cov33
        if len(key) != 2:
            raise NotImplementedError("'Ucov' keys must be a tuple with size=2")
        if isinstance(item, Cov33):
            raise NotImplementedError("Item for 'key={}' must be of type 'Cov33'".format(key))
        super().__setitem__(key, item)

    def get_ucov(self):
        from collections import namedtuple
        Pointer = namedtuple('Pointer', ('beg', 'end'))
        indices = []
        for key in sorted(self):
            i,j = key
            if i != j:
                continue
            indices += [i]*len(self[i,j].xx)
        indices = np.array(indices)
        cov = Cov(np.zeros((len(indices),)*2))
        pntrs = {}
        for i,j in sorted(self):
            # this is for tuple indices
            mask = np.array([(ind == i).all() for ind in indices])
            beg = (mask).nonzero()[0][0]
            end = (mask).nonzero()[0][-1]
            pntrs[i] = Pointer(beg, end)
        for i,j in sorted(self):
            xi = self[i,i].xx
            xj = self[j,j].xx
            C = self[i,j].reinit(xi, xj)
            cov[pntrs[i].beg:pntrs[i].end+1, pntrs[j].beg:pntrs[j].end+1] = C
            cov[pntrs[j].beg:pntrs[j].end+1, pntrs[i].beg:pntrs[i].end+1] = C.T
        return pntrs, cov

    def sampling(self, nsmp):
        from sandy.sandy_input import options
        from sandy.records import Samples
        pntrs, cov = self.get_ucov()
        if "no_correlations" in options:
            if options['no_correlations']:
                cov = cov.empty_off_diagonals()
        samples = cov.sampling(nsmp)
        DICT = {}
        for k,p in pntrs.items():
            xx = self[k,k].xx
            DICT[k] = Samples(xx, samples[p.beg:p.end+1])
        return DICT

    def plot(self):
        """
        Plot the correlation matrix of the union covariance as a pseudocolor
        plot of a 2-D array.
        The colorbar is also added in the figure.
        """
        import matplotlib.pyplot as plt
        pntrs, cov = self.get_ucov()
        fig, ax = plt.subplots()
        corr = cov.corr()
        pcm = ax.matshow(corr, vmin=-1, vmax=1, cmap='bwr', aspect='auto')

        # now, add spacers to separate the different keys
        xtick_name = []
        ytick_name = []
        tick_pos = []

        for mt,p in sorted(pntrs.items()):
            ax.axhline(y=p.beg-0.5, linewidth=1, color='k')
            ax.axvline(x=p.beg-0.5, linewidth=1, color='k')
            delta = p.end - p.beg
            tick_pos.append(p.beg + delta/2)
            ytick_name.append("{}".format(mt))
            xtick_name.append("{}".format(mt))

        plt.xticks(tick_pos, xtick_name)
        plt.yticks(tick_pos, xtick_name)
        cbar = plt.colorbar(pcm)
        cbar.ax.tick_params()
        fig.show()



